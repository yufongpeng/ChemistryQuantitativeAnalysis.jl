"""
    read_calibration(file::String, T, caltable::AbstractAnalysisTable; delim = "\\t")

Read ".mcal" or ".scal" file into julia as `MultipleCalibration` or `SingleCalibration`. 

`T` is sink function for tabular data. `caltable` is the source calibration data. `delim` specifies delimiter for tabular data.

See README.md for the structure of ".mcal" and ".scal" file.
"""
function read_calibration(file::String, T, caltable::AbstractAnalysisTable; delim = "\t")
    endswith(file, ".mcal") || endswith(file, ".scal") || throw(ArgumentError("The file is not a valid calibration directory"))
    config_kw, config_vl = eachline(joinpath(file, "config.txt"))
    config = Dict{Symbol, Any}()
    for (k, v) in zip(split(config_kw, "\t"), split(config_vl, "\t"))
        if k == "type" || k == "zero"
            v = v == "TRUE" || v == "True" || v == "true" || v == ""
        elseif k == "weight"
            v = parse(Float64, v)
        elseif k == "delim"
            v = v == "\\t" ? "\t" : ","
        elseif k == "analyte" || k == "isd"
            v = parse(Int, v)
        else
            continue
        end
        get!(config, Symbol(k), v)
    end
    if endswith(file, ".scal")
        return SingleCalibration(config[:analyte], config[:isd], caltable)
    end
    tbl = CSV.read(joinpath(file, "table.txt"), T; delim = get(config, :delim, delim))
    calibration(tbl, caltable; config...)
end

"""
    read_analysistable(file::String, T; analyte_fn = default_analyte_fn, delim = "\\t")

Read ".tbl" file into julia as `ColumnAnalysisTable` or `RowAnalysisTable`. 

`T` is sink function for tabular data. `analyte_fn` is responsible for converting `analyte` in string type into user defined type. `delim` specifies delimiter for tabular data.

See README.md for the structure of ".tbl" file.
"""
function read_analysistable(file::String, T; analyte_fn = default_analyte_fn, delim = "\t")
    endswith(file, ".tbl") || throw(ArgumentError("The file is not a valid table directory"))
    config = Dict{Symbol, Any}()
    lasttag = nothing
    for i in eachline(joinpath(file, "config.txt"))
        isempty(i) && continue
        i == "\t" && continue
        tag = match(r"\[(.*)\]", i)
        if isnothing(tag)
            isnothing(lasttag) && continue
            v = get!(config, lasttag, i)
            if v != i
                config[lasttag] = push!(vectorize(v), i)
            end
        else
            lasttag = Symbol(first(tag))
        end
    end
    delim = get(config, :delim, delim)
    delim = delim == "\\t" ? "\t" : delim  
    if haskey(config, :ISD)
        isd_map = Symbol(config[:ISD])
        analyte_name = Symbol(config[:Analyte])
        sample_name = Symbol.(filter!(!isempty, vectorize(config[:Sample])))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => Float64), types = Dict(isd_map => Int, analyte_name => String), validate = false)
        isd_map = isd_map in propertynames(tbl) ? isd_map : nothing
        isnothing(isd_map) || replace!(getproperty(tbl, isd_map), missing => 0)
        for i in sample_name
            replace!(getproperty(tbl, i), missing => 0)
        end
        RowAnalysisTable(sample_name, analyte_name, analyte_fn.(getproperty(tbl, analyte_name)), isd_map, tbl)
    else
        sample_name = Symbol(first(split(config[:Sample], "\t")))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => Float64), types = Dict(sample_name => String))
        analyte_name = String[]
        isd_map = Int[]
        for i in config[:Analyte]
            s = split(i, "\t")
            push!(analyte_name, first(s))
            push!(isd_map, (length(s) == 1 || isempty(last(s))) ? 0 : parse(Int, last(s)))
        end
        getproperty(tbl, Symbol(first(analyte_name))) isa AbstractVector{<: Number} || throw(ArgumentError("StackAnalysisTable is not implemented yet."))
        for i in Symbol.(analyte_name)
            replace!(getproperty(tbl, i), missing => 0)
        end
        ColumnAnalysisTable(sample_name, Symbol.(analyte_name), analyte_fn.(analyte_name), isd_map, tbl)
    end
end

"""
    read_sampletable(file::String, T; analyte_fn = default_analyte_fn, delim = "\\t")

Read ".tbl" file into julia as `SampleWrapperTable`. 

`T` is sink function for tabular data. `analyte_fn` is responsible for converting `analyte` in string type into user defined type. `delim` specifies delimiter for tabular data.

See README.md for the structure of ".tbl" file.
"""
function read_sampletable(file::String, T; analyte_fn = default_analyte_fn, delim = "\t")
    endswith(file, ".tbl") || throw(ArgumentError("The file is not a valid table directory"))
    analysistable = read_analysistable(file, T; analyte_fn, delim)
    if "cal_map.txt" in readdir(file)
        cal_map = filter!(!isempty, readlines(joinpath(file, "cal_map.txt")))
        cal_map = [isempty(v) ? i : parse(Int, v) for (i, v) in enumerate(cal_map)]
    else
        cal_map = collect(eachindex(analysistable.analytes))
    end
    for i in findall(x -> in(x, analysistable.isds), analysistable.analytes)
        cal_map[i] = 0
    end
    SampleWrapperTable(cal_map, analysistable)
end

"""
    read_caltable(file::String, T; analyte_fn = default_analyte_fn, delim = "\\t")

Read ".ctbl" file into julia as `CalWrapperTable`. 

`T` is sink function for tabular data. `analyte_fn` is responsible for converting `analyte` in string type into user defined type. `delim` specifies delimiter for tabular data.

See README.md for the structure of ".ctbl" file.
"""
function read_caltable(file::String, T; analyte_fn = default_analyte_fn, delim = "\t")
    endswith(file, ".ctbl") || throw(ArgumentError("The file is not a valid table directory"))
    conc = read_analysistable(joinpath(file, "conc.tbl"), T; analyte_fn, delim)
    if length(conc.samples) > 1
        signal = read_analysistable(joinpath(file, "signal.tbl"), T; analyte_fn, delim)
        level_map = parse.(Int, filter!(!isempty, readlines(joinpath(file, "level_map.txt"))))
    else
        signal = nothing
        level_map = Int[]
    end
    CalWrapperTable(level_map, conc, signal)
end

"""
    read_project.read(file::String, T; analyte_fn = default_analyte_fn)

Read ".pjc" file into julia as `Project`. `T` is sink function for tabular data, `analyte_fn` is responsible for converting `analyte` in string type into user defined type.

See README.md for the structure of ".pjc" file.
"""
function read_project(file::String, T; analyte_fn = default_analyte_fn)
    endswith(file, ".pjc") || throw(ArgumentError("The file is not a valid project directory"))
    config = Dict{Symbol, Any}()
    lasttag = nothing
    for i in eachline(joinpath(file, "config.txt"))
        isempty(i) && continue
        i == "\t" && continue
        tag = match(r"\[(.*)\]", i)
        if isnothing(tag)
            isnothing(lasttag) && continue
            v = get!(config, lasttag, i)
            if v != i
                config[lasttag] = push!(vectorize(v), i)
            end
        else
            lasttag = Symbol(first(tag))
        end
    end 
    delim = config[:delim] == "\\t" ? "\t" : ","    
    caltable = read_caltable(joinpath(file, "cal.ctbl"), T; analyte_fn, delim)
    if !in("calibration", readdir(file)) || isempty(readdir(joinpath(file, "calibration")))
        if isnothing(caltable.signaltable)
            ids = [(i, find_isd(caltable.conctable, analyte)) for (i, analyte) in enumerate(caltable.conctable.analytes)]
            cal = [SingleCalibration(i, j, caltable) for (i, j) in ids if j > -1]
        else
            cal = [calibration(caltable, analyte) for analyte in caltable.conctable.analytes if find_isd(caltable.conctable, analyte) > -1]
        end
    else
        cal = [read_calibration(joinpath(file, "calibration", f), T, caltable; delim) for f in readdir(joinpath(file, "calibration")) if endswith(f, ".mcal") || endswith(f, ".scal")]
    end
    fs = findfirst(==("sample.tbl"), readdir(file))
    fr = findfirst(==("result.tbl"), readdir(file))
    Project(cal, caltable,
        isnothing(fs) ? nothing : read_sampletable(joinpath(file, "sample.tbl"), T; analyte_fn, delim),
        isnothing(fr) ? nothing : read_sampletable(joinpath(file, "result.tbl"), T; analyte_fn, delim)
    )
end

function show(io::IO, cal::MultipleCalibration{A, T}) where {A, T}
    print(io, "Calibration{$A, $(shorten_type_repr(T))} of $(cal.caltable.conctable.analytes[cal.analyte]) with $(length(unique(cal.table.level[cal.table.include]))) levels and $(length(findall(cal.table.include))) points")
end

function show(io::IO, cal::SingleCalibration{A, T}) where {A, T}
    print(io, "Calibration{$A, $(shorten_type_repr(T))} of $(cal.caltable.conctable.analytes[cal.analyte]) with single level")
end

function show(io::IO, ::MIME"text/plain", cal::MultipleCalibration)
    print(io, cal, ":\n")
    print(io, "∘ Analyte: ", cal.caltable.conctable.analytes[cal.analyte], "\n")
    print(io, "∘ ISD: ", cal.isd > 0 ? cal.caltable.conctable.analytes[cal.isd] : nothing, "\n")
    print(io, "∘ Type: ", cal.type ? "linear\n" : "quadratic\n")
    print(io, "∘ (0, 0): ", cal.zero ? "included\n" : "ommitted\n")
    print(io, "∘ Weight: ", weight_repr(cal.weight), "\n")
    print(io, "∘ Formula: ", formula_repr(cal), "\n")
    print(io, "∘ R²: ", r2(cal.model))
end

function show(io::IO, ::MIME"text/plain", cal::SingleCalibration)
    print(io, cal, ":\n")
    print(io, "∘ Analyte: ", cal.caltable.conctable.analytes[cal.analyte], "\n")
    print(io, "∘ ISD: ", cal.isd > 0 ? cal.caltable.conctable.analytes[cal.isd] : nothing, "\n")
    print(io, "∘ Concentration: ", first(get_analyte(cal.caltable.conctable, cal.isd)))
end

function show(io::IO, project::Project{A, T}) where {A, T}
    print(io, string("Project{$A, $(shorten_type_repr(T))} with $(length(project.caltable.conctable.isds)) internal standards out of $(length(project.caltable.conctable.analytes)) analytes"))
end

function show(io::IO, ::MIME"text/plain", project::Project)
    print(io, project, ":\n")
    print(io, "∘ Calibration: \n")
    show(io, MIME"text/plain"(), project.calibration)
    print(io, "\n\n∘ Sample: \n")
    show(io, MIME"text/plain"(), project.sampletable)
    print(io, "\n\n∘ Result: \n")
    show(io, MIME"text/plain"(), project.resulttable)
end

function shorten_type_repr(T)
    t = repr(T)
    length(t) > 50 ? replace(t, r"\{.*\}" => "{...}") : t
end

function show(io::IO, tbl::ColumnAnalysisTable{A, T}) where {A, T}
    print(io, string("ColumnAnalysisTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analytes)) analytes and $(length(tbl.samples)) samples"))
end

function show(io::IO, ::MIME"text/plain", tbl::ColumnAnalysisTable)
    print(io, tbl, ":\n")
    show(io, MIME"text/plain"(), tbl.table)
end

function show(io::IO, tbl::RowAnalysisTable{A, T}) where {A, T}
    print(io, string("RowAnalysisTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analytes)) analytes and $(length(tbl.samples)) samples"))
end

function show(io::IO, ::MIME"text/plain", tbl::RowAnalysisTable)
    print(io, tbl, ":\n")
    show(io, MIME"text/plain"(), tbl.table)
end

function show(io::IO, tbl::SampleWrapperTable{A, T}) where {A, T}
    print(io, string("SampleWrapperTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analytes)) analytes and $(length(tbl.samples)) samples"))
end

function show(io::IO, ::MIME"text/plain", tbl::SampleWrapperTable)
    c = [i > 0 ? tbl.analytes[i] : nothing for i in tbl.cal_map]
    if length(c) > 15
        c = string(join(first(c, 5), ", "), ", ..., ", join(last(c, 5), ", "))
    else
        c = join(c, ", ")
    end
    print(io, tbl, ":\n")
    print(io, "∘ Calibration curve: ", c, "\n")
    show(io, MIME"text/plain"(), tbl.table)
end

function show(io::IO, tbl::CalWrapperTable{A, T}) where {A, T}
    isnothing(tbl.signaltable) && return print(io, string("CalWrapperTable{$A, $(shorten_type_repr(T))} with $(length(tbl.conctable.analytes)) analytes"))
    print(io, string("CalWrapperTable{$A, $(shorten_type_repr(T))} with $(length(tbl.conctable.analytes)) analytes, $(length(tbl.conctable.samples)) levels and $(length(tbl.signaltable.samples)) points."))
end

function show(io::IO, ::MIME"text/plain", tbl::CalWrapperTable)
    print(io, tbl, ":\n")
    print(io, "∘ Level: ", join(tbl.level_map, ", "), "\n")
    print(io, "∘ Concentration: \n")
    show(io, MIME"text/plain"(), tbl.conctable.table)
    print(io, "\n∘ Signal: \n")
    show(io, MIME"text/plain"(), isnothing(tbl.signaltable) ? nothing : tbl.signaltable.table)
end

function write(file::String, tbl::RowAnalysisTable; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[Analyte]\n", tbl.analyte_name, "\n\n[ISD]\n", string(tbl.isd_map), "\n\n[Sample]\n", join(tbl.sample_name, "\n"))
    end
    CSV.write(joinpath(file, "table.txt"), tbl.table; delim)
end

function write(file::String, tbl::ColumnAnalysisTable; delim = "\t")
    mkpath(file)
    config = open(joinpath(file, "config.txt"), "w+")
    Base.write(config, "[Analyte]\n")
    for (i, j) in zip(tbl.analyte_name, tbl.isd_map)
        Base.write(config, i, "\t", string(j), "\n")
    end
    Base.write(config, "\n\n[Sample]\n", tbl.sample_name)
    close(config)
    CSV.write(joinpath(file, "table.txt"), tbl.table; delim)
end

function write(file::String, tbl::SampleWrapperTable; delim = "\t")
    write(file, tbl.analysistable; delim)
    open(joinpath(file, "cal_map.txt"), "w+") do cal_map
        Base.write(cal_map, join(tbl.cal_map, "\n"))
    end
end

function write(file::String, tbl::CalWrapperTable; delim = "\t")
    write(joinpath(file, "conc.tbl"), tbl.conctable; delim)
    isnothing(tbl.signaltable) || write(joinpath(file, "signal.tbl"), tbl.signaltable; delim)
    isempty(tbl.level_map) || open(joinpath(file, "level_map.txt"), "w+") do level_map
        Base.write(level_map, join(tbl.level_map, "\n"))
    end
end

function write(file::String, cal::MultipleCalibration; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "analyte\tisd\ttype\tzero\tweight\n", string(cal.analyte), "\t", string(cal.isd), "\t", string(cal.type), "\t", string(cal.zero), "\t", string(cal.weight))
    end
    CSV.write(joinpath(file, "table.txt"), cal.table; delim)
end
function write(file::String, cal::SingleCalibration; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "analyte\tisd\n", string(cal.analyte), "\t", string(cal.isd))
    end
end

"""
    Calibration.write(file::String, object; delim = "\\t")

Write `object` into ".scal" for `SingleCalibration`, ".mcal" for `MultipleCalibration`, ".tbl" for `AbstractAnalysisTable` and `SampleWrapperTable`, ".ctbl" for `CalWrapperTable`, and ".pjc" for `Project`. 

`delim` specifies delimiter for tabular data.

See README.md for the structure of ".pjc" file.
"""
function write(file::String, project::Project; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[delim]\n", delim == "\t" ? "\\t" : ",")
    end
    write(joinpath(file, "cal.ctbl"), project.caltable; delim)
    mkpath(joinpath(file, "calibration"))
    ft = isnothing(project.caltable.signaltable) ? "scal" : "mcal"
    for (i, c) in enumerate(project.calibration)
        write(joinpath(file, "calibration", "$i.$ft"), c; delim)
    end
    isnothing(project.sampletable) || write(joinpath(file, "sample.tbl"), project.sampletable; delim)
    isnothing(project.resulttable) || write(joinpath(file, "result.tbl"), project.resulttable; delim)
end

"""
    Calibration.read(file::String, T; analyte_fn = default_analyte_fn)

Read ".pjc" file into julia as `Project`. `T` is sink function for tabular data, `analyte_fn` is responsible for converting `analyte` in string type into user defined type.

See README.md for the structure of ".pjc" file.
"""
const read = read_project
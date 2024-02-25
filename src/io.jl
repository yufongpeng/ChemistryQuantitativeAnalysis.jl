function read_config(file::String)
    config = Dictionary{Symbol, Any}()
    lasttag = nothing
    for i in eachline(file)
        isempty(i) && continue
        i == "\t" && continue
        tag = match(r"\[(.*)\]", i)
        if isnothing(tag)
            isnothing(lasttag) && continue
            if haskey(config, lasttag)
                config[lasttag] = push!(vectorize(config[lasttag]), i)
            else
                set!(config, lasttag, i)
            end
        else
            lasttag = Symbol(first(tag))
        end
    end
    if haskey(config, :delim)
        config[:delim] = unescape_string(config[:delim])
    end
    config
end  

"""
    read_calibration(file::String; analytetype = String, numbertype = Float64, delim = '\\t') -> AbstractCalibration{analytetype, numbertype}

Read ".mcal" or ".scal" file into julia as `MultipleCalibration` or `SingleCalibration`. `analytetype` is a concrete type for `analyte`, 
, `numbertype` is the type for numeric data, and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist. 

`numbertype` is forced to `Float64` because of issue of `GLM`.
    
See `ChemistryQuantitativeAnalysis.read` for the requirement of `analytetype`.

See README.md for the structure of ".mcal" and ".scal" file.
"""
function read_calibration(file::String; analytetype = String, numbertype = Float64, delim = '\t')
    endswith(file, ".mcal") || endswith(file, ".scal") || throw(ArgumentError("The file is not a valid calibration directory"))
    numbertype = Float64 # Bug of GLM
    config = read_config(joinpath(file, "config.txt"))
    if endswith(file, ".scal")
        return SingleCalibration((analytetype(config[:analyte]), ), parse(numbertype, config[:conc]))
    end
    tbl = CSV.read(joinpath(file, "table.txt"), Table; delim = get(config, :delim, delim), typemap = Dict(Int => numbertype, Float64 => numbertype), types = Dict(:level => Int))
    calibration(convert.(analytetype, (config[:analyte], config[:isd])), tbl; type = parse(Bool, config[:type]), zero = parse(Bool, config[:zero]), weight = parse(Float64, config[:weight]))
end

"""
    read_datatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t', levelname = :level) -> AbstractDataTable{analytetype, sampletype, numbertype}

Read ".dt" file into julia as `ColumnDataTable` or `RowDataTable`. `T` is the sink function for tabular data, 
`analytetype` is a concrete type for `analyte`, `sampletype` is a concrete type for `sample`, `numbertype` is the type for numeric data, 
and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist. `level` is specifically used for AnalysisMethod, indicating the column representing calibration level; this column should be all integers. 

See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".dt" file.
"""
function read_datatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t', levelname = nothing)
    endswith(file, ".dt") || throw(ArgumentError("The file is not a valid table directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    if config[:Type] == "R"
        analyte_name = Symbol(config[:Analyte])
        sample_name = String.(filter!(!isempty, vectorize(config[:Sample])))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => numbertype), types = Dict(analyte_name => analytetype), validate = false)
        for i in Symbol.(sample_name)
            replace!(getproperty(tbl, i), missing => 0)
        end
        RowDataTable(getproperty(tbl, analyte_name), cqamap(sampletype, sample_name), analyte_name, tbl)
    else
        config[:Type] == "C" || throw(ArgumentError("Other datatables are not implemented yet."))
        sample_name = Symbol(first(split(config[:Sample], "\t")))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => numbertype), types = Dict(sample_name => sampletype, levelname => Int), validate = false)
        analyte_name = String.(filter!(!isempty, vectorize(config[:Analyte])))
        for i in Symbol.(analyte_name)
            replace!(getproperty(tbl, i), missing => 0)
        end
        ColumnDataTable(cqamap(analytetype, analyte_name), getproperty(tbl, sample_name), sample_name, tbl)
    end
end

"""
    read_analysistable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t') -> AnalysisTable{analytetype, sampletype}

Read ".at" file into julia as `AnalysisTable`. `T` is the sink function for tabular data, 
`analytetype` is a concrete type for `analyte`, `sampletype` is a concrete type for `sample`, `numbertype` is the type for numeric data, 
and `delim` specifies delimiter for tabular data if `config[:delim]` in `tables` do not exist. 
    
See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".at" file.
"""
function read_analysistable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t')
    endswith(file, ".at") || throw(ArgumentError("The file is not a valid table directory"))
    files = filter!(f -> endswith(f, ".dt"), readdir(file))
    tables = map(files) do f
        read_datatable(joinpath(file, f), T; analytetype, sampletype, numbertype, delim)
    end
    AnalysisTable(Symbol.(replace.(files, Ref(".dt" => ""), Ref(r"^\d*_" => ""))), tables)
end
"""
    read_method(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t') -> AnalysisMethod{analytetype, <: Table}

Read ".mt" file into julia as `AnalysisMethod`. `T` is the sink function for tabular data, 
`analytetype` is a concrete type for `analyte`, `sampletype` is a concrete type for `sample`, `numbertype` is the type for numeric data, 
and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist. 
    
See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".mt" file.
"""
function read_method(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t')
    endswith(file, ".mt") || throw(ArgumentError("The file is not a valid table directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    signal = Symbol(get(config, :signal, :area))
    analytetable = CSV.read(joinpath(file, "analytetable.txt"), Table)
    analyte = analytetype.(analytetable.analyte)
    isd = replace(analytetable.isd, missing => 0)
    conctable = read_datatable(joinpath(file, "true_concentration.dt"), T; analytetype, sampletype = Int, numbertype, delim)
    if length(sampleobj(conctable)) > 1
        signaltable = read_datatable(joinpath(file, "$signal.dt"), T; analytetype, sampletype, numbertype, delim, levelname = get(config, :levelname, nothing))
        if haskey(config, :levelname) && in(Symbol(config[:levelname]), propertynames(signaltable))
            pointlevel = getproperty(signaltable, Symbol(config[:levelname]))
        elseif haskey(config, :pointlevel)
            pointlevel = parse.(Int, config[:pointlevel])
        else
            nl = length(sampleobj(conctable))
            ns = length(sampleobj(signaltable))
            nr = ns ÷ nl
            rr = ns - nl * nr
            pointlevel = vcat(repeat([1], rr), repeat(1:nl, nr))
        end
        calibration = map(enumerate(analytetable.calibration)) do (i, c)
            ismissing(c) ? i : c
        end
    else
        signaltable = nothing
        pointlevel = [1]
        calibration = isd
    end
    AnalysisMethod(Table(analytetable; analyte, isd, calibration), signal, pointlevel, conctable, signaltable)
end
"""
    read_batch(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t') -> Batch{analytetype}

Read ".batch" file into julia as `Batch`. `T` is the sink function for tabular data, 
`analytetype` is a concrete type for `analyte`, `sampletype` is a concrete type for `sample`, `numbertype` is the type for numeric data, 
and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist. 
    
See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".batch" file.
"""
function read_batch(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t')
    endswith(file, ".batch") || throw(ArgumentError("The file is not a valid batch directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)   
    method = read_method(joinpath(file, "method.mt"), T; analytetype, sampletype, numbertype, delim)
    if !in("calibration", readdir(file)) || isempty(readdir(joinpath(file, "calibration")))
        if isnothing(method.signaltable)
            cal = [SingleCalibration((analyte, ), first(getanalyte(method.conctable, analyte))) for analyte in analyteobj(method.conctable)]
        else
            cal = [calibration(method, analyte) for analyte in analyteobj(method.conctable) if !isisd(analyte, method)]
        end
    else
        cal = [read_calibration(joinpath(file, "calibration", f); delim, analytetype, numbertype) for f in readdir(joinpath(file, "calibration")) if endswith(f, ".mcal") || endswith(f, ".scal")]
    end
    fd = findfirst(==("data.at"), readdir(file))
    Batch(method, cal,
        isnothing(fd) ? nothing : read_analysistable(joinpath(file, "data.at"), T; analytetype, sampletype, numbertype, delim),
    )
end

function print_summary(io::IO, cal::MultipleCalibration{A}) where A
    print(io, "MultipleCalibration{$A} of ", first(cal.analyte), " with ", length(unique(cal.table.level[cal.table.include])), " levels and ", length(findall(cal.table.include)), " points")
end

function print_summary(io::IO, cal::SingleCalibration{A}) where A
    print(io, "SingleCalibration{$A} of ", first(cal.analyte), " with single level")
end

function show(io::IO, ::MIME"text/plain", cal::MultipleCalibration{A, N, T}) where {A, N, T}
    print(io, "MultipleCalibration{$A, $N, $(shorten_type_repr(T))} with ", length(unique(cal.table.level[cal.table.include])), " levels and ", length(findall(cal.table.include)), " points:\n")
    print(io, "∘ Analyte: ", first(cal.analyte), "\n")
    print(io, "∘ ISD: ", last(cal.analyte), "\n")
    print(io, "∘ Type: ", cal.type ? "linear\n" : "quadratic\n")
    print(io, "∘ (0, 0): ", cal.zero ? "included\n" : "ommitted\n")
    print(io, "∘ Weight: ", weight_repr(cal.weight), "\n")
    print(io, "∘ Formula: ", formula_repr(cal), "\n")
    print(io, "∘ R²: ", r2(cal.model))
end

function show(io::IO, ::MIME"text/plain", cal::SingleCalibration{A, N}) where {A, N}
    print(io, "SingleCalibration{$A, $N} with single level:\n")
    print(io, "∘ Analyte (ISD): ", first(cal.analyte), "\n")
    print(io, "∘ Concentration: ", cal.conc)
end

function print_summary(io::IO, batch::Batch{A, T}) where {A, T}
    print(io, "Batch{$A} with ", count(<(0), batch.method.analytetable.isd), " internal standards out of ", length(batch.method.analytetable.analyte), " analytes")
end

function show(io::IO, ::MIME"text/plain", batch::Batch)
    print_summary(io, batch)
    print(io, ":\n∘ Analytes:\n")
    show(io, MIME"text/plain"(), batch.method.analytetable)
    print(io, "\n\n∘ Calibration:")
    if length(batch.calibration) > 10
        for c in @view batch.calibration[1:5]
            print(io, "\n ")
            print_summary(io, c)
        end
        print(io, "\n ⋮")
        for c in @view batch.calibration[end - 4:end]
            print(io, "\n ")
            print_summary(io, c)
        end
    else
        for c in batch.calibration
            print(io, "\n ")
            print_summary(io, c)
        end
    end
    print(io, "\n\n∘ Data:\n")
    show(io, MIME"text/plain"(), batch.data)
end

function shorten_type_repr(T)
    t = repr(T)
    length(t) > 50 ? replace(t, r"\{.*\}" => "{…}") : t
end

function print_summary(io::IO, tbl::ColumnDataTable{A, S, N, T}) where {A, S, N, T}
    print(io, "ColumnDataTable{$A, $S, $N, $(shorten_type_repr(T))} with ", length(analyteobj(tbl)), " analytes and ", length(sampleobj(tbl)), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::ColumnDataTable)
    print_summary(io, tbl)
    println(io, ":")
    show(io, MIME"text/plain"(), table(tbl))
end

function print_summary(io::IO, tbl::RowDataTable{A, S, N, T}) where {A, S, N, T}
    print(io, "RowDataTable{$A, $S, $N, $(shorten_type_repr(T))} with ", length(analyteobj(tbl)), " analytes and ", length(sampleobj(tbl)), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::RowDataTable)
    print_summary(io, tbl)
    println(io, ":")
    show(io, MIME"text/plain"(), table(tbl))
end

function print_summary(io::IO, tbl::AnalysisTable{A, S, T}) where {A, S, T}
    print(io, "AnalysisTable{$A, $S, $(shorten_type_repr(T))} with ", length(analyteobj(tbl)), " analytes and ", length(sampleobj(tbl)), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::AnalysisTable)
    print_summary(io, tbl)
    print(io, ":")
    if length(tbl) > 2
        ks = propertynames(tbl)
        print(io, "\n∘ ", ks[1], " | ")
        show(io, MIME"text/plain"(), tbl[ks[1]])
        print(io, "\n   ⋮")
        print(io, "\n∘ ", ks[end], " | ")
        show(io, MIME"text/plain"(), tbl[ks[end]])
    else
        for (k, v) in pairs(tbl)
            print(io, "\n∘ ", k, " | ")
            show(io, MIME"text/plain"(), v)
        end
    end

end

function print_summary(io::IO, tbl::AnalysisMethod{A}) where A
    if isnothing(tbl.signaltable)
        print(io, "AnalysisMethod{$A} with ", count(<(0), tbl.analytetable.isd), " internal standards out of ", length(tbl.analytetable.analyte), " analytes")
    else
        print(io, "AnalysisMethod{$A} with ", count(<(0), tbl.analytetable.isd), " internal standards out of ", length(tbl.analytetable.analyte), " analytes, ", length(sampleobj(tbl.conctable)), " levels and ", length(sampleobj(tbl.signaltable)), " points")
    end
end

function show(io::IO, ::MIME"text/plain", tbl::AnalysisMethod)
    print_summary(io, tbl)
    print(io, ":\n∘ Analytes:\n")
    show(io, MIME"text/plain"(), tbl.analytetable)
    print(io, "\n\n∘ Level: ", join(tbl.pointlevel, ", "), "\n")
    print(io, "∘ Concentration: \n")
    show(io, MIME"text/plain"(), table(tbl.conctable))
    print(io, "\n\n∘ Signal: \n")
    show(io, MIME"text/plain"(), isnothing(tbl.signaltable) ? nothing : table(tbl.signaltable))
end

function write(file::String, tbl::RowDataTable; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[Type]\nR\n\n[delim]\n", escape_string(string(delim)), "\n\n[Analyte]\n", analytecol(tbl), "\n\n[Sample]\n", join(samplename(tbl), "\n"))
    end
    CSV.write(joinpath(file, "table.txt"), table(tbl); delim)
end

function write(file::String, tbl::ColumnDataTable; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[Type]\nC\n\n[delim]\n", escape_string(string(delim)), "\n\n[Analyte]\n", join(analytename(tbl), "\n"), "\n\n[Sample]\n", samplecol(tbl))
    end
    CSV.write(joinpath(file, "table.txt"), table(tbl); delim)
end

function write(file::String, tbl::AnalysisTable; delim = '\t')
    mkpath(file)
    rm.(readdir(file; join = true); recursive = true)
    for (i, (k, v)) in enumerate(pairs(tbl))
        write(joinpath(file, "$(i - 1)_$k.dt"), v; delim)
    end
end

function write(file::String, tbl::AnalysisMethod; delim = '\t')
    mkpath(file)
    write(joinpath(file, "true_concentration.dt"), tbl.conctable; delim)
    isnothing(tbl.signaltable) || write(joinpath(file, "$(tbl.signal).dt"), tbl.signaltable; delim)
    id = nothing
    if isa(tbl.signaltable, ColumnDataTable)
        id = findfirst(x -> getproperty(tbl.signaltable, x) == tbl.pointlevel, propertynames(tbl.signaltable))
    end
    open(joinpath(file, "config.txt"), "w+") do config
        isnothing(id) ? Base.write(config, "[signal]\n", tbl.signal, "\n\n[delim]\n", escape_string(string(delim)), "\n\n[pointlevel]\n", join(tbl.pointlevel, "\n")) : 
            Base.write(config, "[signal]\n", tbl.signal, "\n\n[delim]\n", escape_string(string(delim)), "\n\n[levelname]\n", propertynames(tbl.signaltable)[id])
    end
    CSV.write(joinpath(file, "analytetable.txt"), tbl.analytetable; delim)
end

function write(file::String, cal::MultipleCalibration; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[analyte]\n", string(first(cal.analyte)), "\n\n[isd]\n", string(last(cal.analyte)), 
                "\n\n[type]\n", string(cal.type), "\n\n[zero]\n", string(cal.zero), "\n\n[weight]\n", string(cal.weight), 
                "\n\n[delim]\n", escape_string(string(delim)))
    end
    CSV.write(joinpath(file, "table.txt"), cal.table; delim)
end
function write(file::String, cal::SingleCalibration; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[analyte]\n", string(first(cal.analyte)), "\n\n[conc]\n", string(cal.conc))
    end
end

"""
    ChemistryQuantitativeAnalysis.write(file::String, object; delim = "\\t")

Write `object` into ".scal" for `SingleCalibration`, ".mcal" for `MultipleCalibration`, ".at" for `AnalysisTable`, ".mt" for `AnalysisMethod`, and ".batch" for `Batch`. 

`delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

See README.md for the structure of all files.
"""
function write(file::String, batch::Batch; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[delim]\n", escape_string(string(delim)))
    end
    write(joinpath(file, "method.mt"), batch.method; delim)
    mkpath(joinpath(file, "calibration"))
    ft = isnothing(batch.method.signaltable) ? "scal" : "mcal"
    for (i, c) in enumerate(batch.calibration)
        write(joinpath(file, "calibration", "$i.$ft"), c; delim)
    end
    isnothing(batch.data) || write(joinpath(file, "data.at"), batch.data; delim)
end

"""
    ChemistryQuantitativeAnalysis.read(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t') -> Batch{analytetype}

Read ".scal" as `SingleCalibration`, ".mcal" as `MultipleCalibration`, ".at" as `AnalysisTable`, ".mt" as `AnalysisMethod`, and ".batch" as `Batch`. 

`T` is the sink function for tabular data, 
`analytetype` is a concrete type for `analyte`, `sampletype` is a concrete type for `sample`, `numbertype` is the type for numeric data, 
and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

For `analytetype` and `sampletype`, `string(cqaconvert(analytetype, x))` and `string(cqaconvert(sampletype, x))` should equal `x` if `x` is a valid string. Additionally, `tryparse` have to be extended for `CSV` parsing:
* `tryparse(::Type{analytetype}, x::AbstractString)` is neccessary for `RowDataTable`.
* `tryparse(::Type{sampletype}, x::AbstractString)` is neccessary for `ColumnDataTable`.

See README.md for the structure of ".batch" file.
"""
function read(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t')
    if endswith(file, ".batch")
        read_batch(file, T; analytetype, sampletype, numbertype, delim)
    elseif endswith(file, ".scal") || endswith(file, ".mcal")
        read_calibration(file; analytetype, numbertype, delim)
    elseif endswith(file, ".at")
        read_analysistable(file, T; analytetype, sampletype, numbertype, delim)
    elseif endswith(file, ".mt")
        read_method(file, T; analytetype, sampletype, numbertype, delim)
    elseif endswith(file, ".dt")
        read_datatable(file, T; analytetype, sampletype, numbertype, delim)
    end
end

"""
    mkbatch(file::String; 
            delim = '\t',
            data_config = Dict{Symbol, Any}(), 
            signal_config = Dict{Symbol, Any}(), 
            conc_config = Dict{Symbol, Any}(), 
            method_config = Dict{Symbol, Any}()
            )

Create a template directory of a batch. See "README.md" for available keys and values for config.

The default table wrapper is `ColumnDataTable`, which default sample column name is "Sample", and default sample column name for levels is "Level".
The default analyte column name for `RowDataTable` is "Analyte".

The priorities of default analytes are `conc_config[:Analyte]`, `signal_config[:Analyte]`, `data_config[:Analyte]`, and lastly `["Analyte1", "Analyte2"]`.
The default samples for data are `["S1", "S2"]`. The default calibration points are `["C1", "C2", "C3", "C4", "C5"]`. The default calibration levels are `[1, 2, 3, 4, 5]`.

"""
function mkbatch(file::String; 
                delim = '\t',
                data_config = Dict{Symbol, Any}(), 
                signal_config = Dict{Symbol, Any}(), 
                conc_config = Dict{Symbol, Any}(), 
                method_config = Dict{Symbol, Any}()
                )
    mkpath(file)
    default_analyte = get(conc_config, :Analyte, get(signal_config, :Analyte, get(data_config, :Analyte, ["Analyte1", "Analyte2"])))
    default_sample = ["S1", "S2"]
    default_point = ["C1", "C2", "C3", "C4", "C5"]
    default_level = [1, 2, 3, 4, 5]
    default_config_c = [
        Dict(:Type => "C", :delim => delim, :Analyte => default_analyte, :Sample => "Sample"), 
        Dict(:Type => "C", :delim => delim, :Analyte => default_analyte, :Sample => "Sample"), 
        Dict(:Type => "C", :delim => delim, :Analyte => default_analyte, :Sample => "Level")
    ]
    default_config_r = [
        Dict(:Type => "R", :delim => delim, :Analyte => "Analyte", :Sample => default_sample), 
        Dict(:Type => "R", :delim => delim, :Analyte => "Analyte", :Sample => default_point), 
        Dict(:Type => "R", :delim => delim, :Analyte => "Analyte", :Sample => default_level)
    ]
    table = Vector{Table}(undef, 3)
    confs = Vector{String}(undef, 3)
    data_config = convert(Dict{Symbol, Any}, data_config)
    signal_config = convert(Dict{Symbol, Any}, signal_config)
    conc_config = convert(Dict{Symbol, Any}, conc_config)
    method_config = convert(Dict{Symbol, Any}, method_config)
    for (i, (config, default_c, default_r)) in enumerate(zip([data_config, signal_config, conc_config], default_config_c, default_config_r))
        default = get(config, :Type, "C") == "C" ? default_c : default_r
        for (k, v) in default
            get!(config, k, v)
        end
        if config[:Type] == "C"
            #table[i] = string(config[:Sample], config[:delim], join(config[:Analyte], config[:delim]))
            config[:Analyte_merge] = join(config[:Analyte], "\n")
            config[:Sample_merge] = config[:Sample]
        else
            #table[i] = string(config[:Analyte], config[:delim], join(config[:Sample], config[:delim]))
            config[:Sample_merge] = join(config[:Sample], "\n")
            config[:Analyte_merge] = config[:Analyte]
        end
        confs[i] = string("[Type]\n", config[:Type], "\n\n[delim]\n", escape_string(string(config[:delim])), 
                                    "\n\n[Analyte]\n", config[:Analyte_merge], 
                                    "\n\n[Sample]\n", config[:Sample_merge])
    end
    if data_config[:Type] == "C"
        #table[1] = string(table[1], "\n", join((string(x, data_config[:delim], join(repeat(["1.0"], length(data_config[:Analyte])), data_config[:delim])) for x in default_sample), "\n"))
        table[1] = Table(; Symbol(data_config[:Sample]) => default_sample, (Symbol.(data_config[:Analyte]) .=> Ref(repeat([1.0], length(default_sample))))...)
    else
        table[1] = Table(; Symbol(data_config[:Analyte]) => default_analyte, (Symbol.(data_config[:Sample]) .=> Ref(repeat([1.0], length(default_analyte))))...)
    end
    if signal_config[:Type] == "R" && conc_config[:Type] == "R"
        default_level = parse.(Int, string.(conc_config[:Sample]))
        default_point = signal_config[:Sample]
        pointlevel = [i > lastindex(default_level) ? default_level[end] : default_level[i] for i in eachindex(default_point)]
    elseif signal_config[:Type] == "R"
        default_point = signal_config[:Sample]
        default_level = collect(eachindex(default_point))
        pointlevel = default_level
    elseif conc_config[:Type] == "R"
        default_level = parse.(Int, string.(conc_config[:Sample]))
        default_point = map(default_level) do s
            string("C", s)            
        end
        pointlevel = default_level
    end
    default_method = if signal_config[:Type] == "C"
        Dict(:signal => "area", :delim => delim, :levelname => "Level")
    else
        Dict(:signal => "area", :delim => delim, :pointlevel => pointlevel)
    end
    for (k, v) in default_method
        get!(method_config, k, v)
    end
    if haskey(method_config, :pointlevel)
        if length(method_config[:pointlevel]) != length(pointlevel) || any(x -> !in(x, default_level), method_config[:pointlevel])
            @warn "Replace poinlevel with default: $(pointlevel)"
            method_config[:pointlevel] = pointlevel
        end
    end
    if signal_config[:Type] == "C"
        confms = string("[signal]\n", method_config[:signal], "\n\n[delim]\n", escape_string(string(method_config[:delim])), 
                        "\n\n[levelname]\n", method_config[:levelname])
        table[2] = Table(; Symbol(signal_config[:Sample]) => default_point, Symbol(method_config[:levelname]) => default_level, (Symbol.(signal_config[:Analyte]) .=> Ref(repeat([1.0], length(default_point))))...)
    else
        confms = string("[signal]\n", method_config[:signal], "\n\n[delim]\n", escape_string(string(method_config[:delim])), 
                        "\n\n[pointlevel]\n", join(method_config[:pointlevel], "\n"))
        table[2] = Table(; Symbol(signal_config[:Analyte]) => default_analyte, (Symbol.(signal_config[:Sample]) .=> Ref(repeat([1.0], length(default_analyte))))...)
    end
    if conc_config[:Type] == "C"
        table[3] = Table(; Symbol(conc_config[:Sample]) => default_level, (Symbol.(conc_config[:Analyte]) .=> Ref(repeat([1.0], length(default_level))))...)
    else
        table[3] = Table(; Symbol(conc_config[:Analyte]) => default_analyte, (Symbol.(conc_config[:Sample]) .=> Ref(repeat([1.0], length(default_analyte))))...)
    end
    signal = method_config[:signal]
    data_name = string(0, "_", signal, ".dt")
    mkpath(joinpath(file, "data.at", data_name))
    write(joinpath(file, "data.at", data_name, "config.txt"), confs[1])
    CSV.write(joinpath(file, "data.at", data_name, "table.txt"), table[1]; delim = data_config[:delim])
    mkpath(joinpath(file, "method.mt", "true_concentration.dt"))
    mkpath(joinpath(file, "method.mt", "$signal.dt"))
    write(joinpath(file, "method.mt", "$signal.dt", "config.txt"), confs[2])
    CSV.write(joinpath(file, "method.mt", "$signal.dt", "table.txt"), table[2]; delim = signal_config[:delim])
    write(joinpath(file, "method.mt", "true_concentration.dt", "config.txt"), confs[3])
    CSV.write(joinpath(file, "method.mt", "true_concentration.dt", "table.txt"), table[3]; delim = conc_config[:delim])
    CSV.write(joinpath(file, "method.mt", "analytetable.txt"), Table(; analyte = default_analyte, isd = repeat([0], length(default_analyte)), calibration = collect(eachindex(default_analyte))); delim = method_config[:delim])
    write(joinpath(file, "method.mt", "config.txt"), confms)
    write(joinpath(file, "config.txt"), string("[delim]\n", escape_string(string(delim))))
end
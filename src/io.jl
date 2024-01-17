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
    read_calibration(file::String; analytetype::Type{A} = String, delim = '\\t') -> AbstractCalibration{A}

Read ".mcal" or ".scal" file into julia as `MultipleCalibration` or `SingleCalibration`. `analytetype` is a concrete type for `analyte` which msut have a method for string input, 
and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

See README.md for the structure of ".mcal" and ".scal" file.
"""
function read_calibration(file::String; analytetype = String, delim = '\t')
    endswith(file, ".mcal") || endswith(file, ".scal") || throw(ArgumentError("The file is not a valid calibration directory"))
    config = read_config(joinpath(file, "config.txt"))
    if endswith(file, ".scal")
        return SingleCalibration((analytetype(config[:analyte]), ), parse(Float64, config[:conc]))
    end
    tbl = CSV.read(joinpath(file, "table.txt"), Table; delim = get(config, :delim, delim))
    calibration(analytetype.((config[:analyte], config[:isd])), tbl; type = parse(Bool, config[:type]), zero = parse(Bool, config[:zero]), weight = parse(Float64, config[:weight]))
end

"""
    read_datatable(file::String, T; analytetype::Type{A} = String, delim = '\\t', levelname = :level) -> AbstractDataTable{A, S <: T}

Read ".dt" file into julia as `ColumnDataTable` or `RowDataTable`. `T` is the sink function for tabular data, `analytetype` is a concrete type for `analyte` which msut have a method for string input, 
and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist. `level` is specifically used for methodtable, indicating the column representing calibration level; this column should be all integers.

See README.md for the structure of ".dt" file.
"""
function read_datatable(file::String, T; analytetype = String, delim = '\t', levelname = nothing)
    endswith(file, ".dt") || throw(ArgumentError("The file is not a valid table directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    if config[:Type] == "R"
        analyte_name = Symbol(config[:Analyte])
        sample_name = Symbol.(filter!(!isempty, vectorize(config[:Sample])))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => Float64), types = Dict(analyte_name => analytetype), validate = false)
        for i in sample_name
            replace!(getproperty(tbl, i), missing => 0)
        end
        RowDataTable(getproperty(tbl, analyte_name), analyte_name, sample_name, tbl)
    else
        sample_name = Symbol(first(split(config[:Sample], "\t")))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => Float64), types = Dict(sample_name => String, levelname => Int), validate = false)
        analyte_name = String.(filter!(!isempty, vectorize(config[:Analyte])))
        config[:Type] == "C" || throw(ArgumentError("StackDataTable is not implemented yet."))
        for i in Symbol.(analyte_name)
            replace!(getproperty(tbl, i), missing => 0)
        end
        ColumnDataTable(analytetype.(analyte_name), sample_name, tbl)
    end
end

"""
    read_analysistable(file::String, T; tabletype = T, analytetype::Type{A} = String, delim = '\\t') -> AnalysisTable{A, S <: T}

Read ".at" file into julia as `AnalysisTable`. `T` is the sink function for tabular data, `tabletype` is `T` parameter in the type signature of `Batch` which determines the underlying table type
, and `analytetype` is a concrete type for `analyte` which msut have a method for string input, and `delim` specifies delimiter for tabular data if `config[:delim]` in `tables` do not exist.

If `tabletype` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".at" file.
"""
function read_analysistable(file::String, T; tabletype = T, analytetype = String, delim = '\t')
    endswith(file, ".at") || throw(ArgumentError("The file is not a valid table directory"))
    files = filter!(f -> endswith(f, ".dt"), readdir(file))
    tables = map(files) do f
        read_datatable(joinpath(file, f), T; analytetype, delim)
    end
    Cons = isnothing(tabletype) ? AnalysisTable : AnalysisTable{tabletype}
    Cons(Symbol.(replace.(files, Ref(".dt" => ""), Ref(r"^\d*_" => ""))), tables)
end
"""
    read_methodtable(file::String, T; tabletype = T, analytetype::Type{A} = String, delim = '\\t') -> MethodTable{A, S <: T}

Read ".mt" file into julia as `MethodTable`. `T` is the sink function for tabular data, `tabletype` is `T` parameter in the type signature of `MethodTable` which determines the underlying table type
, and `analytetype` is a concrete type for `analyte` which msut have a method for string input, and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

If `tabletype` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".mt" file.
"""
function read_methodtable(file::String, T; tabletype = T, analytetype = String, delim = '\t')
    endswith(file, ".mt") || throw(ArgumentError("The file is not a valid table directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    signal = Symbol(get(config, :signal, :area))
    analytetable = CSV.read(joinpath(file, "analytetable.txt"), Table)
    analyte = analytetype.(analytetable.analyte)
    isd = replace(analytetable.isd, missing => 0)
    conctable = read_datatable(joinpath(file, "true_concentration.dt"), T; analytetype, delim)
    if length(conctable.sample) > 1
        signaltable = read_datatable(joinpath(file, "$signal.dt"), T; analytetype, delim, levelname = get(config, :levelname, nothing))
        if haskey(config, :levelname) && in(Symbol(config[:levelname]), propertynames(signaltable))
            pointlevel = getproperty(signaltable, Symbol(config[:levelname]))
        elseif haskey(config, :pointlevel)
            pointlevel = parse.(Int, config[:pointlevel])
        else
            nl = length(conctable.sample)
            ns = length(signal.sample)
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
    Cons = isnothing(tabletype) ? MethodTable : MethodTable{tabletype}
    Cons(Table(analytetable; analyte, isd, calibration), signal, pointlevel, conctable, signaltable)
end
"""
    read_batch(file::String, T; tabletype = T, analytetype = String, delim = '\\t') -> Batch{A, tabletype}

Read ".batch" file into julia as `Batch`. `T` is the sink function for tabular data, `tabletype` is `T` parameter in the type signature of `Batch` which determines the underlying table type
, and `analytetype` is a concrete type for `analyte` which must have a method for string input, and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist

If `tabletype` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".batch" file.
"""
function read_batch(file::String, T; tabletype = T, analytetype = String, delim = '\t')
    endswith(file, ".batch") || throw(ArgumentError("The file is not a valid batch directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)   
    method = read_methodtable(joinpath(file, "method.mt"), T; tabletype, analytetype, delim)
    if !in("calibration", readdir(file)) || isempty(readdir(joinpath(file, "calibration")))
        if isnothing(method.signaltable)
            cal = [SingleCalibration((analyte, ), first(getanalyte(method.conctable, analyte))) for analyte in method.conctable.analyte]
        else
            cal = [calibration(method, analyte) for analyte in method.conctable.analyte if !isisd(method, analyte)]
        end
    else
        cal = [read_calibration(joinpath(file, "calibration", f); delim, analytetype) for f in readdir(joinpath(file, "calibration")) if endswith(f, ".mcal") || endswith(f, ".scal")]
    end
    fd = findfirst(==("data.at"), readdir(file))
    Cons = isnothing(tabletype) ? Batch : Batch{tabletype}
    Cons(method, cal,
        isnothing(fd) ? nothing : read_analysistable(joinpath(file, "data.at"), T; tabletype, analytetype, delim),
    )
end

function print_summary(io::IO, cal::MultipleCalibration{A}) where A
    print(io, "Calibration{$A} of ", first(cal.analyte), " with ", length(unique(cal.table.level[cal.table.include])), " levels and ", length(findall(cal.table.include)), " points")
end

function print_summary(io::IO, cal::SingleCalibration{A}) where A
    print(io, "Calibration{$A} of ", first(cal.analyte), " with single level")
end

function show(io::IO, ::MIME"text/plain", cal::MultipleCalibration)
    print(io, typeof(cal), " with ", length(unique(cal.table.level[cal.table.include])), " levels and ", length(findall(cal.table.include)), "points:\n")
    print(io, "∘ Analyte: ", first(cal.analyte), "\n")
    print(io, "∘ ISD: ", last(cal.analyte), "\n")
    print(io, "∘ Type: ", cal.type ? "linear\n" : "quadratic\n")
    print(io, "∘ (0, 0): ", cal.zero ? "included\n" : "ommitted\n")
    print(io, "∘ Weight: ", weight_repr(cal.weight), "\n")
    print(io, "∘ Formula: ", formula_repr(cal), "\n")
    print(io, "∘ R²: ", r2(cal.model))
end

function show(io::IO, ::MIME"text/plain", cal::SingleCalibration)
    print(io, typeof(cal), " with single level:\n")
    print(io, "∘ Analyte (ISD): ", first(cal.analyte), "\n")
    print(io, "∘ Concentration: ", cal.conc)
end

function print_summary(io::IO, batch::Batch{A, T}) where {A, T}
    print(io, "Batch{$A, $(shorten_type_repr(T))} with ", length(batch.method.analytetable.isd .< 0), " internal standards out of ", length(batch.method.analytetable.analyte), " analytes")
end

function show(io::IO, ::MIME"text/plain", batch::Batch{A, T}) where {A, T}
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
    length(t) > 50 ? replace(t, r"\{.*\}" => "{...}") : t
end

function print_summary(io::IO, tbl::ColumnDataTable{A, T}) where {A, T}
    print(io, "ColumnDataTable{$A, $(shorten_type_repr(T))} with ", length(tbl.analyte), " analytes and ", length(tbl.sample), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::ColumnDataTable)
    print_summary(io, tbl)
    println(io, "::")
    show(io, MIME"text/plain"(), tbl.table)
end

function print_summary(io::IO, tbl::RowDataTable{A, T}) where {A, T}
    print(io, "RowDataTable{$A, $(shorten_type_repr(T))} with ", length(tbl.analyte), " analytes and ", length(tbl.sample), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::RowDataTable)
    print_summary(io, tbl)
    println(io, "::")
    show(io, MIME"text/plain"(), tbl.table)
end

function print_summary(io::IO, tbl::AnalysisTable{A, T}) where {A, T}
    print(io, "AnalysisTable{$A, $(shorten_type_repr(T))} with ", length(tbl.analyte), " analytes and ", length(tbl.sample), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::AnalysisTable)
    print_summary(io, tbl)
    print(io, ":")
    for (k, v) in pairs(tbl.tables)
        print(io, "\n∘ ", k, " | ")
        show(io, MIME"text/plain"(), v)
    end
end

function print_summary(io::IO, tbl::MethodTable{A, T}) where {A, T}
    if isnothing(tbl.signaltable)
        print(io, "MethodTable{$A, $(shorten_type_repr(T))} with ", length(tbl.analytetable.analyte), " analytes")
    else
        print(io, "MethodTable{$A, $(shorten_type_repr(T))} with ", length(tbl.analytetable.analyte), " analytes, ", length(tbl.conctable.sample), " levels and ", length(tbl.signaltable.sample), " points")
    end
end

function show(io::IO, ::MIME"text/plain", tbl::MethodTable)
    print_summary(io, tbl)
    print(io, ":\n∘ Analytes:\n")
    show(io, MIME"text/plain"(), tbl.analytetable)
    print(io, "\n\n∘ Level: ", join(tbl.pointlevel, ", "), "\n")
    print(io, "∘ Concentration: \n")
    show(io, MIME"text/plain"(), tbl.conctable.table)
    print(io, "\n\n∘ Signal: \n")
    show(io, MIME"text/plain"(), isnothing(tbl.signaltable) ? nothing : tbl.signaltable.table)
end

function write(file::String, tbl::RowDataTable; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[Type]\nR\n\n[delim]\n", escape_string(string(delim)), "\n\n[Analyte]\n", tbl.analytecol, "\n\n[Sample]\n", join(tbl.sample, "\n"))
    end
    CSV.write(joinpath(file, "table.txt"), tbl.table; delim)
end

function write(file::String, tbl::ColumnDataTable; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[Type]\nC\n\n[delim]\n", escape_string(string(delim)), "\n\n[Analyte]\n", join(tbl.analytename, "\n"), "\n\n[Sample]\n", tbl.samplecol)
    end
    CSV.write(joinpath(file, "table.txt"), tbl.table; delim)
end

function write(file::String, tbl::AnalysisTable; delim = '\t')
    mkpath(file)
    rm.(readdir(file; join = true); recursive = true)
    for (i, (k, v)) in enumerate(pairs(tbl.tables))
        write(joinpath(file, "$(i - 1)_$k.dt"), v; delim)
    end
end

function write(file::String, tbl::MethodTable; delim = '\t')
    mkpath(file)
    write(joinpath(file, "true_concentration.dt"), tbl.conctable; delim)
    isnothing(tbl.signaltable) || write(joinpath(file, "$(tbl.signal).dt"), tbl.signaltable; delim)
    id = nothing
    if isa(tbl.signaltable, ColumnDataTable)
        id = findfirst(x -> getproperty(tbl.signaltable.table, x) == tbl.pointlevel, propertynames(tbl.signaltable.table))
    end
    open(joinpath(file, "config.txt"), "w+") do config
        isnothing(id) ? Base.write(config, "[signal]\n", tbl.signal, "\n\n[delim]\n", escape_string(string(delim)), "\n\n[pointlevel]\n", join(tbl.pointlevel, "\n")) : 
            Base.write(config, "[signal]\n", tbl.signal, "\n\n[delim]\n", escape_string(string(delim)), "\n\n[levelname]\n", propertynames(tbl.signaltable.table)[id])
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

Write `object` into ".scal" for `SingleCalibration`, ".mcal" for `MultipleCalibration`, ".at" for `AnalysisTable`, ".mt" for `MethodTable`, and ".batch" for `Batch`. 

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
    ChemistryQuantitativeAnalysis.read(file::String, T; tabletype = T, analytetype = String, delim = '\\t') -> Batch{A, tabletype}

Read ".scal" as `SingleCalibration`, ".mcal" as `MultipleCalibration`, ".at" as `AnalysisTable`, ".mt" as `MethodTable`, and ".batch" as `Batch`. 

`T` is the sink function for tabular data, `tabletype` is `T` parameter in the type signature which determines the underlying table type
, `analytetype` is a concrete type for `analyte` which msut have a method for string input, and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

If `tabletype` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".batch" file.
"""
function read(file::String, T; tabletype = T, analytetype = String, delim = '\t')
    if endswith(file, ".batch")
        read_batch(file, T; tabletype, analytetype, delim)
    elseif endswith(file, ".scal") || endswith(file, ".mcal")
        read_calibration(file; analytetype, delim)
    elseif endswith(file, ".at")
        read_analysistable(file, T; tabletype, analytetype, delim)
    elseif endswith(file, ".mt")
        read_methodtable(file, T; tabletype, analytetype, delim)
    elseif endswith(file, ".dt")
        read_datatable(file, T; analytetype, delim)
    end
end

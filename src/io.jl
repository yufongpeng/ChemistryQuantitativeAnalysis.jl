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
        config[:delim] = config[:delim] == "\\t" ? "\t" : config[:delim]
    end
    config
end  

"""
    read_calibration(file::String; delim = "\\t") -> AbstractCalibration{A}

Read ".mcal" or ".scal" file into julia as `MultipleCalibration` or `SingleCalibration`. `analyte_fn` is responsible for converting `analyte` in string type into user defined type `A`, and `delim` specifies delimiter for tabular data.

See README.md for the structure of ".mcal" and ".scal" file.
"""
function read_calibration(file::String; analyte_fn = default_analyte_fn, delim = "\t")
    endswith(file, ".mcal") || endswith(file, ".scal") || throw(ArgumentError("The file is not a valid calibration directory"))
    config = read_config(joinpath(file, "config.txt"))
    if endswith(file, ".scal")
        return SingleCalibration((analyte_fn(config[:analyte]), ), parse(Float64, config[:conc]))
    end
    tbl = CSV.read(joinpath(file, "table.txt"), Table; delim = get(config, :delim, delim))
    calibration(analyte_fn.((config[:analyte], config[:isd])), tbl; type = parse(Bool, config[:type]), zero = parse(Bool, config[:zero]), weight = parse(Float64, config[:weight]))
end

"""
    read_datatable(file::String, T; analyte_fn = default_analyte_fn, delim = "\\t") -> AbstractDataTable{A, S <: T}

Read ".dt" file into julia as `ColumnDataTable` or `RowDataTable`. `T` is the sink function for tabular data, `analyte_fn` is responsible for converting `analyte` in string type into user defined type `A`, and `delim` specifies delimiter for tabular data.

See README.md for the structure of ".dt" file.
"""
function read_datatable(file::String, T; analyte_fn = default_analyte_fn, delim = "\t")
    endswith(file, ".dt") || throw(ArgumentError("The file is not a valid table directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    if config[:Type] == "R"
        analyte_name = Symbol(config[:Analyte])
        sample_name = Symbol.(filter!(!isempty, vectorize(config[:Sample])))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => Float64), types = Dict(analyte_name => String), validate = false)
        for i in sample_name
            replace!(getproperty(tbl, i), missing => 0)
        end
        RowDataTable(sample_name, analyte_name, analyte_fn.(getproperty(tbl, analyte_name)), tbl)
    else
        sample_name = Symbol(first(split(config[:Sample], "\t")))
        tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => Float64), types = Dict(sample_name => String))
        analyte_name = String.(filter!(!isempty, vectorize(config[:Analyte])))

        analyte_name = String[]
        isd_map = Int[]
        for i in config[:Analyte]
            s = split(i, "\t")
            push!(analyte_name, first(s))
            push!(isd_map, (length(s) == 1 || isempty(last(s))) ? 0 : parse(Int, last(s)))
        end
        config[:Type] == "C" || throw(ArgumentError("StackDataTable is not implemented yet."))
        for i in Symbol.(analyte_name)
            replace!(getproperty(tbl, i), missing => 0)
        end
        ColumnDataTable(sample_name, Symbol.(analyte_name), analyte_fn.(analyte_name), tbl)
    end
end

"""
    read_analysistable(file::String, T; table_type = T, analyte_fn = default_analyte_fn, delim = "\\t") -> AnalysisTable{A, S <: T}

Read ".at" file into julia as `AnalysisTable`. `T` is the sink function for tabular data, `table_type` is `T` parameter in the type signature of `Batch` which determines the underlying table type
, and `analyte_fn` is responsible for converting `analyte` in string type into user defined type `A`, and `delim` specifies delimiter for tabular data.

If `table_type` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".at" file.
"""
function read_analysistable(file::String, T; table_type = T, analyte_fn = default_analyte_fn, delim = "\t")
    endswith(file, ".at") || throw(ArgumentError("The file is not a valid table directory"))
    files = filter!(f -> endswith(f, ".dt"), readdir(file))
    tables = map(files) do f
        read_datatable(joinpath(file, f), T; analyte_fn, delim)
    end
    Cons = isnothing(table_type) ? AnalysisTable : AnalysisTable{table_type}
    Cons(Symbol.(replace.(files, Ref(".dt" => ""), Ref(r"^\d*_" => ""))), tables)
end
"""
    read_methodtable(file::String, T; table_type = T, analyte_fn = default_analyte_fn, delim = "\\t") -> MethodTable{A, S <: T}

Read ".mt" file into julia as `MethodTable`. `T` is the sink function for tabular data, `table_type` is `T` parameter in the type signature of `MethodTable` which determines the underlying table type
, and `analyte_fn` is responsible for converting `analyte` in string type into user defined type `A`, and `delim` specifies delimiter for tabular data.

If `table_type` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".mt" file.
"""
function read_methodtable(file::String, T; table_type = T, analyte_fn = default_analyte_fn, delim = "\t")
    endswith(file, ".mt") || throw(ArgumentError("The file is not a valid table directory"))
    config = read_config(joinpath(file, "config.txt"))
    signal = Symbol(get(config, :signal, :area))
    analyte_map = CSV.read(joinpath(file, "analyte_map.txt"), Table)
    analytes = analyte_fn.(analyte_map.analytes)
    isd = replace(analyte_map.isd, missing => 0)
    conctable = read_datatable(joinpath(file, "true_concentration.dt"), T; analyte_fn, delim)
    if length(conctable.samples) > 1
        signaltable = read_datatable(joinpath(file, "$signal.dt"), T; analyte_fn, delim)
        if haskey(config, :level_map)
            level_map = parse.(Int, config[:level_map])
        else
            nl = length(conctable.samples)
            ns = length(signal.samples)
            nr = ns ÷ nl
            rr = ns - nl * nr
            level_map = vcat(repeat([1], rr), repeat(1:nl, nr))
        end
        calibration = map(enumerate(analyte_map.calibration)) do (i, c)
            ismissing(c) ? i : c
        end
    else
        signaltable = nothing
        level_map = [1]
        calibration = isd
    end
    Cons = isnothing(table_type) ? MethodTable : MethodTable{table_type}
    Cons(signal, analytes, isd, calibration, level_map, conctable, signaltable)
end
"""
    read_batch(file::String, T; table_type = T, analyte_fn = default_analyte_fn) -> Batch{A, table_type}

Read ".batch" file into julia as `Batch`. `T` is the sink function for tabular data, `table_type` is `T` parameter in the type signature of `Batch` which determines the underlying table type
, and `analyte_fn` is responsible for converting `analyte` in string type into user defined type `A`.

If `table_type` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".batch" file.
"""
function read_batch(file::String, T; table_type = T, analyte_fn = default_analyte_fn)
    endswith(file, ".batch") || throw(ArgumentError("The file is not a valid batch directory"))
    config = read_config(joinpath(file, "config.txt"))
    delim = config[:delim]   
    method = read_methodtable(joinpath(file, "method.mt"), T; table_type, analyte_fn, delim)
    if !in("calibration", readdir(file)) || isempty(readdir(joinpath(file, "calibration")))
        if isnothing(method.signaltable)
            cal = [SingleCalibration((analyte, ), first(get_analyte(method.conctable, analyte))) for analyte in method.conctable.analytes]
        else
            cal = [calibration(method, analyte) for analyte in method.conctable.analytes if !isisd(method, analyte)]
        end
    else
        cal = [read_calibration(joinpath(file, "calibration", f); delim, analyte_fn) for f in readdir(joinpath(file, "calibration")) if endswith(f, ".mcal") || endswith(f, ".scal")]
    end
    fd = findfirst(==("data.at"), readdir(file))
    Cons = isnothing(table_type) ? Batch : Batch{table_type}
    Cons(method, cal,
        isnothing(fd) ? nothing : read_analysistable(joinpath(file, "data.at"), T; table_type, analyte_fn, delim),
    )
end

function show(io::IO, cal::MultipleCalibration{A}) where A
    print(io, "Calibration{$A} of $(first(cal.analyte)) with $(length(unique(cal.table.level[cal.table.include]))) levels and $(length(findall(cal.table.include))) points")
end

function show(io::IO, cal::SingleCalibration{A}) where A
    print(io, "Calibration{$A} of $(first(cal.analyte)) with single level")
end

function show(io::IO, ::MIME"text/plain", cal::MultipleCalibration)
    print(io, cal, ":\n")
    print(io, "∘ Analyte: ", first(cal.analyte), "\n")
    print(io, "∘ ISD: ", last(cal.analyte), "\n")
    print(io, "∘ Type: ", cal.type ? "linear\n" : "quadratic\n")
    print(io, "∘ (0, 0): ", cal.zero ? "included\n" : "ommitted\n")
    print(io, "∘ Weight: ", weight_repr(cal.weight), "\n")
    print(io, "∘ Formula: ", formula_repr(cal), "\n")
    print(io, "∘ R²: ", r2(cal.model))
end

function show(io::IO, ::MIME"text/plain", cal::SingleCalibration)
    print(io, cal, ":\n")
    print(io, "∘ Analyte (ISD): ", first(cal.analyte), "\n")
    print(io, "∘ Concentration: ", cal.conc)
end

function show(io::IO, batch::Batch{A, T}) where {A, T}
    print(io, string("Batch{$A, $(shorten_type_repr(T))} with $(length(batch.method.analyte_map.isd .< 0)) internal standards out of $(length(batch.method.analyte_map.analytes)) analytes"))
end

function show(io::IO, ::MIME"text/plain", batch::Batch)
    print(io, batch, ":\n")
    print(io, "∘ Analytes:\n")
    show(io, MIME"text/plain"(), batch.method.analyte_map)
    print(io, "\n\n∘ Calibration:\n")
    show(io, MIME"text/plain"(), batch.calibration)
    print(io, "\n\n∘ Data:\n")
    show(io, MIME"text/plain"(), batch.data)
end

function shorten_type_repr(T)
    t = repr(T)
    length(t) > 50 ? replace(t, r"\{.*\}" => "{...}") : t
end

function show(io::IO, tbl::ColumnDataTable{A, T}) where {A, T}
    print(io, string("ColumnDataTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analytes)) analytes and $(length(tbl.samples)) samples"))
end

function show(io::IO, ::MIME"text/plain", tbl::ColumnDataTable)
    print(io, tbl, ":\n")
    show(io, MIME"text/plain"(), tbl.table)
end

function show(io::IO, tbl::RowDataTable{A, T}) where {A, T}
    print(io, string("RowDataTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analytes)) analytes and $(length(tbl.samples)) samples"))
end

function show(io::IO, ::MIME"text/plain", tbl::RowDataTable)
    print(io, tbl, ":\n")
    show(io, MIME"text/plain"(), tbl.table)
end

function show(io::IO, tbl::AnalysisTable{A, T}) where {A, T}
    print(io, string("AnalysisTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analytes)) analytes and $(length(tbl.samples)) samples"))
end

function show(io::IO, ::MIME"text/plain", tbl::AnalysisTable)
    print(io, tbl, ":")
    for (k, v) in pairs(tbl.tables)
        print(io, "\n∘ ", k, ": \n")
        show(io, MIME"text/plain"(), v)
    end
end

function show(io::IO, tbl::MethodTable{A, T}) where {A, T}
    isnothing(tbl.signaltable) && return print(io, string("MethodTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analyte_map.analytes)) analytes"))
    print(io, string("MethodTable{$A, $(shorten_type_repr(T))} with $(length(tbl.analyte_map.analytes)) analytes, $(length(tbl.conctable.samples)) levels and $(length(tbl.signaltable.samples)) points."))
end

function show(io::IO, ::MIME"text/plain", tbl::MethodTable)
    print(io, tbl, ":\n")
    print(io, "∘ Analytes:\n")
    show(io, MIME"text/plain"(), tbl.analyte_map)
    print(io, "\n∘ Level: ", join(tbl.level_map, ", "), "\n")
    print(io, "∘ Concentration: \n")
    show(io, MIME"text/plain"(), tbl.conctable.table)
    print(io, "\n∘ Signal: \n")
    show(io, MIME"text/plain"(), isnothing(tbl.signaltable) ? nothing : tbl.signaltable.table)
end

function write(file::String, tbl::RowDataTable; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[Type]\nR\n\n[Analyte]\n", tbl.analyte_name, "\n\n[Sample]\n", join(tbl.sample_name, "\n"))
    end
    CSV.write(joinpath(file, "table.txt"), tbl.table; delim)
end

function write(file::String, tbl::ColumnDataTable; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[Type]\nC\n\n[Analyte]\n", join(tbl.analyte_name, "\n"), "\n\n[Sample]\n", tbl.sample_name)
    end
    CSV.write(joinpath(file, "table.txt"), tbl.table; delim)
end

function write(file::String, tbl::AnalysisTable; delim = "\t")
    mkpath(file)
    rm.(readdir(file; join = true); recursive = true)
    for (i, (k, v)) in enumerate(pairs(tbl.tables))
        write(joinpath(file, "$(i - 1)_$k.dt"), v; delim)
    end
end

function write(file::String, tbl::MethodTable; delim = "\t")
    mkpath(file)
    write(joinpath(file, "true_concentration.dt"), tbl.conctable; delim)
    isnothing(tbl.signaltable) || write(joinpath(file, "$(tbl.signal).dt"), tbl.signaltable; delim)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[signal]\n", tbl.signal, "\n\n[level_map]\n", join(tbl.level_map, "\n"))
    end
    CSV.write(joinpath(file, "analyte_map.txt"), tbl.analyte_map; delim)
end

function write(file::String, cal::MultipleCalibration; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[analyte]\n", string(first(cal.analyte)), "\n\n[isd]\n", string(last(cal.analyte)), 
                "\n\n[type]\n", string(cal.type), "\n\n[zero]\n", string(cal.zero), "\n\n[weight]\n", string(cal.weight))
    end
    CSV.write(joinpath(file, "table.txt"), cal.table; delim)
end
function write(file::String, cal::SingleCalibration; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[analyte]\n", string(first(cal.analyte)), "\n\n[conc]\n", string(cal.conc))
    end
end

"""
    Calibration.write(file::String, object; delim = "\\t")

Write `object` into ".scal" for `SingleCalibration`, ".mcal" for `MultipleCalibration`, ".at" for `AnalysisTable`, ".mt" for `MethodTable`, and ".batch" for `Batch`. 

`delim` specifies delimiter for tabular data.

See README.md for the structure of all files.
"""
function write(file::String, batch::Batch; delim = "\t")
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[delim]\n", delim == "\t" ? "\\t" : ",")
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
    Calibration.read(file::String, T; table_type = T, analyte_fn = default_analyte_fn) -> Batch{A, table_type}

Read ".batch" file into julia as `Batch`. `T` is the sink function for tabular data, `table_type` is `T` parameter in the type signature of `Batch` which determines the underlying table type, and `analyte_fn` is responsible for converting `analyte` in string type into user defined type `A`.

If `table_type` is set to `nothing`, table type will be determined automatically which may be too restrict when using parameterized table types.

See README.md for the structure of ".batch" file.
"""
const read = read_batch
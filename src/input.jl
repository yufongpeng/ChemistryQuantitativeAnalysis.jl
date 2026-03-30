"""
    ChemistryQuantitativeAnalysis.read(file, T; analytetype = String, sampletype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\\t')

Read file into `ChemistryQuantitativeAnalysis` objects.

* `file::String`: result depending on extention
    * ".ical": `InternalCalibrator`
    * ".ecal": `ExternalCalibrator`
    * ".sdt": `SampleDataTable`
    * ".adt" `AnalyteDataTable`
    * ".at": `AnalysisTable`
    * ".am": `AnalysisMethod`
    * ".batch": `Batch`
* `T` is the sink function for tabular data.
* `analytetype` is a concrete type for `analyte`.
* `sampletype` is a concrete type for `sample.`
* `numbertype` is the type for numeric data, 
* `modeltype`: calibration model type.
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

For any `x` of type `analytetype` or `sampletype`, `x == cqaconvert(type, string(x)))`. 
Additionally, `tryparse` have to be extended for `CSV` parsing:
* `tryparse(::Type{analytetype}, x::AbstractString)` is neccessary for `AnalyteDataTable`.
* `tryparse(::Type{sampletype}, x::AbstractString)` is neccessary for `SampleDataTable`.

See README.md for the structure of ".batch" file.
"""
function read(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\t', kwargs...)
    if endswith(file, ".batch")
        read_batch(file, T; analytetype, sampletype, numbertype, modeltype, delim, kwargs...)
    elseif endswith(file, ".ical") || endswith(file, ".ecal")
        read_calibrator(file; analytetype, numbertype, delim, kwargs...)
    elseif endswith(file, ".at")
        read_analysistable(file, T; analytetype, sampletype, numbertype, delim, kwargs...)
    elseif endswith(file, ".am")
        read_method(file, T; analytetype, sampletype, numbertype, delim, modeltype, kwargs...)
    elseif endswith(file, ".sdt") || endswith(file, ".adt")
        read_datatable(file, T; analytetype, sampletype, numbertype, delim, kwargs...)
    end
end

function read_config(file::String)
    d = dirname(file)
    in(basename(file), readdir(isempty(d) ? pwd() : d)) || throw(Base.IOError(string("read(", file, "): no such file or directory (ENOENT)"), -4058))
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
    read_calibrator(file::String; analytetype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\\t') -> AbstractCalibrator{analytetype, numbertype}

Read file as calibrator.

* `file::String`: result depending on extention
    * ".ical": `InternalCalibrator`
    * ".ecal": `ExternalCalibrator`
* `analytetype` is a concrete type for `analyte`.
* `numbertype` is the type for numeric data, 
* `modeltype`: calibration model type.
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

`numbertype` is forced to `Float64` because of issue of `GLM`.
    
See `ChemistryQuantitativeAnalysis.read` for the requirement of `analytetype`.

See README.md for the structure of ".ecal" and ".ical" file.
"""
function read_calibrator(file::String; analytetype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\t')
    d = dirname(file)
    in(basename(file), readdir(isempty(d) ? pwd() : d)) || throw(Base.IOError(string("read(", file, "): no such file or directory (ENOENT)"), -4058))
    endswith(file, ".ecal") || endswith(file, ".ical") || throw(ArgumentError("The file is not a valid calibrator directory"))
    # numbertype = Float64 # Bug of GLM
    config = read_config(joinpath(file, "config.txt"))
    if endswith(file, ".ical")
        return InternalCalibrator(analytetype(config[:analyte]), parse(numbertype, config[:conc]))
    end
    tbl = CSV.read(joinpath(file, "table.txt"), Table; delim = get(config, :delim, delim), typemap = Dict(Int => numbertype, Float64 => numbertype), types = Dict(:level => Int), stringtype = String)
    model = cqaconvert(modeltype, config[:model])
    externalcalibrate(cqaconvert(analytetype, config[:analyte]), config[:isd] == "nothing" ? nothing : cqaconvert(analytetype, config[:isd]), tbl, model; decode_config(modeltype, config)...)
end

"""
    read_datatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t', levelname = :level) -> AbstractDataTable{analytetype, sampletype, numbertype}

Read file as data table.

* `file::String`: result depending on extention
    * ".sdt": `SampleDataTable`
    * ".adt" `AnalyteDataTable`
* `T` is the sink function for tabular data.
* `analytetype` is a concrete type for `analyte`.
* `sampletype` is a concrete type for `sample.`
* `numbertype` is the type for numeric data, 
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".sdt" and ".adt" file.
"""
function read_datatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t', levelname = nothing)
    endswith(file, ".sdt") ? read_sampledatatable(file, T; analytetype, sampletype, numbertype, delim, levelname) : 
    endswith(file, ".adt") ? read_analytedatatable(file, T; analytetype, sampletype, numbertype, delim) : 
    throw(ArgumentError("The file is not a valid table directory"))
end

"""
    read_sampledatatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t', levelname = :level) -> SampleDataTable{analytetype, sampletype, numbertype}

Read file into `SampleDataTable`.

* `file::String`: it must have extention of ".sdt"
* `T` is the sink function for tabular data.
* `analytetype` is a concrete type for `analyte`.
* `sampletype` is a concrete type for `sample.`
* `numbertype` is the type for numeric data, 
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.
* `level` is specifically used for AnalysisMethod, indicating the column representing calibration level; this column should be all integers. 

See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".sdt" file.
"""
function read_sampledatatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t', levelname = nothing)
    endswith(file, ".sdt") || throw(ArgumentError("The file is not a valid table directory"))
    d = dirname(file)
    in(basename(file), readdir(isempty(d) ? pwd() : d)) || throw(Base.IOError(string("read(", file, "): no such file or directory (ENOENT)"), -4058))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    sample_name = Symbol(first(split(config[:Sample], "\t")))
    tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => numbertype), types = Dict(sample_name => sampletype, levelname => Int), stringtype = String, validate = false)
    analyte_name = String.(filter!(!isempty, vectorize(config[:Analyte])))
    for i in Symbol.(analyte_name)
        replace!(getproperty(tbl, i), missing => 0)
    end
    SampleDataTable(cqamap(analytetype, analyte_name), getproperty(tbl, sample_name), sample_name, tbl)
end

"""
    read_analytedatatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t') -> AnalyteDataTable{analytetype, sampletype, numbertype}

Read file into `AnalyteDataTable`.

* `file::String`: it must have extention of ".adt"
* `T` is the sink function for tabular data.
* `analytetype` is a concrete type for `analyte`.
* `sampletype` is a concrete type for `sample.`
* `numbertype` is the type for numeric data, 
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".adt" file.
"""
function read_analytedatatable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t')
    endswith(file, ".adt") || throw(ArgumentError("The file is not a valid table directory"))
    d = dirname(file)
    in(basename(file), readdir(isempty(d) ? pwd() : d)) || throw(Base.IOError(string("read(", file, "): no such file or directory (ENOENT)"), -4058))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    analyte_name = Symbol(config[:Analyte])
    sample_name = String.(filter!(!isempty, vectorize(config[:Sample])))
    tbl = CSV.read(joinpath(file, "table.txt"), T; delim, typemap = Dict(Int => numbertype), types = Dict(analyte_name => analytetype), stringtype = String, validate = false)
    for i in Symbol.(sample_name)
        replace!(getproperty(tbl, i), missing => 0)
    end
    AnalyteDataTable(getproperty(tbl, analyte_name), cqamap(sampletype, sample_name), analyte_name, tbl)

end

"""
    read_analysistable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\\t') -> AnalysisTable{analytetype, sampletype}

Read file into `AnalysisTable`.

* `file::String`: it must have extention of ".at"
* `T` is the sink function for tabular data.
* `analytetype` is a concrete type for `analyte`.
* `sampletype` is a concrete type for `sample.`
* `numbertype` is the type for numeric data, 
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.
    
See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".at" file.
"""
function read_analysistable(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, delim = '\t')
    endswith(file, ".at") || throw(ArgumentError("The file is not a valid table directory"))
    d = dirname(file)
    in(basename(file), readdir(isempty(d) ? pwd() : d)) || throw(Base.IOError(string("read(", file, "): no such file or directory (ENOENT)"), -4058))
    files = filter!(f -> endswith(f, r".[a,s]dt"), readdir(file))
    tables = map(files) do f
        read_datatable(joinpath(file, f), T; analytetype, sampletype, numbertype, delim)
    end
    AnalysisTable(Symbol.(replace.(files, Ref(r".[a,s]dt" => ""), Ref(r"^\d*_" => ""))), tables)
end
"""
    read_method(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\\t') -> AnalysisMethod{analytetype, <: Table}

Read file into `AnalysisMethod`.

* `file::String`: it must have extention of ".am"
* `T` is the sink function for tabular data.
* `analytetype` is a concrete type for `analyte`.
* `sampletype` is a concrete type for `sample.`
* `numbertype` is the type for numeric data, 
* `modeltype`: calibration model type.
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.
  
See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".am" file.
"""
function read_method(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\t')
    endswith(file, ".am") || throw(ArgumentError("The file is not a valid method directory"))
    d = dirname(file)
    in(basename(file), readdir(isempty(d) ? pwd() : d)) || throw(Base.IOError(string("read(", file, "): no such file or directory (ENOENT)"), -4058))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)
    signal = Symbol(get(config, :signal, :area))
    rel_sig = Symbol(get(config, :rel_sig, :relative_signal))
    est_conc = Symbol(get(config, :est_conc, :estimated_concentration))
    nom_conc = Symbol(get(config, :nom_conc, :nominal_concentration))
    acc = Symbol(get(config, :acc, :accuracy))
    analytetable = CSV.read(joinpath(file, "analytetable.txt"), Table; stringtype = String)
    analyte = analytetype.(analytetable.analyte)
    isd = :isd in propertynames(analytetable) ? replace(analytetable.isd, missing => 0) : zeros(length(analytetable))
    std = :std in propertynames(analytetable) ? map(enumerate(analytetable.std)) do (i, c)
            ismissing(c) ? i : c
        end : nothing
    model = :model in propertynames(analytetable) ? cqaconvert.(modeltype, analytetable.model) : [LinearCalibrator(ConstWeight()) for _ in eachindex(analytetable)]
    model = convert(Vector{modeltype}, model)
    conctable = read_datatable(joinpath(file, in("$nom_conc.sdt", readdir(file)) ? "$nom_conc.sdt" : "$nom_conc.adt"), T; analytetype, sampletype = Int, numbertype, delim)
    if length(sampleobj(conctable)) > 1
        signaltable = read_datatable(joinpath(file, in("$signal.sdt", readdir(file)) ? "$signal.sdt" : "$signal.adt"), T; analytetype, sampletype, numbertype, delim, levelname = get(config, :levelname, nothing))
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
        if isnothing(std)
            std = collect(eachindex(analytetable))
        end
    else
        signaltable = nothing
        pointlevel = [1]
        if isnothing(std)
            std = isd
        end
        # std = [v > 0 ? v : i for (i, v) in enumerate(isd)]
    end
    AnalysisMethod(Table(analytetable; analyte, isd, std, model), signal, rel_sig, est_conc, nom_conc, acc, pointlevel, conctable, signaltable)
end
"""
    read_batch(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\\t') -> Batch{analytetype}

Read file into `Batch`.

* `file::String`: it must have extention of ".batch"
* `T` is the sink function for tabular data.
* `analytetype` is a concrete type for `analyte`.
* `sampletype` is a concrete type for `sample.`
* `numbertype` is the type for numeric data, 
* `modeltype`: calibration model type.
* `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.
   
See `ChemistryQuantitativeAnalysis.read` for the requirements of `analytetype` and `sampletype`.

See README.md for the structure of ".batch" file.
"""
function read_batch(file::String, T; analytetype = String, sampletype = String, numbertype = Float64, modeltype = CalibrationModel, delim = '\t')
    endswith(file, ".batch") || throw(ArgumentError("The file is not a valid batch directory"))
    d = dirname(file)
    in(basename(file), readdir(isempty(d) ? pwd() : d)) || throw(Base.IOError(string("read(", file, "): no such file or directory (ENOENT)"), -4058))
    config = read_config(joinpath(file, "config.txt"))
    delim = get(config, :delim, delim)   
    method = read_method(joinpath(file, "method.am"), T; analytetype, sampletype, numbertype, modeltype, delim)
    if !in("calibrator", readdir(file)) || isempty(readdir(joinpath(file, "calibrator")))
        cal = (AbstractCalibrator{T} where {T <: typeof(method).parameters[1]})[]
        # if isnothing(method.signaltable)
            # cal = [InternalCalibrator(analyte, mean(getanalyte(method.conctable, analyte))) for analyte in method.calanalyte]
        # else
            # cal = [calibrate(method, analyte) for analyte in method.calanalyte]
        # end
    else
        cal = [read_calibrator(joinpath(file, "calibrator", f); delim, analytetype, numbertype, modeltype) for f in readdir(joinpath(file, "calibrator")) if endswith(f, ".ecal") || endswith(f, ".ical")]
    end
    fd = findfirst(==("data.at"), readdir(file))
    batch = Batch(method, cal,
        isnothing(fd) ? nothing : read_analysistable(joinpath(file, "data.at"), T; analytetype, sampletype, numbertype, delim),
    )
    for analyte in batch.analyte
        j = findfirst(x -> ==(x.analyte, analyte), batch.calibrator)
        (isnothing(j) || batch.calibrator[j] isa InternalCalibrator) && continue
        if batch.method.analytetable.model[j] != batch.calibrator[j].model
            model_calibrator!(batch.calibrator[j], batch.method.analytetable.method[j])
        end
    end
    batch
end

"""
    decode_config(::Type{<: CalibrationModel}, config::Dictionary) -> NamedTuple

Decode config for keyword arguments of `mkcalmodel` for certain calibration model type.
"""
function decode_config(::Type{<: CalibrationModel}, config::Dictionary)
    (; )
end

"""
    cqaconvert(::Type{T}, x::S) -> T
    cqaconvert(fn::Function, x::S)

Call constructor, `parse`, or `fn`. Extend this function if defining `T(x)` is not safe.

Return
* `parse(T, x)` if `T <: Number` and `S <: Union{AbstractString, AbstractChar}`.
* `fn(x)` if `fn` is a `Function`.
* `T(x)` otherwise.

Note that for any `x` of type `T`, `x == cqaconvert(T, string(x)))`.
"""
cqaconvert(::Type{T}, x) where T = T(x)::T
cqaconvert(::Type{T}, x::T) where T = x
cqaconvert(::Type{T}, x::Union{AbstractString, AbstractChar}) where {T <: Number} = parse(T, x)::T
cqaconvert(fn::T, x) where {T <: Function} = fn(x)

"""
    cqamap(::Type{T}, x::AbstractVector{S}) -> AbstractVector{T}
    cqamap(fn::Function, x::AbstractVector)

Convert `x` to a vector of desired type. If direct construction is not possible, i.e. neither `T <: S` nor `S <: T`, it applys `cqaconvert` on every elements.

For a function, `T` is inferred by `cqatype(fn)` first and then `cqatype(fn, v)` for the returned vector `v` to avoid abstract element type. 
"""
cqamap(::Type{T}, x::AbstractVector{T}) where T = x
cqamap(::Type{T}, x::AbstractVector{<: T}) where T = Vector{T}(x)
cqamap(::Type{T}, x::AbstractVector{S}) where {S, T <: S} = Vector{T}(x)
cqamap(::Type{T}, x::AbstractVector) where T = typedmap(e -> cqaconvert(T, e), T, x)
function cqamap(fn::T, x::AbstractVector) where {T <: Function}
    R = cqatype(fn)
    if isabstracttype(R)
        v = map(e -> cqaconvert(fn, e), x)
        cqamap(cqatype(fn, v), v)
    else
        typedmap(e -> cqaconvert(fn, e), R, x)
    end
end

"""
    cqatype(fn::Function)
    cqatype(fn::Function, x::AbstractVector) -> Type

Default returned type of `cqaconvert(fn, x)`. By default, the union of types of each element is used to avoid abstract element type. Extend this function if neccessary.
"""
cqatype(fn::T) where {T <: Function} = Any
cqatype(fn::T, v::AbstractVector) where {T <: Function} = Union{typeof.(v)...}

function cqaconvert(::Type{T}, x::AbstractString) where {T <: ComposedWeight}
    m = match(r"^[a-z,A-Z,\{\}]*\(\)$", x)
    isnothing(m) && throw(ArgumentError("Cannot convert $x to type $T"))
    r = eval(Meta.parse(x))
    r isa T ? r : throw(ArgumentError("Cannot convert $x to type $T"))
end

function cqaconvert(::Type{T}, x::AbstractString) where {T <: CalibrationModel}
    m = match(r"^CalibrationModel\{[a-z,A-Z,\d]*\}\([a-z,A-Z,\d]*\(.*\)\)$", x)
    m = isnothing(m) ? match(r"^[a-z,A-Z,\d]*Calibrator\([a-z,A-Z,\d]*\(.*\)\)$", x) : m
    isnothing(m) && throw(ArgumentError("Cannot convert $x to type $T"))
    r = eval(Meta.parse(x))
    r isa T ? r : throw(ArgumentError("Cannot convert $x to type $T"))
end
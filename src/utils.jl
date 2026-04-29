"""
    typedmap(::Type{T}, c...) 
    typedmap(f, ::Type{T}, c...)

Transform collection `c` by applying `f` to each element. For multiple collection arguments, apply `f` elementwise, and stop when any of them  
is exhausted. The element type of returned collection will be forced to `T`. In addtion, `typedmap(T, c...)` is equivalent to `typedmap(T, T, c...)`.
"""
typedmap(::Type{T}, iters...) where T = collect(T, Base.Generator(T, iters...))
typedmap(fn::F, ::Type{T}, iters...) where {F, T} = collect(T, Base.Generator(fn, iters...))

"""
    stdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return calibration standard of `analyte` based on `method`.
"""
function stdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = getanalyteid(method, analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = method.analytetable.std[aid]
    (iid > 0 ? method.analytetable.analyte[iid] : nothing)::Union{A, Nothing}
end
"""
    isdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return internal standard of `analyte` based on `method`.
"""
function isdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = getanalyteid(method, analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = method.analytetable.isd[aid]
    (iid > 0 ? method.analytetable.analyte[iid] : iid == 0 ? nothing : analyte)::Union{A, Nothing}
end
"""
    isstd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return if `analyte` is a internal standard based on `method`.
"""
function isstd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    !isnothing(findfirst(==(aid), method.analytetable.std))
end
"""
    isisd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return if `analyte` is a internal standard based on `method`.
"""
function isisd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    method.analytetable.isd[aid] < 0
end

# getter method
"""
    table(dt::AbstractDataTable)

Get field `:table`. 
"""
table(dt::AbstractDataTable) = getfield(dt, :table)
"""
    tables(dt::AnalysisTable)

Get field `:tables`. 
"""
tables(at::AnalysisTable) = getfield(at, :tables)
"""
    table_convert(::Type{T}, data::AbstractDataTable)

Convert internal table to type `T`.
"""
table_convert(::Type{D}, data::AbstractDataTable{A, S, N, T}) where {D, A, S, N, T <: D} = data
table_convert(::Type{D}, data::SampleDataTable) where D = 
    SampleDataTable(analyteobj(data), sampleobj(data), idcol(data), D(table(data)))
table_convert(::Type{D}, data::AnalyteDataTable) where D = 
    AnalyteDataTable(analyteobj(data), sampleobj(data), idcol(data), D(table(data)))
"""
    analyteobj(dt::AbstractDataTable{A}) -> Vector{A}
    analyteobj(at::AnalysisTable{A}) -> Vector{A}

Equivalent to `getfield(dt, :analyte)` or `getfield(at, :analyte)`.
"""
analyteobj(dt::AbstractDataTable) = getfield(dt, :analyte)
analyteobj(at::AnalysisTable) = getfield(at, :analyte)
"""
    sampleobj(dt::AbstractDataTable{A, S}) -> Vector{S}
    sampleobj(at::AnalysisTable{A, S, T}) -> Vector{S}

Equivalent to `getfield(dt, :sample)` or `getfield(at, :sample)`.
"""
sampleobj(dt::AbstractDataTable) = getfield(dt, :sample)
sampleobj(at::AnalysisTable) = getfield(at, :sample)
"""
    analytename(dt::AbstractDataTable) -> Vector{Symbol}
    analytename(at::AnalysisTable) -> Vector{Symbol}

Equivalent to `Symbol.(analyteobj(dt))` or `Symbol.(analyteobj(at))`.
"""
analytename(dt::AbstractDataTable) = Symbol.(analyteobj(dt))
analytename(at::AnalysisTable) = Symbol.(analyteobj(at))
"""
    samplename(dt::AbstractDataTable) -> Vector{Symbol}
    samplename(at::AnalysisTable) -> Vector{Symbol}

Equivalent to `Symbol.(sampleobj(dt))` or `Symbol.(sampleobj(at))`.
"""
samplename(dt::AbstractDataTable) = Symbol.(sampleobj(dt))
samplename(at::AnalysisTable) = Symbol.(sampleobj(at))
"""
    samplecol(dt::SampleDataTable) -> Symbol

Equivalent to `getfield(dt, :samplecol)`.
"""
samplecol(dt::SampleDataTable) = getfield(dt, :samplecol)
"""
    analytecol(dt::SampleDataTable) -> Symbol

Equivalent to `getfield(dt, :analytecol)`.
"""
analytecol(dt::AnalyteDataTable) = getfield(dt, :analytecol)
"""
    idcol(dt::AbstractDataTable) -> Symbol

Unified interface for identity column. Equivalent to `samplecol(dt::SampleDataTable)` or `analytecol(dt::SampleDataTable)`.
"""
idcol(dt::SampleDataTable) = samplecol(dt)
idcol(dt::AnalyteDataTable) = analytecol(dt)

"""
    findanalyte(dt::AbstractDataTable{A}, analyte::B) where {A, B <: A}
    findanalyte(dt::AbstractDataTable, analyte::Symbol)

Return the index of the first element of `analyteobj(dt)` or `analytename(dt)` for which the element equals to `analyte`.
"""
findanalyte(dt::AbstractDataTable{A}, analyte::B) where {A, B <: A} = findfirst(==(analyte), analyteobj(dt))
findanalyte(dt::AbstractDataTable, analyte::Symbol) = findfirst(==(analyte), analytename(dt))

"""
    findsample(dt::AbstractDataTable{A, S}, sample::T) where {A, S, T <: S}
    findsample(dt::AbstractDataTable, sample::Symbol)

Return the index of the first element of `sampleobj(dt)` or `samplename(dt)` for which the element equals to `sample`.
"""
findsample(dt::AbstractDataTable{A, S}, sample::T) where {A, S, T <: S} = findfirst(==(sample), sampleobj(dt))
findsample(dt::AbstractDataTable, sample::Symbol) = findfirst(==(sample), samplename(dt))

"""
    getanalyte(dt::AnalyteDataTable, id::Int)
    getanalyte(dt::SampleDataTable, id::Int)
    getanalyte(dt::AnalyteDataTable, analyte)
    getanalyte(dt::SampleDataTable, analyte)

Get data belonging to `analyte` or `analyteobj(dt)[id]` as a `Vector`. For `AnalyteDataTable`, a new vector is created; mutating this vector will not change the value in `dt`.
"""
getanalyte(dt::AnalyteDataTable{A, S, N}, id::Int) where {A, S, N} = [getproperty(dt, p)[id] for p in samplename(dt)]::Vector{N}
function getanalyte(dt::AnalyteDataTable{A, S, N}, analyte) where {A, S, N}
    id = findanalyte(dt, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    [getproperty(dt, p)[id] for p in samplename(dt)]::Vector{N}
end
getanalyte(dt::SampleDataTable{A, S, N}, id::Int) where {A, S, N} = getproperty(dt, analytename(dt)[id])::AbstractVector{N}
function getanalyte(dt::SampleDataTable{A, S, N}, analyte) where {A, S, N}
    id = findanalyte(dt, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    getproperty(dt, Symbol(analyteobj(dt)[id]))::AbstractVector{N}
end

"""
    getsample(dt::AnalyteDataTable, id::Int)
    getsample(dt::SampleDataTable, id::Int)
    getsample(dt::AnalyteDataTable, sample)
    getsample(dt::SampleDataTable, sample)

Get data belonging to `sample` or `sampleobj(dt)[id]` as a `Vector`. For `SampleDataTable`, a new vector is created; mutating this vector will not change the value in `dt`.
"""
getsample(dt::AnalyteDataTable{A, S, N}, id::Int) where {A, S, N} = getproperty(dt, Symbol(sampleobj(dt)[id]))::AbstractVector{N}
function getsample(dt::AnalyteDataTable{A, S, N}, sample) where {A, S, N}
    id = findsample(dt, sample)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    getproperty(dt, Symbol(sampleobj(dt)[id]))::AbstractVector{N}
end
getsample(dt::SampleDataTable{A, S, N}, id::Int) where {A, S, N} = [getproperty(dt, p)[id] for p in analytename(dt)]::Vector{N}
function getsample(dt::SampleDataTable{A, S, N}, sample) where {A, S, N}
    id = findsample(dt, sample)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    [getproperty(dt, p)[id] for p in analytename(dt)]::Vector{N}
end

"""
    findcalibrator(batch::Batch, analyte)
    findcalibrator(calibrator::AbstractVector, analyte)
    findcalibrator(calibrator::AbstractVector, analyte::Symbol)

Return the index of the first element of `calibrator` or `batch.calibrator` for which the element is the calibrator of `analyte`.
"""
findcalibrator(batch::Batch, analyte) = findcalibrator(batch.calibrator, analyte)
function findcalibrator(calibrator::AbstractVector{A}, analyte) where {A <: AbstractCalibrator}
    findfirst(x -> x.analyte == analyte, calibrator)
end
function findcalibrator(calibrator::AbstractVector{A}, analyte::Symbol) where {A <: AbstractCalibrator}
    findfirst(x -> Symbol(x.analyte) == analyte, calibrator)
end

"""
    getcalibrator(batch::Batch, analyte)
    getcalibrator(calibrator::AbstractVector, analyte)
    getcalibrator(calibrator::AbstractVector, analyte::Symbol)

Get calibrator of `analyte` from `calibrator` or `batch.calibrator`.
"""
getcalibrator(batch::Batch, analyte) = getcalibrator(batch.calibrator, analyte)
function getcalibrator(calibrator::AbstractVector{A}, analyte) where {A <: AbstractCalibrator}
    id = findcalibrator(calibrator, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte does not have a calibration curve."))
    calibrator[id]
end
getcalibrator(batch::Batch, id::Int) = getcalibrator(batch.calibrator, id)
getcalibrator(calibrator::AbstractVector{A}, id::Int) where {A <: AbstractCalibrator} = calibrator[id]

"""
    dynamic_range(cal::AbstractCalibrator)

Return dynamic range as a `Tuple` (lloq, uloq).
"""
dynamic_range(cal::AbstractCalibrator) = (lloq(cal), uloq(cal))

"""
    signal_range(cal::AbstractCalibrator)

Return theoretical signal of dynamic range as a `Tuple` (lloq, uloq).
"""
signal_range(cal::AbstractCalibrator) = (predict(cal, collect(dynamic_range(cal)))..., )

"""
    lloq(cal::AbstractCalibrator)

Return lower limit of quantification.
"""
function lloq(cal::ExternalCalibrator) 
    i = findfirst(cal.table.include)
    isnothing(i) ? eltype(cal.table.x)(NaN) : cal.table.x[i]
end
lloq(cal::InternalCalibrator{A, N}) where {A, N} = N(0)

"""
    uloq(cal::AbstractCalibrator)

Return upper limit of quantification.
"""
function uloq(cal::ExternalCalibrator) 
    i = findlast(cal.table.include)
    isnothing(i) ? eltype(cal.table.x)(NaN) : cal.table.x[i]
end
uloq(cal::InternalCalibrator{A, N}) where {A, N} = N(Inf)

"""
    signal_lloq(cal::AbstractCalibrator)

Return theoretical signal of lower limit of quantification.
"""
signal_lloq(cal::AbstractCalibrator) = only(predict(cal, [lloq(cal)]))

"""
    signal_uloq(cal::AbstractCalibrator)

Return theoretical signal of upper limit of quantification.
"""
signal_uloq(cal::AbstractCalibrator) = only(predict(cal, [uloq(cal)]))

"""
    blank(cal::AbstractCalibrator)
    blank(model::CalibrationModel)

Blank signal of `cal`.
"""
blank(cal::ExternalCalibrator) = blank(cal.model)
blank(cal::InternalCalibrator) = 0
blank(model::CalibrationModel{Proportional}) = 0
blank(model::CalibrationModel{Linear}) = 0
blank(model::CalibrationModel{QuadraticOrigin}) = 0
blank(model::CalibrationModel{Quadratic}) = 0
blank(model::CalibrationModel{Logarithmic}) = -Inf
blank(model::CalibrationModel{Exponential}) = 0
blank(model::CalibrationModel{Power}) = 0

"""
    signal_lob(cal::AbstractCalibrator)

Theoretical limit of blank (LOB) in signal space.
"""
function signal_lob(cal::ExternalCalibrator)
    zero_id = findall(==(0), cal.table.level)
    if isempty(zero_id) && isnan(lloq(cal))
        eltype(cal.table.y)(NaN)
    elseif isempty(zero_id)
        blank(cal) + 1.645 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]]) 
    else 
        mean(cal.table.y[zero_id]) + 1.645 * std(cal.table.y[zero_id])
    end
end
signal_lob(cal::InternalCalibrator) = blank(cal)

"""
    signal_lod(cal::AbstractCalibrator)

Theoretical limit of detection (LOD) in signal space.
"""
function signal_lod(cal::ExternalCalibrator) 
    zero_id = findall(==(0), cal.table.level)
    if isempty(zero_id) && isnan(lloq(cal))
        eltype(cal.table.y)(NaN)
    elseif isempty(zero_id)
        blank(cal) + 3.3 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]])
    else 
        signal_lob(cal) + 1.645 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]])
    end
end
signal_lod(cal::InternalCalibrator) = signal_lob(cal)

"""
    signal_loq(cal::AbstractCalibrator)

Theoretical limit of quantification (LOQ) in signal space.
"""
function signal_loq(cal::ExternalCalibrator) 
    if isnan(lloq(cal))
        eltype(cal.table.y)(NaN)
    else
        blank(cal) + 10 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]])
    end
end
signal_loq(cal::InternalCalibrator) = signal_lod(cal)

"""
    getweight(cal::InternalCalibrator)
    getweight(cal::ExternalCalibrator)
    getweight(model::CalibrationModel)

Get weight object from `cal` or `model`. 
"""
getweight(cal::InternalCalibrator) = ConstWeight()
getweight(cal::ExternalCalibrator) = getweight(cal.model)
getweight(model::CalibrationModel) = model.weight

"""
    findoutofrange(batch::Batch[, analyte]; data = batch.method.est_conc, limit = dynamic_range, dev = (0.2, 0.15)) -> Vector{Pair}

Find all samples that have quantity `data` out of acceptable range for each analyte. The return vector is analyte-id vector pairs. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable range.
* `dev::Tuple`: allowed devience of range.
"""
findoutofrange(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = dynamic_range, dev = (0.2, 0.15)) = _findoutofrange(batch, analyte, data, isoutofrange, limit, dev)

"""
    findunderrange(batch::Batch[, analyte]; data = batch.method.est_conc, limit = lloq, dev = 0.2) -> Vector{Pair}

Find all samples that have quantity `data` under acceptable range for each analyte. The return vector is analyte-id vector pairs. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable lower limit.
* `dev::Real`: allowed devience of limit.
"""
findunderrange(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = lloq, dev = 0.2) = _findoutofrange(batch, analyte, data, isunderrange, limit, dev)

"""
    findoverrange(batch::Batch[, analyte]; data = batch.method.est_conc, limit = uloq, dev = 0.15) -> Vector{Pair}

Find all samples that have quantity `data` over acceptable range for each analyte. The return vector is analyte-id vector pairs. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable upper limit.
* `dev::Tuple`: allowed devience of limit.
"""
findoverrange(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = uloq, dev = 0.15) = _findoutofrange(batch, analyte, data, isoverrange, limit, dev)

"""
    markoutofrange(batch::Batch[, analyte]; data = batch.method.est_conc, limit = dynamic_range, dev = (0.2, 0.15), value = (x, y) -> (NaN, NaN)) -> AbstractDataTable

Mark all samples that have quantity `data` out of acceptable range with `value` for each analyte. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable range.
* `dev::Tuple`: allowed devience of range.
* `value::Function`: 2-arg functions returning marked values for each analyte. The first argument is calibrator; the second is sample data vector.
"""
markoutofrange(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = dynamic_range, dev = (0.2, 0.15), value = (x, y) -> (NaN, NaN)) = 
    _markoutofrange(batch, analyte, data, (isunderrange, isoverrange), limit, dev, value)

"""
    markunderrange(batch::Batch[, analyte]; data = batch.method.est_conc, limit = lloq, dev = 0.2, value = (x, y) -> NaN) -> AbstractDataTable

Mark all samples that have quantity `data` under acceptable range with `value` for each analyte. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable lower limit.
* `dev::Real`: allowed devience of limit.
* `value::Function`: 2-arg function returning marked value for each analyte. The first argument is calibrator; the second is sample data vector.
"""
markunderrange(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = lloq, dev = 0.2, value = (x, y) -> NaN) = 
    _markoutofrange(batch, analyte, data, isunderrange, limit, dev, value)

"""
    markoverrange(batch::Batch; data = batch.method.est_conc, limit = uloq, dev = 0.15, value = (x, y) -> NaN) -> AbstractDataTable

Mark all samples that have quantity `data` over acceptable range with `value` for each analyte. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable upper limit.
* `dev::Tuple`: allowed devience of limit.
* `value::Function`: 2-arg function returning marked value for each analyte. The first argument is calibrator; the second is sample data vector.
"""
markoverrange(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = uloq, dev = 0.15, value = (x, y) -> NaN) = 
    _markoutofrange(batch, analyte, data, isoverrange, limit, dev, value)

"""
    markoutofrange!(batch::Batch[, analyte]; data = batch.method.est_conc, limit = dynamic_range, dev = (0.2, 0.15), value = (x, y) -> (NaN, NaN)) -> Batch

Mark all samples that have quantity `data` out of acceptable range with `value` for each analyte, and returns the updated `batch`. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable range.
* `dev::Tuple`: allowed devience of range.
* `value::Function`: 2-arg functions returning marked values for each analyte. The first argument is calibrator; the second is sample data vector.
"""
function markoutofrange!(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = dynamic_range, dev = (0.2, 0.15), value = (x, y) -> (NaN, NaN)) 
    set!(batch.data, data, markoutofrange(batch, analyte; data, limit, dev, value))
    batch
end

"""
    markunderrange!(batch::Batch[, analyte]; data = batch.method.est_conc, limit = lloq, dev = 0.2, value = (x, y) -> NaN) -> Batch

Mark all samples that have quantity `data` under acceptable range with `value` for each analyte, and returns the updated `batch`. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable lower limit.
* `dev::Real`: allowed devience of limit.
* `value::Function`: 2-arg function returning marked value for each analyte. The first argument is calibrator; the second is sample data vector.
"""
function markunderrange!(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = dynamic_range, dev = 0.2, value = (x, y) -> NaN) 
    set!(batch.data, data, markunderrange(batch, analyte; data, limit, dev, value))
    batch
end

"""
    markoverrange!(batch::Batch[, analyte]; data = batch.method.est_conc, limit = uloq, dev = 0.15, value = (x, y) -> NaN) -> Batch

Mark all samples that have quantity `data` over acceptable range with `value` for each analyte, and returns the updated `batch`. 

* `analyte` specifys target analyte; otherwise, all analytes are included.
* `data::Symbol`: quantity to be evaluated.
* `limit::Function` takes each calibrator and returns acceptable lower limit.
* `dev::Real`: allowed devience of limit.
* `value::Function`: 2-arg function returning marked value for each analyte. The first argument is calibrator; the second is sample data vector.
"""
function markoverrange!(batch::Batch, analyte = nothing; data = batch.method.est_conc, limit = dynamic_range, dev = 0.15, value = (x, y) -> NaN) 
    set!(batch.data, data, markoverrange(batch, analyte; data, limit, dev, value))
    batch
end
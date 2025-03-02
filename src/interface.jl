Tables.istable(::Type{<: AbstractDataTable}) = true
Tables.rowaccess(::Type{<: AbstractDataTable}) = true
Tables.columnaccess(::Type{<: AbstractDataTable}) = true
Tables.columns(dt::AbstractDataTable) = Tables.columns(table(dt))
Tables.rows(dt::AbstractDataTable) = Tables.rows(table(dt))
Tables.getcolumn(dt::AbstractDataTable, property) = Tables.getcolumn(table(dt), property)
Tables.columnnames(dt::AbstractDataTable) = Tables.columnnames(table(dt))
Base.getproperty(dt::AbstractDataTable, property::Symbol) = Base.getproperty(table(dt), property)
Base.propertynames(dt::AbstractDataTable) = Base.propertynames(table(dt))
Base.eltype(dt::AbstractDataTable) = Base.eltype(table(dt))
Base.length(dt::AbstractDataTable) = size(table(dt), 1)
Base.iterate(dt::SampleDataTable, st = 1) = st > length(dt) ? nothing : (first(Base.iterate(table(dt), st)), st + 1)
Base.iterate(dt::AnalyteDataTable, st = 1) = st > length(dt) ? nothing : (first(Base.iterate(table(dt), st)), st + 1)
Base.getindex(dt::AbstractDataTable, x...) = Base.getindex(table(dt), x...)
Base.getindex(dt::SampleDataTable{A, S}, analyte::B, sample::T) where {A, S, B <: A, T <: S} = getanalyte(dt, analyte)[findsample(dt, sample)]
Base.getindex(dt::AnalyteDataTable{A, S}, analyte::B, sample::T) where {A, S, B <: A, T <: S} = getsample(dt, sample)[findanalyte(dt, analyte)]
Base.setindex!(dt::AbstractDataTable, value, keys...) = Base.setindex!(table(dt), value, keys...)
Base.setindex!(dt::SampleDataTable{A, S, N}, value::N, analyte::B, sample::T) where {A, S, B <: A, T <: S, N} = (setindex!(getanalyte(dt, analyte), value, findsample(dt, sample)); dt)
Base.setindex!(dt::SampleDataTable{A, S, N}, value, analyte::B, sample::T) where {A, S, B <: A, T <: S, N} = (setindex!(getanalyte(dt, analyte), convert(N, value), findsample(dt, sample)); dt)
Base.setindex!(dt::AnalyteDataTable{A, S, N}, value::N, analyte::B, sample::T) where {A, S, B <: A, T <: S, N} = (setindex!(getsample(dt, sample), value, findanalyte(dt, analyte)); dt)
Base.setindex!(dt::AnalyteDataTable{A, S, N}, value, analyte::B, sample::T) where {A, S, B <: A, T <: S, N} = (setindex!(getsample(dt, sample), convert(N, value), findanalyte(dt, analyte)); dt)
Base.copy(dt::SampleDataTable) = SampleDataTable(copy(analyteobj(dt)), copy(sampleobj(dt)), samplecol(dt), copy(table(dt)))
Base.copy(dt::AnalyteDataTable) = AnalyteDataTable(copy(analyteobj(dt)), copy(sampleobj(dt)), analytecol(dt), copy(table(dt)))

function Base.getindex(calibration::AbstractVector{<: AbstractCalibration{A}}, analyte::B) where {A, B <: A}
    id = findfirst(x -> first(x.analyte) == analyte, calibration)
    isnothing(id) && throw(ArgumentError("Analyte $analyte does not have a calibration curve."))
    calibration[id]
end
function Base.setindex!(calibration::AbstractVector{<: AbstractCalibration{A}}, v::S, analyte::B) where {A, B <: A, S <: AbstractCalibration{B}}
    id = findfirst(x -> first(x.analyte) == analyte, calibration)
    isnothing(id) && throw(ArgumentError("Analyte $analyte does not have a calibration curve."))
    calibration[id] = v
    calibration
end
Base.copy(cal::MultipleCalibration) = MultipleCalibration(cal.analyte, cal.type, cal.zero, cal.weight, cal.formula, cal.table, cal.model)
Base.copy(cal::SingleCalibration) = SingleCalibration(cal.analyte, cal.conc)

Base.isassigned(at::AnalysisTable, i::Symbol) = Base.isassigned(tables(at), i)
Dictionaries.isinsertable(at::AnalysisTable) = true
Dictionaries.issetable(at::AnalysisTable) = true
Dictionaries.set!(at::AnalysisTable, i, v) = (set!(tables(at), i, v); at)
Dictionaries.unset!(at::AnalysisTable, i) = (unset!(tables(at), i); at)
Base.insert!(at::AnalysisTable, i, v) = (insert!(tables(at), i, v); at)
Base.get!(at::AnalysisTable, i, v) = get!(tables(at), i, v)
Base.get(at::AnalysisTable, i, v) = get(tables(at), i, v)
Base.delete!(at::AnalysisTable, i) = (delete!(tables(at), i); at)
Base.getproperty(at::AnalysisTable, property::Symbol) = tables(at)[property]
Base.propertynames(at::AnalysisTable) = Tuple(keys(tables(at)))
Base.getindex(at::AnalysisTable, i) = getindex(tables(at), i)
Base.setindex!(at::AnalysisTable, v, i) = (setindex!(tables(at), v, i); at)
Base.pairs(at::AnalysisTable) = pairs(tables(at))
Base.keys(at::AnalysisTable) = keys(tables(at))
Base.values(at::AnalysisTable) = values(tables(at))
Base.haskey(at::AnalysisTable, i) = haskey(tables(at), i)
Base.iterate(at::AnalysisTable, s...) = Base.iterate(tables(at), s...)
Base.length(at::AnalysisTable) = length(tables(at))
Base.copy(at::AnalysisTable) = AnalysisTable(copy(analyteobj(at)), copy(sampleobj(at)), Dictionary(collect(keys(tables(at))), values(tables(at))))
function Base.getproperty(tbl::AnalysisMethod, p::Symbol)
    if p == :analyte
        getfield(tbl, :analytetable).analyte
    elseif p == :isd
        getfield(tbl, :analytetable).analyte[getfield(tbl, :analytetable).isd .< 0]
    elseif p == :nonisd
        getfield(tbl, :analytetable).analyte[getfield(tbl, :analytetable).isd .>= 0]
    elseif p == :point
        s = getfield(tbl, :signaltable)
        isnothing(s) ? s : sampleobj(s)
    elseif p == :level
        sampleobj(getfield(tbl, :conctable))
    else
        getfield(tbl, p)
    end
end
Base.propertynames(tbl::AnalysisMethod) = (:analytetable, :signal, :rel_sig, :est_conc, :true_conc, :acc, :pointlevel, :conctable, :signaltable, :analyte, :isd, :nonisd, :point, :level)

function Base.getproperty(batch::Batch, p::Symbol)
    if p == :analyte
        getfield(batch, :method).analyte
    elseif p == :isd
        getfield(batch, :method).isd
    elseif p == :nonisd
        getfield(batch, :method).nonisd
    elseif p == :point
        getfield(batch, :method).point
    elseif p == :level
        getfield(batch, :method).level
    else
        getfield(batch, p)
    end
end
Base.propertynames(batch::Batch) = (:method, :calibration, :data, :analyte, :isd, :nonisd, :point, :level)

struct EachAnalyte{T}
    table::T
end

struct EachSample{T}
    table::T
end

"""
    eachanalyte(dt::AbstractDataTable)

Create an iterator which gets data belonging to each analyte as a `Vector`. For `AnalyteDataTable`, new vectors are created; mutating these vector will not change the value in `dt`.
"""
eachanalyte(dt::AbstractDataTable) = EachAnalyte(dt)
"""
    eachsample(dt::AbstractDataTable)

Create an iterator which gets data belonging to each sample as a `Vector`. For `SampleDataTable`, new vectors are created; mutating these vector will not change the value in `dt`.
"""
eachsample(dt::AbstractDataTable) = EachSample(dt)

Base.eltype(it::EachAnalyte{<: AbstractDataTable{A, S, N}}) where {A, S, N} = AbstractVector{N}
Base.eltype(it::EachSample{<: AbstractDataTable{A, S, N}}) where {A, S, N} = AbstractVector{N}
Base.length(it::EachAnalyte) = Base.length(analyteobj(it.table))
Base.length(it::EachSample) = Base.length(sampleobj(it.table))
Base.iterate(it::EachAnalyte{<: SampleDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(table(it.table), analytename(it.table)[st]), st + 1)
Base.iterate(it::EachAnalyte{<: AnalyteDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(table(it.table), p)[st] for p in samplename(it.table)], st + 1)
Base.iterate(it::EachSample{<: SampleDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(table(it.table), p)[st] for p in analytename(it.table)], st + 1)
Base.iterate(it::EachSample{<: AnalyteDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(table(it.table), samplename(it.table)[st]), st + 1)

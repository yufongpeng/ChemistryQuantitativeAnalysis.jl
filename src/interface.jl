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
Base.iterate(dt::ColumnDataTable, st = 1) = st > length(dt) ? nothing : (first(Base.iterate(table(dt), st)), st + 1)
Base.iterate(dt::RowDataTable, st = 1) = st > length(dt) ? nothing : (first(Base.iterate(table(dt), st)), st + 1)
Base.getindex(dt::AbstractDataTable, x...) = Base.getindex(table(dt), x...)
Base.setindex!(dt::AbstractDataTable, value, keys...) = Base.setindex!(table(dt), value, keys...)

function Base.getindex(calibration::AbstractVector{<: AbstractCalibration{A}}, analyte::B) where {A, B <: A}
    id = findfirst(x -> first(x.analyte) == analyte, calibration)
    isnothing(id) && throw(ArgumentError("Analyte $analyte does not have a calibration curve."))
    calibration[id]
end
function Base.setindex!(calibration::AbstractVector{<: AbstractCalibration{A}}, v::S, analyte::B) where {A, B <: A, S <: AbstractCalibration{B}}
    id = findfirst(x -> first(x.analyte) == analyte, calibration)
    isnothing(id) && throw(ArgumentError("Analyte $analyte does not have a calibration curve."))
    calibration[id] = v
end

Dictionaries.set!(at::AnalysisTable, i, v) = set!(tables(at), i, v)
Dictionaries.unset!(at::AnalysisTable, i) = unset!(tables(at), i)
Base.insert!(at::AnalysisTable, i, v) = insert!(tables(at), i, v)
Base.get!(at::AnalysisTable, i, v) = get!(tables(at), i, v)
Base.get(at::AnalysisTable, i, v) = get(tables(at), i, v)
Base.delete!(at::AnalysisTable, i) = delete!(tables(at), i)
Base.getproperty(at::AnalysisTable, property::Symbol) = tables(at)[property]
Base.propertynames(at::AnalysisTable) = Tuple(keys(tables(at)))
Base.getindex(at::AnalysisTable, i) = getindex(tables(at), i)
Base.pairs(at::AnalysisTable) = pairs(tables(at))
Base.keys(at::AnalysisTable) = keys(tables(at))
Base.haskey(at::AnalysisTable, i) = haskey(tables(at), i)
Base.iterate(at::AnalysisTable, st = 1) = st > length(at) ? nothing : (first(Base.iterate(tables(at), st)), st + 1)
Base.length(at::AnalysisTable) = length(tables(at))

function Base.getproperty(tbl::MethodTable, p::Symbol)
    if p == :analyte
        getfield(tbl, :analytetable).analyte
    elseif p == :isd
        getfield(tbl, :analytetable).analyte[getfield(tbl, :analytetable).isd .< 0]
    elseif p == :nonisd
        getfield(tbl, :analytetable).analyte[getfield(tbl, :analytetable).isd .>= 0]
    elseif p == :point
        s = getfield(tbl, :signaltable)
        isnothing(s) ? s : s.sample
    elseif p == :level
        getfield(tbl, :conctable).sample
    else
        getfield(tbl, p)
    end
end
Base.propertynames(tbl::MethodTable) = (:analytetable, :signal, :pointlevel, :conctable, :signaltable, :analyte, :isd, :nonisd, :point, :level)

function Base.getproperty(batch::Batch, p::Symbol)
    if p == :analyte
        getfield(batch, :method).analyte
    elseif p == :isd
        getfield(batch, :method).analyte[getfield(batch, :method).isd .< 0]
    elseif p == :nonisd
        getfield(batch, :method).analyte[getfield(batch, :method).isd .>= 0]
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

Create an iterator which gets data belonging to each analyte as a `Vector`. For `RowDataTable`, new vectors are created; mutating these vector will not change the value in `dt`.
"""
eachanalyte(dt::AbstractDataTable) = EachAnalyte(dt)
"""
    eachsample(dt::AbstractDataTable)

Create an iterator which gets data belonging to each sample as a `Vector`. For `ColumnDataTable`, new vectors are created; mutating these vector will not change the value in `dt`.
"""
eachsample(dt::AbstractDataTable) = EachSample(dt)

Base.eltype(it::EachAnalyte) = Base.eltype(it.table)
Base.eltype(it::EachSample) = Base.eltype(it.table)
Base.length(it::EachAnalyte) = Base.length(analyteobj(it.table))
Base.length(it::EachSample) = Base.length(sampleobj(it.table))
Base.iterate(it::EachAnalyte{<: ColumnDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(table(it.table), analytename(it.table)[st]), st + 1)
Base.iterate(it::EachAnalyte{<: RowDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(table(it.table), p)[st] for p in samplename(it.table)], st + 1)
Base.iterate(it::EachSample{<: ColumnDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(table(it.table), p)[st] for p in analytename(it.table)], st + 1)
Base.iterate(it::EachSample{<: RowDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(table(it.table), samplename(it.table)[st]), st + 1)

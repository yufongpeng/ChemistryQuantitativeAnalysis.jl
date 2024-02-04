Tables.istable(::Type{<: AbstractDataTable}) = true
Tables.rowaccess(::Type{<: AbstractDataTable}) = true
Tables.columnaccess(::Type{<: AbstractDataTable}) = true

Tables.columns(dt::AbstractDataTable) = Tables.columns(dt.table)
Tables.rows(dt::AbstractDataTable) = Tables.rows(dt.table)
Base.eltype(dt::AbstractDataTable) = Base.eltype(dt.table)
Base.length(dt::AbstractDataTable) = size(dt.table, 1)
Base.iterate(dt::ColumnDataTable, st = 1) = st > length(dt) ? nothing : (first(Base.iterate(dt.table, st)), st + 1)
Base.iterate(dt::RowDataTable, st = 1) = st > length(dt) ? nothing : (first(Base.iterate(dt.table, st)), st + 1)
Base.getindex(dt::AbstractDataTable, x...) = Base.getindex(dt.table, x...)
Base.setindex!(dt::AbstractDataTable, value, keys...) = Base.setindex!(dt.table, value, keys...)
function Base.getindex(calibration::AbstractVector{<: AbstractCalibration{A}}, analyte::B) where {A, B <: A}
    id = findfirst(x -> first(x.analyte) == analyte, calibration)
    isnothing(id) && throw(ArgumentError("Analyte $analyte does not have a calibration curve."))
    calibration[id]
end
Dictionaries.set!(at::AnalysisTable, i, v) = set!(at.tables, i, v)

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
Base.length(it::EachAnalyte) = Base.length(it.table.analyte)
Base.length(it::EachSample) = Base.length(it.table.sample)
Base.iterate(it::EachAnalyte{<: ColumnDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(it.table.table, it.table.analytename[st]), st + 1)
Base.iterate(it::EachAnalyte{<: RowDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(it.table.table, p)[st] for p in it.table.samplename], st + 1)
Base.iterate(it::EachSample{<: ColumnDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(it.table.table, p)[st] for p in it.table.analytename], st + 1)
Base.iterate(it::EachSample{<: RowDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(it.table.table, it.table.samplename[st]), st + 1)

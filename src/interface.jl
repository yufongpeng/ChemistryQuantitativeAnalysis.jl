Tables.istable(::Type{<: AbstractDataTable}) = true
Tables.rowaccess(::Type{<: AbstractDataTable}) = true
Tables.columnaccess(::Type{<: AbstractDataTable}) = true

Tables.columns(tbl::AbstractDataTable) = Tables.columns(tbl.table)
Tables.rows(tbl::AbstractDataTable) = Tables.rows(tbl.table)
Base.eltype(tbl::AbstractDataTable) = Base.eltype(tbl.table)
Base.length(tbl::AbstractDataTable) = size(tbl.table, 1)
Base.iterate(tbl::ColumnDataTable, st = 1) = st > length(tbl) ? nothing : (first(Base.iterate(tbl.table, st)), st + 1)
Base.iterate(tbl::RowDataTable, st = 1) = st > length(tbl) ? nothing : (first(Base.iterate(tbl.table, st)), st + 1)
Base.getindex(tbl::AbstractDataTable, x...) = Base.getindex(tbl.table, x...)
Base.setindex!(tbl::AbstractDataTable, value, keys...) = Base.setindex!(tbl.table, value, keys...)
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

eachanalyte(tbl::AbstractDataTable) = EachAnalyte(tbl)
eachsample(tbl::AbstractDataTable) = EachSample(tbl)
Base.eltype(it::EachAnalyte) = Base.eltype(it.table)
Base.eltype(it::EachSample) = Base.eltype(it.table)
Base.length(it::EachAnalyte) = Base.length(it.table.analyte)
Base.length(it::EachSample) = Base.length(it.table.sample)
Base.iterate(it::EachAnalyte{<: ColumnDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(it.table.table, it.table.analytename[st]), st + 1)
Base.iterate(it::EachAnalyte{<: RowDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(it.table.table, p)[st] for p in it.table.samplename], st + 1)
Base.iterate(it::EachSample{<: ColumnDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(it.table.table, p)[st] for p in it.table.analytename], st + 1)
Base.iterate(it::EachSample{<: RowDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(it.table.table, it.table.samplename[st]), st + 1)

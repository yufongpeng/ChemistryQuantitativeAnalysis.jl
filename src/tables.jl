Tables.istable(::Type{RowDataTable}) = true
Tables.istable(::Type{ColumnDataTable}) = true
Tables.rowaccess(::Type{RowDataTable}) = true
Tables.rowaccess(::Type{ColumnDataTable}) = true
Tables.columnaccess(::Type{RowDataTable}) = true
Tables.columnaccess(::Type{ColumnDataTable}) = true

Tables.columns(tbl::ColumnDataTable) = Tables.columns(tbl.table)
Tables.columns(tbl::RowDataTable) = Tables.columns(tbl.table)

Tables.rows(tbl::ColumnDataTable) = Tables.rows(tbl.table)
Tables.rows(tbl::RowDataTable) = Tables.rows(tbl.table)
Base.eltype(tbl::ColumnDataTable) = Base.eltype(tbl.table)
Base.eltype(tbl::RowDataTable) = Base.eltype(tbl.table)
Base.length(tbl::ColumnDataTable) = size(tbl.table, 1)
Base.length(tbl::RowDataTable) = size(tbl.table, 1)
Base.iterate(tbl::ColumnDataTable, st = 1) = st > length(tbl) ? nothing : (first(Base.iterate(tbl.table, st)), st + 1)
Base.iterate(tbl::RowDataTable, st = 1) = st > length(tbl) ? nothing : (first(Base.iterate(tbl.table, st)), st + 1)

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
Base.iterate(it::EachAnalyte{<: RowDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(dt.table, p)[st] for p in dt.samplename], st + 1)
Base.iterate(it::EachSample{<: ColumnDataTable}, st = 1) = st > length(it) ? nothing : ([getproperty(dt.table, p)[st] for p in dt.analytename], st + 1)
Base.iterate(it::EachSample{<: RowDataTable}, st = 1) = st > length(it) ? nothing : (getproperty(it.table.table, it.table.samplename[st]), st + 1)

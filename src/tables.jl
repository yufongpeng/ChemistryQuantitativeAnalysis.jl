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
Base.length(tbl::ColumnDataTable) = Base.length(tbl.table)
Base.length(tbl::RowDataTable) = Base.length(tbl.table)
Base.iterate(tbl::ColumnDataTable, st = 1) = st > length(tbl) ? nothing : (first(Base.iterate(tbl.table, st)), st + 1)
Base.iterate(tbl::RowDataTable, st = 1) = st > length(tbl) ? nothing : (first(Base.iterate(tbl.table, st)), st + 1)

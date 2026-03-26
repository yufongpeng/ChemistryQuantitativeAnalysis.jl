"""
    SampleDataTable(tbl::AnalyteDataTable{A, S, N, T}, samplecol::Symbol, tablesink::TypeOrFn = T)
    SampleDataTable(samplecol::Symbol, tablesink::TypeOrFn, tbl::AnalyteDataTable)
    SampleDataTable(samplecol::Symbol, tbl::AnalyteDataTable)

Convert `AnalyteDataTable` to `SampleDataTable` with `samplecol` as the column name of sample.
"""
SampleDataTable(tbl::AnalyteDataTable{A, S, N, T}, samplecol::Symbol, tablesink::TypeOrFn = T) where {A, S, N, T} = 
    SampleDataTable(analyteobj(tbl), samplecol, tablesink((; (samplecol => sampleobj(tbl), (analytename(tbl) .=> eachanalyte(tbl))...)...)))
SampleDataTable(samplecol::Symbol, tablesink::TypeOrFn, tbl::AnalyteDataTable) = 
    SampleDataTable(tbl, samplecol, tablesink)
SampleDataTable(samplecol::Symbol, tbl::AnalyteDataTable) = 
    SampleDataTable(tbl, samplecol)

"""
    AnalyteDataTable(tbl::SampleDataTable{A, S, N, T}, analytecol::Symbol, tablesink::TypeOrFn = T)
    AnalyteDataTable(analytecol::Symbol, tablesink::TypeOrFn, tbl::SampleDataTable)
    AnalyteDataTable(analytecol::Symbol, tbl::SampleDataTable)

Convert `SampleDataTable` to `AnalyteDataTable` with `analytecol` as the column name of analyte.
"""
AnalyteDataTable(tbl::SampleDataTable{A, S, N, T}, analytecol::Symbol, tablesink::TypeOrFn = T) where {A, S, N, T} = 
    AnalyteDataTable(analytecol, sampleobj(tbl), tablesink((; (analytecol => analyteobj(tbl), (samplename(tbl) .=> eachsample(tbl))...)...)))
AnalyteDataTable(analytecol::Symbol, tablesink::TypeOrFn, tbl::SampleDataTable) = 
    AnalyteDataTable(tbl, analytecol, tablesink)
AnalyteDataTable(analytecol::Symbol, tbl::SampleDataTable) = 
    AnalyteDataTable(tbl, analytecol)
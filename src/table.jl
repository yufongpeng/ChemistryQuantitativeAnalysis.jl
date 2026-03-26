"""
    SampleDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, `N` determines numeric value type, and `T` determines table type.

# Fields
* `analyte::Vector{A}`: analytes in user-defined types. Use `analyteobj` to get this field.
* `sample::Vector{S}`: samples in user-defined types, i.e., `getproperty(table(dt), samplecol(dt))`. Use `sampleobj` to get this field.
* `samplecol::Symbol`: the column name that each element is sample name. Use `samplecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct SampleDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}
    analyte::Vector{A}
    sample::Vector{S}
    samplecol::Symbol
    table::T
    function SampleDataTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, samplecol::Symbol, table::T) where {A, S, T}
        allunique(sample) || throw("`sample` should have unique elemnts.")
        tp = collect(map(eltype, columns(table)))
        ps = propertynames(table)
        idp = findfirst(==(Symbol(first(analyte))), ps)
        isnothing(idp) && throw(ArgumentError("`analyte[1]` is not in the table"))
        N = tp[idp]
        for i in 2:lastindex(analyte)
            idp = findfirst(==(Symbol(analyte[i])), ps)
            isnothing(idp) && throw(ArgumentError("`analyte[$i]` is not in the table"))
            N == tp[idp] || throw(ArgumentError(string("The element type of `analyte[$i]` is `", tp[idp], "`; convert it to `", N, "`.")))
        end
        new{A, S, N, T}(convert(Vector{A}, analyte), convert(Vector{S}, sample), samplecol, table)
    end
end
"""
    SampleDataTable(analyte::AbstractVector, samplecol::Symbol, table)
    SampleDataTable(table, analytename::AbstractVector, samplecol::Symbol)

An interface equivalent to `SampleDataTable(analyte, getproperty(table, samplecol), samplecol, table)`.
"""
SampleDataTable(analyte::AbstractVector, samplecol::Symbol, table) = 
    SampleDataTable(analyte, getproperty(table, samplecol), samplecol, table)
SampleDataTable(table, analyte::AbstractVector, samplecol::Symbol) = 
    SampleDataTable(analyte, samplecol, table)
"""
    SampleDataTable(analytetype::TypeOrFn, samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))
    SampleDataTable(table, analytetype::TypeOrFn, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))
    SampleDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))
    SampleDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))
    

Higher level interfaces for `SampleDataTable{A}`.

* `analytename::AbstractVector`: the column names of `table` that are analyte names. It will be converted to `AbstractVector{String}` before conversion to `AbstractVector{A}` (`string.(analytename)`).
* `samplecol::Symbol`: the column name that each element is sample name.
* `analytetype::Union{A, Function}`: After first conversion of `analytename`, it convert `analytename` to `AbstractVector{A}` using `cqamap(analytetype, analytename)`. See `cqaconvert` for the requirement of `analytetype`.
"""
SampleDataTable(analytetype::TypeOrFn, samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(cqamap(analytetype, string.(analytename)), samplecol, table)
SampleDataTable(table, analytetype::TypeOrFn, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(analytetype, samplecol, table; analytename)
SampleDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(string.(analytename), samplecol, table)
SampleDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(samplecol, table; analytename)
"""
    AnalyteDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}

Tabular data wrapper indicates part of columns represent analyte, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, `N` determines numeric value type, and `T` determines table type.

# Fields
* `analyte::Vector{A}`: analytes in user-defined types, i.e., `getproperty(table(dt), analytecol(dt))`. Use `analyteobj` to get this field.
* `sample::Vector{S}`: samples in user-defined types. Use `sampleobj` to get this field.
* `analytecol::Symbol`: the column name that each element is analyte name. Use `analytecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct AnalyteDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}
    analyte::Vector{A}
    sample::Vector{S}
    analytecol::Symbol
    table::T
    function AnalyteDataTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, analytecol::Symbol, table::T) where {A, S, T}
        allunique(analyte) || throw("`analyte` should have unique elemnts.")
        tp = collect(map(eltype, columns(table)))
        ps = propertynames(table)
        idp = findfirst(==(Symbol(first(sample))), ps)
        isnothing(idp) && throw(ArgumentError("`sample[1]` is not in the table"))
        N = tp[idp]
        for i in 2:lastindex(sample)
            idp = findfirst(==(Symbol(sample[i])), ps)
            isnothing(idp) && throw(ArgumentError("`sample[$i]` is not in the table"))
            N == tp[idp] || throw(ArgumentError(string("The element type of `sample[$i]` is `", tp[idp], "`; convert it to `", N, "`.")))
        end
        new{A, S, N, T}(convert(Vector{A}, analyte), convert(Vector{S}, sample), analytecol, table)
    end
end
"""
    AnalyteDataTable(analytecol::Symbol, samplename::AbstractVector, table)
    AnalyteDataTable(table, analytecol::Symbol, samplename::AbstractVector)
    
An interface equivalent to `AnalyteDataTable(getproperty(table, analytecol), samplename, analytecol, table)`.
"""
AnalyteDataTable(analytecol::Symbol, samplename::AbstractVector, table) = 
    AnalyteDataTable(getproperty(table, analytecol), samplename, analytecol, table)
AnalyteDataTable(table, analytecol::Symbol, samplename::AbstractVector) = 
    AnalyteDataTable(analytecol, samplename, table)
"""
    AnalyteDataTable(analytecol::Symbol, sampletype::TypeOrFn, table; samplename = setdiff(propertynames(table), [analytecol]))
    AnalyteDataTable(table, analytecol::Symbol, sampletype::TypeOrFn; samplename = setdiff(propertynames(table), [analytecol]))
    AnalyteDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol]))
    AnalyteDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol]))

Higher level interfaces for `AnalyteDataTable{A, S}`.

* `analytecol::Symbol`: the column name that each element is analyte name.
* `samplename::AbstractVector`: the column names of `table` that are sample names. It will be converted to `Vector{String}` before conversion to `AbstractVector{S}` (`string.(samplename)`).
* `sampletype::Union{S, Function}`: After first conversion of `sampletype`, it converts `samplename` to `AbstractVector{S}` using `cqamap(sampletype, sampletype)`. See `cqaconvert` for the requirement of `sampletype`.
"""
AnalyteDataTable(analytecol::Symbol, sampletype::TypeOrFn, table; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, cqamap(sampletype, string.(samplename)), table)
AnalyteDataTable(table, analytecol::Symbol, sampletype::TypeOrFn; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, sampletype, table; samplename)
AnalyteDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, string.(samplename), table)
AnalyteDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, table; samplename)

"""
    AnalysisTable{A, S, T <: AbstractDataTable{A, S}}

Wrapper of multiple tables representing different types of values. `A` determines analyte type, `S` determines sample type, and `T` determines datatable type. It is implemented as a `Dictionary{Symbol, T <: AbstractDataTable{A, S}}`, but unlike dictionaries, each datatable can also be extracted using `getproperty`. 

# Fields
* `analyte::Vector{A}`: analytes in user-defined types. Use `analyteobj` to get this field.
* `sample::Vector{S}`: samples in user-defined types. Use `sampleobj` to get this field.
* `tables::Dictionary{Symbol, T}`: a dictionary mapping data type to datatable. Use `tables` to get this field.

# Properties
All keys of `tables`.
"""
struct AnalysisTable{A, S, T <: AbstractDataTable{A, S}}
    analyte::Vector{A}
    sample::Vector{S}
    tables::Dictionary{Symbol, T}
    AnalysisTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, tables::Dictionary{Symbol, T}) where {A, S, T <: AbstractDataTable{A, S}} = new{A, S, T}(convert(Vector{A}, analyte), convert(Vector{S}, sample), tables)
end

"""
    AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S}})

A `Dictionary`-like constructor for `AnalysisTable` from two iterable inputs `keys` and `tables`. The first value of `keys` will be the index for the first value of `tables`.
"""
function AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S}}) where {A, S}
    allequal(analyteobj(table) for table in tables) || throw(ArgumentError("Tables should have identical analyte"))
    allequal(sampleobj(table) for table in tables) || throw(ArgumentError("Tables should have identical sample"))
    AnalysisTable(deepcopy(analyteobj(first(tables))), deepcopy(sampleobj(first(tables))), Dictionary(keys, tables))
end
"""
    analysistable(iter)

A `dictionary`-like function for constructing `AnalysisTable` from an iterable iter of key-value `Pair`s (or other iterables of two elements, such as a two-tuples). Keys should be `Symbol`s, and values should be `AbstractDataTable`s.
"""
function analysistable(iter)
    allequal(analyteobj(last(kv)) for kv in iter) || throw(ArgumentError("Tables should have identical analyte"))
    allequal(sampleobj(last(kv)) for kv in iter) || throw(ArgumentError("Tables should have identical sample"))
    v = last(first(iter))
    AnalysisTable(deepcopy(analyteobj(v)), deepcopy(sampleobj(v)), dictionary(iter))
end
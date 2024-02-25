module ChemistryQuantitativeAnalysis

using GLM, CSV, TypedTables, LinearAlgebra, Dictionaries, ThreadsX, Tables
export ColumnDataTable, RowDataTable, 
    AnalysisTable, analysistable, AnalysisMethod, Batch, 
    MultipleCalibration, SingleCalibration, calibration, update_calibration!,
    analyteobj, sampleobj, analytename, samplename, analytecol, samplecol,
    inv_predict, set_inv_predict, set_inv_predict!, update_inv_predict!,
    relative_signal, set_relative_signal, set_relative_signal!, update_relative_signal!,
    quantification, set_quantification, set_quantification!, update_quantification!,
    accuracy, set_accuracy, set_accuracy!, update_accuracy!,
    isdof, isisd, 
    findanalyte, getanalyte, findsample, getsample, eachanalyte, eachsample,
    dynamic_range, lloq, uloq, signal_range, signal_lloq, signal_uloq, 
    formula_repr, weight_repr, weight_value, formula_repr_utf8, weight_repr_utf8, format_number, mkbatch, 
    typedmap

import Base: getproperty, propertynames, show, write, eltype, length, iterate, 
        getindex, setindex!, insert!, get!, delete!, get, 
        pairs, keys, values, haskey, isassigned, 
        copy
import Dictionaries: set!, unset!, isinsertable, issettable
import Tables: istable, rowaccess, rows, columnaccess, columns
    
abstract type AbstractCalibration{A, N} end
abstract type AbstractDataTable{A, S, N, T} end
# abstract type AbstractAnalysisTable{A, T, S} end
const TypeOrFn = Union{<: Type, <: Function}

"""
    ColumnDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, `N` determines numeric value type, and `T` determines table type.

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types. Use `analyteobj` to get this field.
* `sample`: `Vector{S}`, samples in user-defined types, i.e., `getproperty(table(dt), samplecol(dt))`. Use `sampleobj` to get this field.
* `samplecol`: `Symbol`, the column name that each element is sample name. Use `samplecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct ColumnDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}
    analyte::Vector{A}
    sample::Vector{S}
    samplecol::Symbol
    table::T
    function ColumnDataTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, samplecol::Symbol, table::T) where {A, S, T}
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
    ColumnDataTable(analyte::AbstractVector, samplecol::Symbol, table)
    ColumnDataTable(table, analytename::AbstractVector, samplecol::Symbol)

An interface equivalent to `ColumnDataTable(analyte, getproperty(table, samplecol), samplecol, table)`.
"""
ColumnDataTable(analyte::AbstractVector, samplecol::Symbol, table) = 
    ColumnDataTable(analyte, getproperty(table, samplecol), samplecol, table)
ColumnDataTable(table, analyte::AbstractVector, samplecol::Symbol) = 
    ColumnDataTable(analyte, samplecol, table)
"""
    ColumnDataTable(analytetype::TypeOrFn, samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))
    ColumnDataTable(table, analytetype::TypeOrFn, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))
    ColumnDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))
    ColumnDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))
    

Higher level interfaces for `ColumnDataTable{A}`.

* `analytename`: `AbstractVector`, the column names of `table` that are analyte names. It will be converted to `AbstractVector{String}` before conversion to `AbstractVector{A}` (`string.(analytename)`).
* `samplecol`: `Symbol`, the column name that each element is sample name.
* `analytetype`: type `A` or a function. After first conversion of `analytename`, it convert `analytename` to `AbstractVector{A}` using `cqamap(analytetype, analytename)`. See `cqaconvert` for the requirement of `analytetype`.
"""
ColumnDataTable(analytetype::TypeOrFn, samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol])) = 
    ColumnDataTable(cqamap(analytetype, string.(analytename)), samplecol, table)
ColumnDataTable(table, analytetype::TypeOrFn, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol])) = 
    ColumnDataTable(analytetype, samplecol, table; analytename)
ColumnDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol])) = 
    ColumnDataTable(string.(analytename), samplecol, table)
ColumnDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol])) = 
    ColumnDataTable(samplecol, table; analytename)
"""
    RowDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}

Tabular data wrapper indicates part of columns represent analyte, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, `N` determines numeric value type, and `T` determines table type.

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types, i.e., `getproperty(table(dt), analytecol(dt))`. Use `analyteobj` to get this field.
* `sample`: `Vector{S}`, samples in user-defined types. Use `sampleobj` to get this field.
* `analytecol`: `Symbol`, the column name that each element is analyte name. Use `analytecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct RowDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}
    analyte::Vector{A}
    sample::Vector{S}
    analytecol::Symbol
    table::T
    function RowDataTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, analytecol::Symbol, table::T) where {A, S, T}
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
    RowDataTable(analytecol::Symbol, samplename::AbstractVector, table)
    RowDataTable(table, analytecol::Symbol, samplename::AbstractVector)
    
An interface equivalent to `RowDataTable(getproperty(table, analytecol), samplename, analytecol, table)`.
"""
RowDataTable(analytecol::Symbol, samplename::AbstractVector, table) = 
    RowDataTable(getproperty(table, analytecol), samplename, analytecol, table)
RowDataTable(table, analytecol::Symbol, samplename::AbstractVector) = 
    RowDataTable(analytecol, samplename, table)
"""
    RowDataTable(analytecol::Symbol, sampletype::TypeOrFn, table; samplename = setdiff(propertynames(table), [analytecol]))
    RowDataTable(table, analytecol::Symbol, sampletype::TypeOrFn; samplename = setdiff(propertynames(table), [analytecol]))
    RowDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol]))
    RowDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol]))

Higher level interfaces for `RowDataTable{A, S}`.

* `analytecol`: `Symbol`, the column name that each element is analyte name.
* `samplename`: `AbstractVector`, the column names of `table` that are sample names. It will be converted to `Vector{String}` before conversion to `AbstractVector{S}` (`string.(samplename)`).
* `sampletype`: type `S` or a function. After first conversion of `sampletype`, it converts `samplename` to `AbstractVector{S}` using `cqamap(sampletype, sampletype)`. See `cqaconvert` for the requirement of `sampletype`.
"""
RowDataTable(analytecol::Symbol, sampletype::TypeOrFn, table; samplename = setdiff(propertynames(table), [analytecol])) = 
    RowDataTable(analytecol, cqamap(sampletype, string.(samplename)), table)
RowDataTable(table, analytecol::Symbol, sampletype::TypeOrFn; samplename = setdiff(propertynames(table), [analytecol])) = 
    RowDataTable(analytecol, sampletype, table; samplename)
RowDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol])) = 
    RowDataTable(analytecol, string.(samplename), table)
RowDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol])) = 
    RowDataTable(analytecol, table; samplename)

"""
    AnalysisTable{A, S, T <: AbstractDataTable{A, S}}

Wrapper of multiple tables representing different types of values. `A` determines analyte type, `S` determines sample type, and `T` determines datatable type. It is implemented as a `Dictionary{Symbol, T <: AbstractDataTable{A, S}}`, but unlike dictionaries, each datatable can also be extracted using `getproperty`. 

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types. Use `analyteobj` to get this field.
* `sample`: `Vector{S}`, samples in user-defined types. Use `sampleobj` to get this field.
* `tables`: `Dictionary{Symbol, T}`, a dictionary mapping data type to datatable. Use `tables` to get this field.

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

"""
    AnalysisMethod{A, M <: Table, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}

A type containing analytes settings, and source calibration data. `A` determines analyte type.

# Fields
* `analytetable`: `M` contaning three columns.
    * `analyte`: `AbstractVector{A}`, analytes in user-defined types.
    * `isd`: `AbstractVector{Int}`, index of internal standard. `0` means no internal standard, and `-1` means the analyte itself is a internal standard.
    * `calibration`: index of analyte for calibration curve. `-1` means the analyte itself is a internal standard, so it will not be put into any calibration curve.
* `signal`: `Symbol`, data type for quantification, e.g. `:area`.
* `pointlevel`: `Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.
* `conctable`: `C <: AbstractDataTable{<: A, Int}` containing concentration data for each level. Samples must integers. One level indicates using `SingleCalibration`.
* `signaltable`: `D <: AbstractDataTable{A}` containig signal for each point. It can be `nothing` if signal data is unecessary.

# Properties
* `analyte`: `AbstractVector{A}`, analytes in user-defined types, identical to `analytetable.analyte`.
* `isd`: `AbstractVector{A}`, analytes which are internal standards.
* `nonisd`: `AbstractVector{A}`, analytes which are not internal standards.
* `point`: `AbstractVector{S}`, calibration points, identical to `sampleobj(signaltable)`. If `signaltable` is `nothing`, this value is `nothing` too.
* `level`: `AbstractVector{Int}`, calibration levels, identical to `sampleobj(conctable)`.
"""
struct AnalysisMethod{A, M <: Table, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}
    analytetable::M
    signal::Symbol
    pointlevel::Vector{Int}
    conctable::C
    signaltable::D
    function AnalysisMethod(
                    analytetable::M,
                    signal::Symbol,
                    pointlevel::AbstractVector{Int}, 
                    conctable::AbstractDataTable{B, Int}, 
                    signaltable::AbstractDataTable{C, S}
                ) where {B, M <: Table, C, S}
        cp = propertynames(analytetable)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analytetable"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analytetable"))
        :calibration in cp || throw(ArgumentError("Column `:calibration` is required in `analytetable"))
        A = eltype(analytetable.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable $B should be a subtype of $A"))
        C <: A || throw(ArgumentError("Analyte type of signaltable $C should be a subtype of $A"))
        for a in analyteobj(conctable)
            a in analytetable.analyte || throw(ArgumentError("Analyte `$a` is not in the `analytetable`."))
        end
        if length(sampleobj(conctable)) > 1
            length(pointlevel) == length(sampleobj(signaltable)) || throw(ArgumentError("The length of `pointlevel` is different from that of `sampleobj(table)`."))
            for a in analyteobj(conctable)
                a in analyteobj(signaltable) || throw(ArgumentError("Analyte `$a` is not in the `signatable`."))
            end
        end
        new{A, M, typeof(conctable), typeof(signaltable)}(analytetable, signal, convert(Vector{Int}, pointlevel), conctable, signaltable)
    end
    function AnalysisMethod(
                    analytetable::M,
                    signal::Symbol,
                    pointlevel::AbstractVector{Int}, 
                    conctable::AbstractDataTable{B, Int}, 
                    signaltable::Nothing
                ) where {B, M <: Table}
        cp = propertynames(analytetable)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analytetable"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analytetable"))
        :calibration in cp || throw(ArgumentError("Column `:calibration` is required in `analytetable"))
        A = eltype(analytetable.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable should be a subtype of $A"))
        for a in analyteobj(conctable)
            a in analytetable.analyte || throw(ArgumentError("Analyte `$a` is not in the `analytetable`."))
        end
        new{A, M, typeof(conctable), Nothing}(analytetable, signal, convert(Vector{Int}, pointlevel), conctable, signaltable)
    end
end

"""
    AnalysisMethod(conctable::AbstractDataTable{A, Int}, signaltable::ColumnDataTable, signal::Symbol, levelname::Symbol; kwargs...)
    AnalysisMethod(conctable::AbstractDataTable{A, Int}, signaltable::Nothing, signal::Symbol; kwargs...)
    AnalysisMethod(conctable::AbstractDataTable, signaltable::Union{AbstractDataTable, Nothing}, signal::Symbol, pointlevel::AbstractVector{Int}; kwargs...)

Convenient contructors for `AnalysisMethod` which make constructing a `Table` optional. `kwargs` will be columns in `analytetable`; when `analyte`, `isd` and `calibration` are not provided, it will use analyte in `conctable`. 
    
`levelname` is the column name for `pointlevel`. See `AnalysisMethod` for details of other arguments.
"""
AnalysisMethod(conctable::AbstractDataTable{A, Int}, signaltable::ColumnDataTable, signal::Symbol, levelname::Symbol; kwargs...) where A = 
    AnalysisMethod(conctable, signaltable, signal, getproperty(signaltable, levelname)::AbstractVector{Int}; kwargs...)
AnalysisMethod(conctable::AbstractDataTable{A, Int}, signaltable::Nothing, signal::Symbol; kwargs...) where A = AnalysisMethod(conctable, signaltable, signal, Int[]; kwargs...)
function AnalysisMethod(conctable::AbstractDataTable, 
                    signaltable::Union{AbstractDataTable, Nothing},
                    signal::Symbol,
                    pointlevel::AbstractVector{Int}; kwargs...
            )
    analyte = get(kwargs, :analyte, analyteobj(conctable))
    isd = get(kwargs, :isd, zeros(Int, length(analyteobj(conctable))))
    calibration = get(kwargs, :calibration, collect(eachindex(analyteobj(conctable))))
    analytetable = Table(; analyte, isd, calibration, kwargs...)    
    AnalysisMethod(analytetable, signal, pointlevel, conctable, signaltable)
end
@deprecate MethodTable AnalysisMethod

"""
    MultipleCalibration{A, N, T <: Table} <: AbstractCalibration{A, N}

A mutable type holding all data and settings for a calibration curve. `A` determines analyte type, and `N` determines numeric type. Only `Float64` is supported because of issue of `GLM`.

# Fields
* `analyte`: `Tuple{A, Any}`. First element is the analyte being quantified, and the second element is its internal standard for which `nothing` indicates no internal standard.
* `type`: `Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`).
* `zero`: `Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`).
* `weight`: `Float64` represents the exponential applying to each element of `x` as a weighting vector.
* `formula`: `FormulaTerm`, the formula for fitting calibration curve.
* `table`: `TypedTable.Table`, the cleaned up calibration data, containing 7 columns.
    * `id`: Point name
    * `level`: The index of concentration level. The larger, the higher concentraion it represents.
    * `y`: Signal or relative signal
    * `x`: True concentraion
    * `x̂`: Predicted concentration
    * `accuracy`: Accuracy, i.e. `x̂/x`.
    * `include`: Whether this point is included or not
* `model`: `GLM` object.
"""
mutable struct MultipleCalibration{A, N, T <: Table} <: AbstractCalibration{A, N}
    analyte::Tuple{A, Any}
    type::Bool
    zero::Bool
    weight::N
    formula::FormulaTerm
    table::T
    model
    function MultipleCalibration(analyte::Tuple{A, B}, type, zero, weight, formula, table::T, model) where {A, T, B}
        N = eltype(table.x)
        N == eltype(table.y) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column y is `", eltype(table.y), "`.")))
        N == eltype(table.x̂) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column x̂ is `", eltype(table.x̂), "`.")))
        N == eltype(table.accuracy) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column accuracy is `", eltype(table.accuracy), "`.")))
        new{A, N, T}(analyte, type, zero, N(weight), formula, table, model)
    end
end

"""
    SingleCalibration{A, N} <: AbstractCalibration{A, N}

A mutable type holding all data for single point calibration.

# Fields
* `analyte`: `Tuple{A}` is the analyte with known concentration (internal standard).
* `conc`: `Float64`, concentration of analyte.
"""
mutable struct SingleCalibration{A, N} <: AbstractCalibration{A, N}
    analyte::Tuple{A}
    conc::N
end

"""
    Batch{A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}, D <: Union{AnalysisTable{<: A}, Nothing}}

A type representing a batch for quantitative analysis. `A` determines analyte type, and `T` determines table type.

# Fields
* `method`: `M`, method.
* `calibration`: `C`, calibration curves.
* `data`: `D`, data for analysis, .

# Properties
* `analyte`: `AbstractVector{A}`, analytes in user-defined types, identical to `method.analyte`.
* `isd`: `AbstractVector{<: A}`, analytes which are internal standards, identical to `method.isd`.
* `nonisd`: `AbstractVector{<: A}`, analytes which are not internal standards, identical to `method.nonisd`.
* `point`: `AbstractVector` or `Nothing`, calibration points, identical to `method.point`.
* `level`: `AbstractVector{Int}`, calibration levels, identical to `method.level`.

# Constructors
* `Batch(method, calibration, data = nothing)`
"""
struct Batch{A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}, D <: Union{AnalysisTable{<: A}, Nothing}}
    method::M
    calibration::C
    data::D
    Batch(method::M, calibration::C, data::D) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}, D <: AnalysisTable{<: A}} = new{A, M, C, D}(method, calibration, data)
    Batch(method::M, calibration::C, data::Nothing) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}} = new{A, M, C, Nothing}(method, calibration, data)
    Batch(method::M, calibration::C) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}} = new{A, M, C, Nothing}(method, calibration, nothing)
end

"""
    Batch(method::AnalysisMethod, data = nothing; type = true, zero = false, weight = 0)

Construct a `Batch` from `method`, and optionally `data` with specified calibration parameters. See "MultipleCalibration" for detail description of keyword arguments.
"""
function Batch(method::AnalysisMethod{A}, data = nothing;
                type = true, 
                zero = false, 
                weight = 0
            ) where A
    Batch(
        method,
        if length(sampleobj(method.conctable)) > 1 
            map(analyteobj(method.conctable)) do analyte
                calibration(method, analyte; type, zero, weight)
            end
        else
            map(analyteobj(method.conctable)) do analyte
                SingleCalibration((analyte, ), first(getanalyte(method.conctable, analyte)))
            end
        end
        ,
        data
    )
end
"""
    Batch(batch::Batch, at::AnalysisTable)

Construct a new batch from an old batch and a `AnalysisTable` as new data.
"""
Batch(batch::Batch, at::AnalysisTable) = Batch(batch.method, batch.calibration, at)

include("interface.jl")
include("utils.jl")
include("cal.jl")
include("io.jl")

"""
    ColumnDataTable(tbl::RowDataTable{A, S, N, T}, samplecol::Symbol, tablesink::TypeOrFn = T)
    ColumnDataTable(samplecol::Symbol, tablesink::TypeOrFn, tbl::RowDataTable)
    ColumnDataTable(samplecol::Symbol, tbl::RowDataTable)

Convert `RowDataTable` to `ColumnDataTable` with `samplecol` as the column name of sample.
"""
ColumnDataTable(tbl::RowDataTable{A, S, N, T}, samplecol::Symbol, tablesink::TypeOrFn = T) where {A, S, N, T} = 
    ColumnDataTable(analyteobj(tbl), samplecol, tablesink((; (samplecol => sampleobj(tbl), (analytename(tbl) .=> eachanalyte(tbl))...)...)))
ColumnDataTable(samplecol::Symbol, tablesink::TypeOrFn, tbl::RowDataTable) = 
    ColumnDataTable(tbl, samplecol, tablesink)
ColumnDataTable(samplecol::Symbol, tbl::RowDataTable) = 
    ColumnDataTable(tbl, samplecol)

"""
    RowDataTable(tbl::ColumnDataTable{A, S, N, T}, analytecol::Symbol, tablesink::TypeOrFn = T)
    RowDataTable(analytecol::Symbol, tablesink::TypeOrFn, tbl::ColumnDataTable)
    RowDataTable(analytecol::Symbol, tbl::ColumnDataTable)

Convert `ColumnDataTable` to `RowDataTable` with `analytecol` as the column name of analyte.
"""
RowDataTable(tbl::ColumnDataTable{A, S, N, T}, analytecol::Symbol, tablesink::TypeOrFn = T) where {A, S, N, T} = 
    RowDataTable(analytecol, sampleobj(tbl), tablesink((; (analytecol => analyteobj(tbl), (samplename(tbl) .=> eachsample(tbl))...)...)))
RowDataTable(analytecol::Symbol, tablesink::TypeOrFn, tbl::ColumnDataTable) = 
    RowDataTable(tbl, analytecol, tablesink)
RowDataTable(analytecol::Symbol, tbl::ColumnDataTable) = 
    RowDataTable(tbl, analytecol)

end
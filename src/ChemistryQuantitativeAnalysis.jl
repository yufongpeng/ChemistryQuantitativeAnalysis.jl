module ChemistryQuantitativeAnalysis

using GLM, CSV, TypedTables, LinearAlgebra, Dictionaries, ThreadsX, Tables
import Tables: istable, rowaccess, rows, columnaccess, columns
export MultipleCalibration, SingleCalibration, 
    ColumnDataTable, RowDataTable, AnalysisTable, MethodTable,
    analyteobj, sampleobj, analytename, samplename, analytecol, samplecol,
    Batch, calibration, update_calibration!,
    dynamic_range, lloq, uloq, signal_range, signal_lloq, signal_uloq, accuracy, accuracy!, set_accuracy, set_accuracy!, update_accuracy!,
    inv_predict, inv_predict!, inv_predict_accuracy!, set_inv_predict, set_inv_predict!, update_inv_predict!,
    relative_signal, set_relative_signal, set_relative_signal!, update_relative_signal!,
    quantification, quantification!, set_quantification, set_quantification!, update_quantification!,
    findanalyte, getanalyte, findsample, getsample, set_isd!, eachanalyte, eachsample,
    formula_repr, weight_repr, weight_value, formula_repr_utf8, weight_repr_utf8, format_number, mkbatch, typedmap

import Base: getproperty, propertynames, show, write, eltype, length, iterate, getindex, setindex!, insert!, get!, delete!, get, pairs, keys, haskey
import Dictionaries: set!, unset!
    
abstract type AbstractCalibration{A} end
abstract type AbstractDataTable{A, S, T, V, U} end
# abstract type AbstractAnalysisTable{A, T, S} end

"""
    ColumnDataTable{A, S, T, V <: AbstractVector{A}, U <: AbstractVector{S}} <: AbstractDataTable{A, S, T, V, U}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, and `T` determines table type.

# Fields
* `analyte`: `V <: AbstractVector{A}`, analytes in user-defined types. Use `analyteobj` to get this field.
* `sample`: `U <: AbstractVector{S}`, samples in user-defined types, i.e., `getproperty(table(dt), samplecol(dt))`. Use `sampleobj` to get this field.
* `samplecol`: `Symbol`, the column name that each element is sample name. Use `samplecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct ColumnDataTable{A, S, T, V <: AbstractVector{A}, U <: AbstractVector{S}} <: AbstractDataTable{A, S, T, V, U}
    analyte::V
    sample::U
    samplecol::Symbol
    table::T
    function ColumnDataTable(analyte::V, sample::U, samplecol::Symbol, table::T) where {T, V <: AbstractVector, U <: AbstractVector}
        sample === getproperty(table, samplecol) || throw(ArgumentError("`sample` is not identical to `getproperty(table, sample)`."))
        new{eltype(analyte), eltype(sample), T, V, U}(analyte, sample, samplecol, table)
    end
end

"""
    ColumnDataTable(analytename::AbstractVector{Symbol}, samplecol::Symbol, table; analytetype = String)
    ColumnDataTable(analytename::AbstractVector{B}, samplecol::Symbol, table; analytetype = String)
    ColumnDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]), analytetype = String)
    ColumnDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]), analytetype = String)

User-friendly contructors for `ColumnDataTable{A}`.

* `analytename`: `AbstractVector{Symbol}`, the column names of `table` that are analyte names. It will be converted to `Vector{String}` before conversion to `AbstractVector{A}`. Vector of other type is also valid and will not be converted to `Vector{String}` first.
* `samplecol`: `Symbol`, the column name that each element is sample name.
* `analytetype`: type `A` or a function. After first conversion of `analytename`, it will be converted to `AbstractVector{A}` using `cqamap(analytetype, analytename)`. See `cqaconvert` for the requirement of `analytetype`.
"""
ColumnDataTable(analytename::AbstractVector{Symbol}, samplecol::Symbol, table; analytetype = String) = 
    ColumnDataTable(cqamap(analytetype, string.(analytename)), getproperty(table, samplecol), samplecol, table)
ColumnDataTable(analytename::AbstractVector, samplecol::Symbol, table; analytetype = String) = 
    ColumnDataTable(cqamap(analytetype, analytename), getproperty(table, samplecol), samplecol, table)
ColumnDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]), analytetype = String) = 
    ColumnDataTable(analytename, samplecol, table; analytetype)
ColumnDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]), analytetype = Strin) = 
    ColumnDataTable(analytename, samplecol, table; analytetype)

"""
    RowDataTable{A, S, T, V <: AbstractVector{A}, U <: AbstractVector{S}} <: AbstractDataTable{A, S, T, V, U}

Tabular data wrapper indicates part of columns represent analyte, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, and `T` determines table type.

# Fields
* `analyte`: `V <: AbstractVector{A}`, analytes in user-defined types, i.e., `getproperty(table(dt), analytecol(dt))`. Use `analyteobj` to get this field.
* `sample`: `U <: AbstractVector{S}`, samples in user-defined types. Use `sampleobj` to get this field.
* `analytecol`: `Symbol`, the column name that each element is analyte name. Use `analytecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct RowDataTable{A, S, T, V <: AbstractVector{A}, U <: AbstractVector{S}} <: AbstractDataTable{A, S, T, V, U}
    analyte::V
    sample::U
    analytecol::Symbol
    table::T
    function RowDataTable(analyte::V, sample::U, analytecol::Symbol, table::T) where {T, V <: AbstractVector, U <: AbstractVector}
        analyte === getproperty(table, analytecol) || throw(ArgumentError("`analyte` is not identical to `getproperty(table, analytecol)`."))
        new{eltype(analyte), eltype(sample), T, V, U}(analyte, sample, analytecol, table)
    end
end

"""
    RowDataTable(analytecol::Symbol, samplename::AbstractVector{Symbol}, table; sampletype = String)
    RowDataTable(analytecol::Symbol, samplename::AbstractVector{B}, table; sampletype = String)
    RowDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol]), sampletype = String)
    RowDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol]), sampletype = String)

User-friendly contructors for `RowDataTable{A, S}`.

* `analytecol`: `Symbol`, the column name that each element is analyte name.
* `samplename`: `AbstractVector{Symbol}`, the column names of `table` that are sample names. It will be converted to `Vector{String}` before conversion to `AbstractVector{S}`. Vector of other type is also valid and will not be converted to `Vector{String}` first.
* `sampletype`: type `S` or a function. After first conversion of `sampletype`, it will be converted to `AbstractVector{S}` using `cqamap(sampletype, sampletype)`. See `cqaconvert` for the requirement of `sampletype`.
"""
RowDataTable(analytecol::Symbol, samplename::AbstractVector{Symbol}, table; sampletype = String) = 
    RowDataTable(getproperty(table, analytecol), cqamap(sampletype, string.(samplename)), analytecol, table)
RowDataTable(analytecol::Symbol, samplename::AbstractVector, table; sampletype = String) = 
    RowDataTable(getproperty(table, analytecol), cqamap(sampletype, samplename), analytecol, table)
RowDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol]), sampletype = String) = 
    RowDataTable(analytecol, samplename, table; sampletype)
RowDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol]), sampletype = String) = 
    RowDataTable(analytecol, samplename, table; sampletype)

"""
    AnalysisTable{A, S, T}

Tabular data wrapper of multiple tables for analysis. `A` determines analyte type, `S` determines sample type, and `T` determines table type.

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types. Use `analyteobj` to get this field.
* `sample`: `Vector{S}`, samples in user-defined types. Use `sampleobj` to get this field.
* `tables`: `Dictionary{Symbol, <: AbstractDataTable{A, S, <: T}}`, a dictionary mapping data type to datatable. Use `tables` to get this field.

# Properties
All keys of `tables`.

# Constructors
* `AnalysisTable(analyte, sample, tables)`
* `AnalysisTable{T}(analyte, sample, tables)`
"""
struct AnalysisTable{A, S, T}
    analyte::Vector{A}
    sample::Vector{S}
    tables::Dictionary{Symbol, <: AbstractDataTable{A}}
    AnalysisTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, tables::Dictionary{Symbol, <: AbstractDataTable{A, S, T}}) where {A, S, T} = new{A, S, T}(convert(Vector, analyte), convert(Vector, sample), tables)
    AnalysisTable{T}(analyte::AbstractVector{A}, sample::AbstractVector{S}, tables::Dictionary{Symbol, <: AbstractDataTable{A, S, <: T}}) where {A, S, T} = new{A, S, T}(convert(Vector, analyte), convert(Vector, sample), tables)
end

"""
    AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, T}})
    AnalysisTable{T}(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, <: T}})

A `Dictionary`-like constructor for `AnalysisTable`.
"""
AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S, T}}) where {A, S, T} = AnalysisTable{T}(keys, tables)
function AnalysisTable{T}(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S, <: T}}) where {A, S, T}
    allequal(analyteobj(table) for table in tables) || throw(ArgumentError("Tables should have identical analyte"))g
    allequal(sampleobj(table) for table in tables) || throw(ArgumentError("Tables should have identical sample"))
    AnalysisTable{T}(deepcopy(analyteobj(first(tables))), deepcopy(sampleobj(first(tables))), Dictionary(keys, tables))
end

"""
    MethodTable{A, T, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}

Tabular data wrapper for calibration method. `A` determines analyte type, and `T` determines table type.

# Fields
* `analytetable`: `Table` contaning three columns.
    * `analyte`: `AbstractVector{A}`, analytes in user-defined types.
    * `isd`: `AbstractVector{Int}`, index of internal standard. `0` means no internal standard, and `-1` means the analyte itself is a internal standard.
    * `calibration`: index of analyte for calibration curve. `-1` means the analyte itself is a internal standard, so it will not be put into any calibration curve.
* `signal`: `Symbol`, data type for quantification, e.g. `:area`.
* `pointlevel`: `Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.
* `conctable`: `C <: AbstractDataTable{A, Int, <: T}` containing concentration data for each level. Sample names must be symbol or string of integers for multiple levels. One level indicates using `SingleCalibration`.
* `signaltable`: `D <: AbstractDataTable{A, S, <: T}` containig signal for each point. It can be `nothing` if signal data is unecessary.

# Properties
* `analyte`: `AbstractVector{A}`, analytes in user-defined types, identical to `analytetable.analyte`.
* `isd`: `AbstractVector{A}`, analytes which are internal standards.
* `nonisd`: `AbstractVector{A}`, analytes which are not internal standards.
* `point`: `AbstractVector{S}`, calibration points, identical to `sampleobj(signaltable)`. If `signaltable` is `nothing`, this value is `nothing` too.
* `level`: `AbstractVector{Int}`, calibration levels, identical to `sampleobj(conctable)`.

# Constructors
* `MethodTable(analytetable, signal, pointlevel, conctable, signaltable = nothing)`
* `MethodTable{T}(analytetable, signal, pointlevel, conctable, signaltable = nothing)`
"""
struct MethodTable{A, T, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}
    analytetable::Table
    signal::Symbol
    pointlevel::Vector{Int}
    conctable::C
    signaltable::D
    function MethodTable{T}(
                        analytetable::Table,
                        signal::Symbol,
                        pointlevel::AbstractVector{Int}, 
                        conctable::AbstractDataTable{B, Int, <: T}, 
                        signaltable::AbstractDataTable{C, S, <: T}) where {B, C, S, T}
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
        new{A, T, typeof(conctable), typeof(signaltable)}(analytetable, signal, convert(Vector, pointlevel), conctable, signaltable)
    end
    function MethodTable{T}(
                        analytetable::Table,
                        signal::Symbol,
                        pointlevel::AbstractVector{Int}, 
                        conctable::AbstractDataTable{B, Int, <: T}, 
                        signaltable::Nothing) where {B, T}
        cp = propertynames(analytetable)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analytetable"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analytetable"))
        :calibration in cp || throw(ArgumentError("Column `:calibration` is required in `analytetable"))
        A = eltype(analytetable.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable should be a subtype of $A"))
        for a in analyteobj(conctable)
            a in analytetable.analyte || throw(ArgumentError("Analyte `$a` is not in the `analytetable`."))
        end
        new{A, T, typeof(conctable), Nothing}(analytetable, signal, convert(Vector, pointlevel), conctable, signaltable)
    end
end

MethodTable(
            analytetable::Table,
            signal::Symbol,
            pointlevel::AbstractVector{Int}, 
            conctable::AbstractDataTable{B, Int, T}, 
            signaltable::Nothing) where {B, T} = MethodTable{T}(analytetable, signal, pointlevel, conctable, signaltable)

MethodTable(
            analytetable::Table,
            signal::Symbol,
            pointlevel::AbstractVector{Int}, 
            conctable::AbstractDataTable{B, Int, T}, 
            signaltable::AbstractDataTable{C, S, D}) where {B, C, T, S, D} = MethodTable{promote_type(T, D)}(analytetable, signal, pointlevel, conctable, signaltable)

"""
    MethodTable(conctable::AbstractDataTable, signaltable::Union{AbstractDataTable, Nothing}, signal, pointlevel = []; kwargs...)
    MethodTable{T}(conctable::AbstractDataTable, signaltable::Union{AbstractDataTable, Nothing}, signal, pointlevel = []; kwargs...)
    MethodTable(conctable::AbstractDataTable, signaltable::ColumnDataTable, signal::Symbol, levelname::Symbol; kwargs...)

User-friendly contructors for `MethodTable`. `kwargs` will be columns in `analytetable`; when `analyte`, `isd` and `calibration` are not provided, it will use analyte in `conctable`.
"""
MethodTable(conctable::AbstractDataTable{A, Int, T}, 
            signaltable::AbstractDataTable{B, S, D},
            signal::Symbol,
            pointlevel::AbstractVector{Int}; kwargs...) where {A, B, T, S, D} = MethodTable{promote_type(T, D)}(conctable, signaltable, signal, pointlevel; kwargs...)
MethodTable(conctable::AbstractDataTable{A, Int, T}, 
            signaltable::Nothing,
            signal::Symbol,
            pointlevel::AbstractVector{Int} = Int[]; kwargs...) where {A, T} = MethodTable{A, T}(conctable, signaltable, signal, pointlevel; kwargs...)
MethodTable(conctable::AbstractDataTable{A, Int, T}, 
            signaltable::ColumnDataTable{B, S, D},
            signal::Symbol,
            levelname::Symbol; kwargs...) where {A, B, T, S, D} = MethodTable{promote_type(T, D)}(conctable, signaltable, signal, getproperty(signaltable, levelname); kwargs...)
function MethodTable{T}(conctable::AbstractDataTable, 
                    signaltable::Union{AbstractDataTable, Nothing},
                    signal::Symbol,
                    pointlevel::AbstractVector{Int}; kwargs...
            ) where T
    analytetable = if isempty(kwargs)
        Table(; analyte = analyteobj(conctable))
    else
        Table(; kwargs...)
    end
    if !in(:analyte, propertynames(analytetable))
        analytetable = Table(analytetable; analyte = analyteobj(conctable))
    end       
    if !in(:isd, propertynames(analytetable))
        analytetable = Table(analytetable; isd = zeros(Int, length(analyteobj(conctable))))
    end
    if !in(:calibration, propertynames(analytetable))
        analytetable = Table(analytetable; calibration = collect(eachindex(analyteobj(conctable))))
    end
    analytetable = Table((; analyte = analytetable.analyte, isd = analytetable.isd, calibration = analytetable.calibration), analytetable)    
    MethodTable{T}(analytetable, signal, pointlevel, conctable, signaltable)
end
      
"""
    MultipleCalibration{A} <: AbstractCalibration{A}

A type holding all data for a calibration curve.

# Fields
* `analyte`: `Tuple{A, Any}`. First element is the analyte being quantified, and the second element is its internal standard for which `nothing` indicates no internal standard.
* `type`: `Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`).
* `zero`: `Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`).
* `weight`: `Float64` represents the exponential applying to each element of `x` as a weighting vector.
* `formula`: `FormulaTerm`, the formula for fitting calibration curve.
* `table`: `TypedTable.Table`, the clean up calibration data, containing 7 columns.
    * `id`: Point name
    * `level`: The index of concentration level. The larger, the higher concentraion it represents.
    * `y`: Signal or relative signal
    * `x`: True concentraion
    * `x̂`: Predicted concentration
    * `accuracy`: Accuracy, i.e. `x̂/x`.
    * `include`: Whether this point is included or not
* `model`: `GLM` object.
"""
mutable struct MultipleCalibration{A} <: AbstractCalibration{A}
    analyte::Tuple{A, Any}
    type::Bool
    zero::Bool
    weight::Float64
    formula::FormulaTerm
    table::Table
    model
end

"""
    SingleCalibration{A} <: AbstractCalibration{A}

A type holding all data for single point calibration.

# Fields
* `analyte`: `Tuple{A}` is the analyte with known concentration (internal standard).
* `conc`: `Float64`, concentration of analyte.
"""
mutable struct SingleCalibration{A} <: AbstractCalibration{A}
    analyte::Tuple{A}
    conc::Float64
end

"""
    Batch{A, T}

A type represents a batch for quantitative analysis. `A` determines analyte type, and `T` determines table type.

# Fields
* `method`: `MethodTable{A, <: T}`, method.
* `calibration`: `AbstractVector{MultipleCalibration{<: A}}` or `AbstractVector{SingleCalibration{<: A}}`.
* `data`: Data for analysis, `AnalysisTable{A, S, <: T}` or `Nothing`.

# Properties
* `analyte`: `AbstractVector{A}`, analytes in user-defined types, identical to `method.analytetable.analyte`.
* `isd`: `AbstractVector{<: A}`, analytes which are internal standards.
* `nonisd`: `AbstractVector{<: A}`, analytes which are not internal standards.
* `point`: `AbstractVector{X}` or `Nothing`, calibration points, identical to `method.point`.
* `level`: `AbstractVector{Int}`, calibration levels, identical to `method.level`.

# Constructors
* `Batch(method, calibration, data = nothing)`
* `Batch{T}(method, calibration, data = nothing)`
"""
mutable struct Batch{A, T}
    method::MethodTable{A}
    calibration::AbstractVector{<: AbstractCalibration}
    data::Union{AnalysisTable, Nothing}
    Batch{T}(method::MethodTable{A, <: T}, calibration::AbstractVector{<: AbstractCalibration{<: A}}, data::AnalysisTable{<: A, <: S, <: T}) where {A, S, T} = new{A, T}(method, calibration, data)
    Batch{T}(method::MethodTable{A, <: T}, calibration::AbstractVector{<: AbstractCalibration{<: A}}, data::Nothing) where {A, T} = new{A, T}(method, calibration, data)
    Batch{T}(method::MethodTable{A, <: T}, calibration::AbstractVector{<: AbstractCalibration{<: A}}) where {A, T} = new{A, T}(method, calibration, nothing)
    Batch(method::MethodTable{A, T}, calibration::AbstractVector{<: AbstractCalibration{<: A}}, data::AnalysisTable{<: A, S, D}) where {A, T, S, D} = new{A, promote_type(T, D)}(method, calibration, data)
    Batch(method::MethodTable{A, T}, calibration::AbstractVector{<: AbstractCalibration{<: A}}, data::Nothing) where {A, T} = new{A, T}(method, calibration, data)
    Batch(method::MethodTable{A, T}, calibration::AbstractVector{<: AbstractCalibration{<: A}}) where {A, T} = new{A, T}(method, calibration, nothing)
end

"""
    Batch(method, data = nothing; type = true, zero = false, weight = 0)
    Batch{T}(method, data = nothing; type = true, zero = false, weight = 0)

A user-friendly constructor for `Batch` using `method`, and optionally `data` with specified calibration parameters. See "MultipleCalibration" for detail description of keyword arguments.
"""
Batch(method::MethodTable{A, T}, data = nothing;
                type = true, 
                zero = false, 
                weight = 0
                ) where {A, T} = Batch{T}(method, data; type, zero, weight)
Batch(method::MethodTable{A, T}, data::AnalysisTable{<: A, S, D};
                type = true, 
                zero = false, 
                weight = 0
                ) where {A, T, S, D} = Batch{promote_type(T, D)}(method, data; type, zero, weight)
function Batch{T}(method::MethodTable{A, <: T}, data = nothing;
                type = true, 
                zero = false, 
                weight = 0
                ) where {A, T}
    Batch{T}(
        method,
        if length(sampleobj(method.conctable)) > 1 
            analyteobj(method.conctable) == analyteobj(method.signaltable) ? map(eachindex(analyteobj(method.conctable))) do i
                calibration(method, i; type, zero, weight)
            end : map(analyteobj(method.conctable)) do analyte
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

include("interface.jl")
include("utils.jl")
include("cal.jl")
include("io.jl")

"""
    ColumnDataTable(tbl::RowDataTable{A, S, T}, samplecol::Symbol, tablesink = T)

Convert `RowDataTable` to `ColumnDataTable` with `samplecol` as the column name of sample.
"""
ColumnDataTable(tbl::RowDataTable{A, S, T}, samplecol::Symbol, tablesink = T) where {A, S, T} = 
    ColumnDataTable(analyteobj(tbl), samplecol, tablesink((; (samplecol => sampleobj(tbl), (analytename(tbl) .=> eachanalyte(tbl))...)...)); analytetype = A)

"""
    RowDataTable(tbl::ColumnDataTable{A, T}, analytecol::Symbol, tablesink = T)

Convert `RowDataTable` to `ColumnDataTable` with `analytecol` as the column name of sample.
"""
RowDataTable(tbl::ColumnDataTable{A, S, T}, analytecol::Symbol, tablesink = T) where {A, S, T} = 
    RowDataTable(analytecol, sampleobj(tbl), tablesink((; (analytecol => analyteobj(tbl), (samplename(tbl) .=> eachsample(tbl))...)...)); sampletype = S)

end
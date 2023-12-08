module ChemistryQuantitativeAnalysis

using GLM, CSV, TypedTables, LinearAlgebra, Dictionaries, ThreadsX
export MultipleCalibration, SingleCalibration, 
    ColumnDataTable, RowDataTable, AnalysisTable, MethodTable,
    Batch, calibration, update_calibration!,
    read_calibration, read_datatable, read_analysistable, read_methodtable, read_batch,
    cal_range, lloq, uloq, accuracy, accuracy!, set_accuracy, set_accuracy!, update_accuracy!,
    inv_predict, inv_predict!, inv_predict_accuracy!, set_inv_predict, set_inv_predict!, update_inv_predict!,
    relative_signal, set_relative_signal, set_relative_signal!, update_relative_signal!,
    quantification, quantification!, set_quantification, set_quantification!, update_quantification!,
    find_analyte, get_analyte, find_sample, get_sample, set_isd!,
    formula_repr, weight_repr, weight_value, formula_repr_utf8, weight_repr_utf8, format_number

import Base: getproperty, show, write
    
abstract type AbstractCalibration{A} end
abstract type AbstractDataTable{A, T} end
abstract type AbstractAnalysisTable{A, T} end

"""
    ColumnDataTable{A, T} <: AbstractDataTable{A, T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples. `A` determines analyte type, and `T` determines table type.

# Fields
* `sample_name`: `Symbol`, the column name that each element is sample name.
* `config`: `Table` with 2 columns, `analyte_name`, `analytes` which are properties of this type. To add new analytes, user should modify this table.
* `table`: Tabular data of type `T`.

# Properties
* `analyte_name`: `Vector{Symbol}`, the column names of `table` that are analytes names.
* `analytes`: `Vector{A}`, analytes in user-defined types.
* `samples`: `Vector{Symbol}`, sample names.
"""
struct ColumnDataTable{A, T} <: AbstractDataTable{A, T}
    sample_name::Symbol
    config::Table{NamedTuple{(:analyte_name, :analytes), Tuple{Symbol, A}}, 1, NamedTuple{(:analyte_name, :analytes), Tuple{Vector{Symbol}, Vector{A}}}}
    table::T
end

"""
    ColumnDataTable(sample_name::Symbol, analyte_name::Vector{Symbol}, analytes::Vector, table)
    ColumnDataTable(sample_name::Symbol, analyte_name::Vector, table; analyte_fn = default_analyte_fn)
    ColumnDataTable(sample_name::Symbol, table; analyte_name = setdiff(propertynames(table), [sample_name]), analyte_fn = default_analyte_fn)
    ColumnDataTable(table, sample_name::Symbol; analyte_name = setdiff(propertynames(table), [sample_name]), analyte_fn = default_analyte_fn)

User-friendly contructors for `ColumnDataTable`. See the documentation of the type for detail description.
"""
function ColumnDataTable(sample_name::Symbol, analyte_name::Vector{Symbol}, analytes::Vector, table)
    length(analyte_name) == length(analytes) || throw(ArgumentError("Arguments `analyte_name` and `analytes` should have the same length."))
    ColumnDataTable(sample_name, Table(; analyte_name, analytes), table)
end

ColumnDataTable(sample_name::Symbol, analyte_name::Vector, table; analyte_fn = default_analyte_fn) = 
    ColumnDataTable(sample_name, Symbol.(analyte_name), analyte_fn.(string.(analyte_name)), table)
ColumnDataTable(sample_name::Symbol, table; analyte_name = setdiff(propertynames(table), [sample_name]), analyte_fn = default_analyte_fn) = 
    ColumnDataTable(sample_name, analyte_name, table; analyte_fn)
ColumnDataTable(table, sample_name::Symbol; analyte_name = setdiff(propertynames(table), [sample_name]), analyte_fn = default_analyte_fn) = 
    ColumnDataTable(sample_name, analyte_name, table; analyte_fn) 

function getproperty(tbl::ColumnDataTable, property::Symbol)
    if property in [:analyte_name, :analytes] 
        getproperty(getfield(tbl, :config), property)
    elseif property == :samples
        Symbol.(getproperty(tbl.table, tbl.sample_name))
    elseif !in(property, fieldnames(ColumnDataTable))
        getproperty(getfield(tbl, :table), property)
    else
        getfield(tbl, property)
    end
end

"""
    RowDataTable{A, T} <: AbstractDataTable{A, T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples. `A` determines analyte type, and `T` determines table type.

# Fields
* `sample_name`: `Vector{Symbol}`, the column names that are sample names.
* `analyte_name`: `Symbol`, the column name that each element is analyte name.
* `analytes`: `Vector{A}`, analytes in user-defined types.
* `table`: Tabular data of type `T`.

# Properties
* `samples`: `Vector{Symbol}`, sample names.
"""
struct RowDataTable{A, T} <: AbstractDataTable{A, T}
    sample_name::Vector{Symbol}
    analyte_name::Symbol
    analytes::Vector{A}
    table::T
    function RowDataTable(sample_name::Vector{Symbol}, analyte_name::Symbol, analytes::Vector{A}, table::T) where {A, T}
        length(analytes) == length(getproperty(table, analyte_name)) || throw(ArgumentError("The length of `analytes` is different from that of `table`."))
        new{A, T}(sample_name, analyte_name, analytes, table)
    end
end

"""
    RowDataTable(sample_name::Vector{Symbol}, analyte_name::Symbol, table; analyte_fn = default_analyte_fn)
    RowDataTable(analyte_name::Symbol, table; sample_name = setdiff(propertynames(table), [analyte_name]), analyte_fn = default_analyte_fn)
    RowDataTable(table, analyte_name::Symbol; sample_name = setdiff(propertynames(table), [analyte_name]), analyte_fn = default_analyte_fn)

User-friendly contructors for `RowDataTable`. See the documentation of the type for detail description.
"""
RowDataTable(sample_name::Vector{Symbol}, analyte_name::Symbol, table; analyte_fn = default_analyte_fn) = 
    RowDataTable(sample_name, analyte_name, analyte_fn.(string.(getproperty(table, analyte_name))), table)
RowDataTable(analyte_name::Symbol, table; sample_name = setdiff(propertynames(table), [analyte_name]), analyte_fn = default_analyte_fn) = 
    RowDataTable(sample_name, analyte_name, table; analyte_fn)
RowDataTable(table, analyte_name::Symbol; sample_name = setdiff(propertynames(table), [analyte_name]), analyte_fn = default_analyte_fn) = 
    RowDataTable(sample_name, analyte_name, table; analyte_fn)

function getproperty(tbl::RowDataTable{A}, property::Symbol) where A
    if property == :samples
        getfield(tbl, :sample_name)
    elseif !in(property, fieldnames(RowDataTable))
        getproperty(getfield(tbl, :table), property)
    else
        getfield(tbl, property)
    end
end

"""
    AnalysisTable{A, T} <: AbstractAnalysisTable{A, T}

Tabular data wrapper of multiple tables for analysis. `A` determines analyte type, and `T` determines table type.

# Fields
* `analytes`: `Vector{A}`, analytes in user-defined types.
* `samples`: `Vector{Symbol}`, sample names.
* `tables`: `Dictionary{Symbol, <: AbstractDataTable{A, <: T}}`, a dictionary mapping data type to datatable.

# Properties
All keys of `tables` are available.

# Constructors
* `AnalysisTable(analytes, samples, tables)`
* `AnalysisTable{T}(analytes, samples, tables)`
"""
struct AnalysisTable{A, T} <: AbstractAnalysisTable{A, T}
    analytes::Vector{A}
    samples::Vector{Symbol}
    tables::Dictionary{Symbol, <: AbstractDataTable{A}}
    AnalysisTable(analytes::Vector{A}, samples::Vector{Symbol}, tables::Dictionary{Symbol, <: AbstractDataTable{A, T}}) where {A, T} = new{A, T}(analytes, samples, tables)
    AnalysisTable{T}(analytes::Vector{A}, samples::Vector{Symbol}, tables::Dictionary{Symbol, <: AbstractDataTable{A, <: T}}) where {A, T} = new{A, T}(analytes, samples, tables)
end

"""
    AnalysisTable(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, T}})
    AnalysisTable{T}(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, <: T}})

A `Dictionary`-like constructor for `AnalysisTable`.
"""
AnalysisTable(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, T}}) where {A, T} = AnalysisTable{T}(keys, tables)
function AnalysisTable{T}(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, <: T}}) where {A, T}
    allequal(table.analytes for table in tables) || throw(ArgumentError("Tables should have identical analytes"))
    allequal(table.samples for table in tables) || throw(ArgumentError("Tables should have identical samples"))
    AnalysisTable{T}(deepcopy(first(tables).analytes), deepcopy(first(tables).samples), Dictionary(keys, tables))
end

function getproperty(tbl::AnalysisTable, p::Symbol)
    if p in fieldnames(AnalysisTable)
        getfield(tbl, p)
    else
        get(getfield(tbl, :tables), p, nothing)
    end
end

"""
    MethodTable{A, T}

Tabular data wrapper for calibration method. `A` determines analyte type, and `T` determines table type.

# Fields
* `signal`: `Symbol`, data type for quantification, e.g. `:area`.
* `analyte_map`: `Table` contaning three columns.
    * `analytes`: `Vector{A}`, analytes in user-defined types.
    * `isd`: `Vector{Int}`, index of internal standard. `0` means no internal standard, and `-1` means the analyte itself is a internal standard.
    * `calibration`: index of analyte for calibration curve. `-1` means the analyte itself is a internal standard, so it will not be put into any calibration curve.
* `level_map`: `Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.
* `conctable`: `AbstractDataTable{A, <: T}` containing concentration data for each level. Sample names must be symbol or string of integers for multiple levels. One level indicates using `SingleCalibration`.
* `signaltable`: `AbstractDataTable{A, <: T}` containig signal for each point. It can be `nothing` if signal data is unecessary.

# Properties
* `analytes`: `Vector{A}`, analytes in user-defined types, identical to `analyte_map.analytes`.
* `isds`: `Vector{A}`, analytes which are internal standards.
* `nonisds`: `Vector{A}`, analytes which are not internal standards.
* `points`: `Vector{Symbol}`, calibration points, identical to `signaltable.samples`.
* `levels`: `Vector{Symbol}`, calibration levels, identical to `conctable.samples`.

# Constructors
* `MethodTable(signal, analyte_map, level_map, conctable, signaltable = nothing)`
* `MethodTable{T}(signal, analyte_map, level_map, conctable, signaltable = nothing)`
"""
struct MethodTable{A, T} <: AbstractAnalysisTable{A, T}
    signal::Symbol
    analyte_map::Table{NamedTuple{(:analytes, :isd, :calibration), Tuple{A, Int, Int}}, 1, NamedTuple{(:analytes, :isd, :calibration), Tuple{Vector{A}, Vector{Int}, Vector{Int}}}}
    level_map::Vector{Int}
    conctable::AbstractDataTable
    signaltable::Union{AbstractDataTable, Nothing}
    function MethodTable{T}(
                        signal::Symbol,
                        analyte_map::Table{NamedTuple{(:analytes, :isd, :calibration), Tuple{A, Int, Int}}, 1, NamedTuple{(:analytes, :isd, :calibration), Tuple{Vector{A}, Vector{Int}, Vector{Int}}}},
                        level_map::Vector{Int}, 
                        conctable::AbstractDataTable{<: A, <: T}, 
                        signaltable::AbstractDataTable{<: A, <: T}) where {A, T}
        for a in conctable.analytes
            a in analyte_map.analytes || throw(ArgumentError("Analyte `$a` is not in the `analyte_map`."))
        end
        if length(conctable.samples) > 1
            length(level_map) == length(signaltable.samples) || throw(ArgumentError("The length of `level_map` is different from that of `table.samples`."))
            for a in conctable.analytes
                a in signaltable.analytes || throw(ArgumentError("Analyte `$a` is not in the `signatable`."))
            end
        end
        new{A, T}(signal, analyte_map, level_map, conctable, signaltable)
    end
    function MethodTable{T}(
                        signal::Symbol,
                        analyte_map::Table{NamedTuple{(:analytes, :isd, :calibration), Tuple{A, Int, Int}}, 1, NamedTuple{(:analytes, :isd, :calibration), Tuple{Vector{A}, Vector{Int}, Vector{Int}}}},
                        level_map::Vector{Int}, 
                        conctable::AbstractDataTable{<: A, <: T}, 
                        signaltable::Nothing) where {A, T}
        for a in conctable.analytes
            a in analyte_map.analytes || throw(ArgumentError("Analyte `$a` is not in the `analyte_map`."))
        end
        new{A, T}(signal, analyte_map, level_map, conctable, signaltable)
    end
end

MethodTable(
            signal::Symbol,
            analyte_map::Table{NamedTuple{(:analytes, :isd, :calibration), Tuple{A, Int, Int}}, 1, NamedTuple{(:analytes, :isd, :calibration), Tuple{Vector{A}, Vector{Int}, Vector{Int}}}},
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{<: A, T}, 
            signaltable::Nothing) where {A, T} = MethodTable{T}(signal, analyte_map, level_map, conctable, signaltable)

MethodTable(
            signal::Symbol,
            analyte_map::Table{NamedTuple{(:analytes, :isd, :calibration), Tuple{A, Int, Int}}, 1, NamedTuple{(:analytes, :isd, :calibration), Tuple{Vector{A}, Vector{Int}, Vector{Int}}}},
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{<: A, T}, 
            signaltable::AbstractDataTable{<: A, S}) where {A, T, S} = MethodTable{promote_type(T, S)}(signal, analyte_map, level_map, conctable, signaltable)

function getproperty(tbl::MethodTable, p::Symbol)
    if p == :analytes
        getfield(tbl, :analyte_map).analytes
    elseif p == :isds
        getfield(tbl, :analyte_map).analytes[getfield(tbl, :analyte_map).isd .< 0]
    elseif p == :nonisds
        getfield(tbl, :analyte_map).analytes[getfield(tbl, :analyte_map).isd .>= 0]
    elseif p == :points
        s = getfield(tbl, :signaltable)
        isnothing(s) ? s : s.samples
    elseif p == :levels
        getfield(tbl, :conctable).samples
    else
        getfield(tbl, p)
    end
end

"""
    MethodTable(signal, analytes, isd, calibration, level_map, conctable, signaltable)
    MethodTable(signal, level_map, conctable, signaltable)
    MethodTable)T}(signal, analytes, isd, calibration, level_map, conctable, signaltable)
    MethodTable{T}(signal, level_map, conctable, signaltable)

User-friendly contructors for `MethodTable`. When `analytes`, `isd` and `calibration` are not provided, it will use analytes in `conctable`.
"""
MethodTable(signal::Symbol,
            analytes::Vector{A}, isd::Vector{Int}, calibration::Vector{Int},
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{<: A, T}, 
            signaltable::Union{AbstractDataTable{<: A, S}, Nothing}) where {A, T, S} = MethodTable(signal, Table(; analytes, isd, calibration), level_map, conctable, signaltable)
MethodTable(signal::Symbol,
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{<: A, T}, 
            signaltable::Union{AbstractDataTable{<: A, S}, Nothing}) where {A, T, S} = MethodTable(signal, Table(; analytes = conctable.analytes, isd = zeros(Int, length(conctable.analytes)), calibration = collect(eachindex(conctable.analytes))), level_map, conctable, signaltable) 
MethodTable{T}(signal::Symbol,
            analytes::Vector{A}, isd::Vector{Int}, calibration::Vector{Int},
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{<: A, <: T},
            signaltable::Union{AbstractDataTable{<: A, <: T}, Nothing}) where {A, T} = MethodTable{T}(signal, Table(; analytes, isd, calibration), level_map, conctable, signaltable)
MethodTable{T}(signal::Symbol,
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{<: A, <: T}, 
            signaltable::Union{AbstractDataTable{<: A, <: T}, Nothing}) where {A, T} = MethodTable{T}(signal, Table(; analytes = conctable.analytes, isd = zeros(Int, length(conctable.analytes)), calibration = collect(eachindex(conctable.analytes))), level_map, conctable, signaltable)       
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
* `calibration`: `Vector{MultipleCalibration{<: A}}` or `Vector{SingleCalibration{<: A}}`.
* `data`: Data for analysis, `AnalysisTable{A, <: T}` or `Nothing`.

# Properties
* `analytes`: `Vector{A}`, analytes in user-defined types, identical to `method.analyte_map.analytes`.
* `isds`: `Vector{<: A}`, analytes which are internal standards.
* `nonisds`: `Vector{<: A}`, analytes which are not internal standards.
* `points`: `Vector{Symbol}`, calibration points, identical to `method.signaltable.samples`.
* `levels`: `Vector{Symbol}`, calibration levels, identical to `method.conctable.samples`.

# Constructors
* `Batch(method, calibration, data = nothing)`
* `Batch{T}(method, calibration, data = nothing)`
"""
mutable struct Batch{A, T}
    method::MethodTable{A}
    calibration::Vector{<: AbstractCalibration}
    data::Union{AnalysisTable, Nothing}
    Batch{T}(method::MethodTable{A, <: T}, calibration::Vector{<: AbstractCalibration{<: A}}, data::AnalysisTable{<: A, <: T}) where {A, T} = new{A, T}(method, calibration, data)
    Batch{T}(method::MethodTable{A, <: T}, calibration::Vector{<: AbstractCalibration{<: A}}, data::Nothing) where {A, T} = new{A, T}(method, calibration, data)
    Batch{T}(method::MethodTable{A, <: T}, calibration::Vector{<: AbstractCalibration{<: A}}) where {A, T} = new{A, T}(method, calibration, nothing)
    Batch(method::MethodTable{A, T}, calibration::Vector{<: AbstractCalibration{<: A}}, data::AnalysisTable{<: A, S}) where {A, T, S} = new{A, promote_type(T, S)}(method, calibration, data)
    Batch(method::MethodTable{A, T}, calibration::Vector{<: AbstractCalibration{<: A}}, data::Nothing) where {A, T} = new{A, T}(method, calibration, data)
    Batch(method::MethodTable{A, T}, calibration::Vector{<: AbstractCalibration{<: A}}) where {A, T} = new{A, T}(method, calibration, nothing)
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
Batch(method::MethodTable{A, T}, data::AnalysisTable{<: A, S};
                type = true, 
                zero = false, 
                weight = 0
                ) where {A, T, S} = Batch{promote_type(T, S)}(method, data; type, zero, weight)
function Batch{T}(method::MethodTable{A, <: T}, data = nothing;
                type = true, 
                zero = false, 
                weight = 0
                ) where {A, T}
    Batch{T}(
        method,
        length(method.conctable.samples) > 1 ? map(method.conctable.analytes) do analyte
            calibration(method, analyte; type, zero, weight)
        end : map(method.conctable.analytes) do analyte
            SingleCalibration((analyte, ), first(get_analyte(method.conctable, analyte)))
        end
        ,
        data
    )
end

function getproperty(batch::Batch, p::Symbol)
    if p == :analytes
        getfield(batch, :method).analytes
    elseif p == :isds
        getfield(batch, :method).analytes[getfield(batch, :method).isd .< 0]
    elseif p == :nonisds
        getfield(batch, :method).analytes[getfield(batch, :method).isd .>= 0]
    elseif p == :points
        getfield(batch, :method).points
    elseif p == :levels
        getfield(batch, :method).levels
    else
        getfield(batch, p)
    end
end

include("utils.jl")
include("cal.jl")
include("io.jl")

end

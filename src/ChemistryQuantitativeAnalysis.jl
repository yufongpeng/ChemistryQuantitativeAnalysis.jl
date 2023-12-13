module ChemistryQuantitativeAnalysis

using GLM, CSV, TypedTables, LinearAlgebra, Dictionaries, ThreadsX, Tables
import Tables: istable, rowaccess, rows, columnaccess, columns
export MultipleCalibration, SingleCalibration, 
    ColumnDataTable, RowDataTable, AnalysisTable, MethodTable,
    Batch, calibration, update_calibration!,
    read_calibration, read_datatable, read_analysistable, read_methodtable, read_batch,
    cal_range, lloq, uloq, accuracy, accuracy!, set_accuracy, set_accuracy!, update_accuracy!,
    inv_predict, inv_predict!, inv_predict_accuracy!, set_inv_predict, set_inv_predict!, update_inv_predict!,
    relative_signal, set_relative_signal, set_relative_signal!, update_relative_signal!,
    quantification, quantification!, set_quantification, set_quantification!, update_quantification!,
    findanalyte, getanalyte, findsample, getsample, set_isd!,
    formula_repr, weight_repr, weight_value, formula_repr_utf8, weight_repr_utf8, format_number

import Base: getproperty, show, write, eltype, length, iterate
    
abstract type AbstractCalibration{A} end
abstract type AbstractDataTable{A, T} end
abstract type AbstractAnalysisTable{A, T} end

"""
    ColumnDataTable{A, T} <: AbstractDataTable{A, T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples. `A` determines analyte type, and `T` determines table type.

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types.
* `samplecol`: `Symbol`, the column name that each element is sample name.
* `table`: Tabular data of type `T`.

# Properties
* `analytename`: `Vector{Symbol}`, the column names of `table` that are analyte names.
* `samplename`: `Vector{Symbol}`, sample names in type symbol.
* `sample`: `Vector{String}`, sample names in type string.
"""
struct ColumnDataTable{A, T} <: AbstractDataTable{A, T}
    analyte::Vector{A}
    samplecol::Symbol
    table::T
end

"""
    ColumnDataTable(analyte_name::Vector{Symbol}, samplecol::Symbol, table; analyte_fn = default_analyte_fn)
    ColumnDataTable(samplecol::Symbol, table; analyte_name = setdiff(propertynames(table), [samplecol]), analyte_fn = default_analyte_fn)
    ColumnDataTable(table, samplecol::Symbol; analyte_name = setdiff(propertynames(table), [samplecol]), analyte_fn = default_analyte_fn)

User-friendly contructors for `ColumnDataTable`. See the documentation of the type for detail description.
"""
ColumnDataTable(analyte_name::Vector{Symbol}, samplecol::Symbol, table; analyte_fn = default_analyte_fn) = 
    ColumnDataTable(analyte_fn.(string.(analyte_name)), samplecol, table)
ColumnDataTable(samplecol::Symbol, table; analyte_name = setdiff(propertynames(table), [samplecol]), analyte_fn = default_analyte_fn) = 
    ColumnDataTable(analyte_name, samplecol, table; analyte_fn)
ColumnDataTable(table, samplecol::Symbol; analyte_name = setdiff(propertynames(table), [samplecol]), analyte_fn = default_analyte_fn) = 
    ColumnDataTable(analyte_name, samplecol, table; analyte_fn)

function getproperty(tbl::ColumnDataTable, property::Symbol)
    if property == :sample
        string.(getproperty(tbl.table, tbl.samplecol))
    elseif property == :samplename
        Symbol.(getproperty(tbl.table, tbl.samplecol))
    elseif property == :analytename
        Symbol.(getfield(tbl, :analyte))
    elseif !in(property, fieldnames(ColumnDataTable))
        getproperty(getfield(tbl, :table), property)
    else
        getfield(tbl, property)
    end
end

"""
    RowDataTable{A, T} <: AbstractDataTable{A, T}

Tabular data wrapper indicates part of columns represent analyte, and all rows reprsent samples. `A` determines analyte type, and `T` determines table type.

# Fields
* `analyte`: `Vector{A}`, analyte in user-defined types.
* `analytecol`: `Symbol`, the column name that each element is analyte name.
* `samplename`: `Vector{Symbol}`, the column names that are sample names.
* `table`: Tabular data of type `T`.

# Properties
* `analytename`: `Symbol`, the column name that each element is analyte name.
* `sample`: `Vector{String}`, sample names in type string.
"""
struct RowDataTable{A, T} <: AbstractDataTable{A, T}
    analyte::Vector{A}
    analytecol::Symbol
    samplename::Vector{Symbol}
    table::T
    function RowDataTable(analyte::Vector{A}, analytecol::Symbol, sample_name::Vector{Symbol}, table::T) where {A, T}
        length(analyte) == length(getproperty(table, analytecol)) || throw(ArgumentError("The length of `analyte` is different from that of `table`."))
        new{A, T}(analyte, analytecol, sample_name, table)
    end
end

"""
    RowDataTable(analytecol::Symbol, sample_name::Vector{Symbol}, table; analyte_fn = default_analyte_fn)
    RowDataTable(analytecol::Symbol, table; sample_name = setdiff(propertynames(table), [analytecol]), analyte_fn = default_analyte_fn)
    RowDataTable(table, analytecol::Symbol; sample_name = setdiff(propertynames(table), [analytecol]), analyte_fn = default_analyte_fn)

User-friendly contructors for `RowDataTable`. See the documentation of the type for detail description.
"""
RowDataTable(analytecol::Symbol, sample_name::Vector{Symbol}, table; analyte_fn = default_analyte_fn) = 
    RowDataTable(analyte_fn.(string.(getproperty(table, analytecol))), analytecol, sample_name, table)
RowDataTable(analytecol::Symbol, table; sample_name = setdiff(propertynames(table), [analytecol]), analyte_fn = default_analyte_fn) = 
    RowDataTable(analytecol, sample_name, table; analyte_fn)
RowDataTable(table, analytecol::Symbol; sample_name = setdiff(propertynames(table), [analytecol]), analyte_fn = default_analyte_fn) = 
    RowDataTable(analytecol, sample_name, table; analyte_fn)

function getproperty(tbl::RowDataTable{A}, property::Symbol) where A
    if property == :sample
        string.(getfield(tbl, :samplename))
    elseif property == :analytename
        Symbol.(getfield(tbl, :analyte))
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
* `analyte`: `Vector{A}`, analyte in user-defined types.
* `sample`: `Vector{String}`, sample names.
* `tables`: `Dictionary{Symbol, <: AbstractDataTable{A, <: T}}`, a dictionary mapping data type to datatable.

# Properties
All keys of `tables` are available.

# Constructors
* `AnalysisTable(analyte, sample, tables)`
* `AnalysisTable{T}(analyte, sample, tables)`
"""
struct AnalysisTable{A, T} <: AbstractAnalysisTable{A, T}
    analyte::Vector{A}
    sample::Vector{String}
    tables::Dictionary{Symbol, <: AbstractDataTable{A}}
    AnalysisTable(analyte::Vector{A}, sample::Vector{String}, tables::Dictionary{Symbol, <: AbstractDataTable{A, T}}) where {A, T} = new{A, T}(analyte, sample, tables)
    AnalysisTable{T}(analyte::Vector{A}, sample::Vector{String}, tables::Dictionary{Symbol, <: AbstractDataTable{A, <: T}}) where {A, T} = new{A, T}(analyte, sample, tables)
end

"""
    AnalysisTable(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, T}})
    AnalysisTable{T}(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, <: T}})

A `Dictionary`-like constructor for `AnalysisTable`.
"""
AnalysisTable(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, T}}) where {A, T} = AnalysisTable{T}(keys, tables)
function AnalysisTable{T}(keys::Vector{Symbol}, tables::Vector{<: AbstractDataTable{A, <: T}}) where {A, T}
    allequal(table.analyte for table in tables) || throw(ArgumentError("Tables should have identical analyte"))
    allequal(table.sample for table in tables) || throw(ArgumentError("Tables should have identical sample"))
    AnalysisTable{T}(deepcopy(first(tables).analyte), deepcopy(first(tables).sample), Dictionary(keys, tables))
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
* `analyte_map`: `Table` contaning three columns.
    * `analyte`: `Vector{A}`, analytes in user-defined types.
    * `isd`: `Vector{Int}`, index of internal standard. `0` means no internal standard, and `-1` means the analyte itself is a internal standard.
    * `calibration`: index of analyte for calibration curve. `-1` means the analyte itself is a internal standard, so it will not be put into any calibration curve.
* `signal`: `Symbol`, data type for quantification, e.g. `:area`.
* `level_map`: `Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.
* `conctable`: `AbstractDataTable{A, <: T}` containing concentration data for each level. Sample names must be symbol or string of integers for multiple levels. One level indicates using `SingleCalibration`.
* `signaltable`: `AbstractDataTable{A, <: T}` containig signal for each point. It can be `nothing` if signal data is unecessary.

# Properties
* `analyte`: `Vector{A}`, analytes in user-defined types, identical to `analyte_map.analyte`.
* `isd`: `Vector{A}`, analytes which are internal standards.
* `nonisd`: `Vector{A}`, analytes which are not internal standards.
* `point`: `Vector{Symbol}`, calibration points, identical to `signaltable.sample`.
* `level`: `Vector{Symbol}`, calibration levels, identical to `conctable.sample`.

# Constructors
* `MethodTable(signal, analyte_map, level_map, conctable, signaltable = nothing)`
* `MethodTable{T}(signal, analyte_map, level_map, conctable, signaltable = nothing)`
"""
struct MethodTable{A, T} <: AbstractAnalysisTable{A, T}
    analyte_map::Table
    signal::Symbol
    level_map::Vector{Int}
    conctable::AbstractDataTable
    signaltable::Union{AbstractDataTable, Nothing}
    function MethodTable{T}(
                        analyte_map::Table,
                        signal::Symbol,
                        level_map::Vector{Int}, 
                        conctable::AbstractDataTable{B, <: T}, 
                        signaltable::AbstractDataTable{C, <: T}) where {B, C, T}
        cp = propertynames(analyte_map)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analyte_map"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analyte_map"))
        :calibration in cp || throw(ArgumentError("Column `:calibration` is required in `analyte_map"))
        A = eltype(analyte_map.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable should be a subtype of $A"))
        C <: A || throw(ArgumentError("Analyte type of signaltable should be a subtype of $A"))
        for a in conctable.analyte
            a in analyte_map.analyte || throw(ArgumentError("Analyte `$a` is not in the `analyte_map`."))
        end
        if length(conctable.sample) > 1
            length(level_map) == length(signaltable.sample) || throw(ArgumentError("The length of `level_map` is different from that of `table.sample`."))
            for a in conctable.analyte
                a in signaltable.analyte || throw(ArgumentError("Analyte `$a` is not in the `signatable`."))
            end
        end
        new{A, T}(analyte_map, signal, level_map, conctable, signaltable)
    end
    function MethodTable{T}(
                        analyte_map::Table,
                        signal::Symbol,
                        level_map::Vector{Int}, 
                        conctable::AbstractDataTable{B, <: T}, 
                        signaltable::Nothing) where {B, T}
        cp = propertynames(analyte_map)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analyte_map"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analyte_map"))
        :calibration in cp || throw(ArgumentError("Column `:calibration` is required in `analyte_map"))
        A = eltype(analyte_map.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable should be a subtype of $A"))
        for a in conctable.analyte
            a in analyte_map.analyte || throw(ArgumentError("Analyte `$a` is not in the `analyte_map`."))
        end
        new{A, T}(analyte_map, signal, level_map, conctable, signaltable)
    end
end

MethodTable(
            analyte_map::Table,
            signal::Symbol,
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{B, T}, 
            signaltable::Nothing) where {B, T} = MethodTable{T}(analyte_map, signal, level_map, conctable, signaltable)

MethodTable(
            analyte_map::Table,
            signal::Symbol,
            level_map::Vector{Int}, 
            conctable::AbstractDataTable{B, T}, 
            signaltable::AbstractDataTable{C, S}) where {B, C, T, S} = MethodTable{promote_type(T, S)}(analyte_map, signal, level_map, conctable, signaltable)

function getproperty(tbl::MethodTable, p::Symbol)
    if p == :analyte
        getfield(tbl, :analyte_map).analyte
    elseif p == :isd
        getfield(tbl, :analyte_map).analyte[getfield(tbl, :analyte_map).isd .< 0]
    elseif p == :nonisd
        getfield(tbl, :analyte_map).analyte[getfield(tbl, :analyte_map).isd .>= 0]
    elseif p == :point
        s = getfield(tbl, :signaltable)
        isnothing(s) ? s : s.sample
    elseif p == :level
        getfield(tbl, :conctable).sample
    else
        getfield(tbl, p)
    end
end

"""
    MethodTable(conctable::AbstractDataTable, signaltable::Union{AbstractDataTable, Nothing}, signal, level_map = []; kwargs...)
    MethodTable{T}(conctable::AbstractDataTable, signaltable::Union{AbstractDataTable, Nothing}, signal, level_map = []; kwargs...)

User-friendly contructors for `MethodTable`. `kwargs` will be columns in `analyte_map`; when `analyte`, `isd` and `calibration` are not provided, it will use analyte in `conctable`.
"""
MethodTable(conctable::AbstractDataTable{T}, 
            signaltable::AbstractDataTable{S},
            signal::Symbol,
            level_map::Vector{Int}; kwargs...) where {T, S} = MethodTable{promote_type(T, S)}(conctable, signaltable, signal, level_map; kwargs...)
MethodTable(conctable::AbstractDataTable{T}, 
            signaltable::Nothing,
            signal::Symbol,
            level_map::Vector{Int} = Int[]; kwargs...) where T = MethodTable{T}(conctable, signaltable, signal, level_map; kwargs...)
function MethodTable{T}(conctable::AbstractDataTable, 
                    signaltable::Union{AbstractDataTable, Nothing},
                    signal::Symbol,
                    level_map::Vector{Int}; kwargs...
            ) where T
    analyte_map = if isempty(kwargs)
        Table(; analyte = conctable.analyte)
    else
        Table(; kwargs...)
    end
    if !in(:analyte, propertynames(analyte_map))
        analyte_map = Table(analyte_map; analyte = conctable.analyte)
    end       
    if !in(:isd, propertynames(analyte_map))
        analyte_map = Table(analyte_map; isd = zeros(Int, length(conctable.analyte)))
    end
    if !in(:calibration, propertynames(analyte_map))
        analyte_map = Table(analyte_map; calibration = collect(eachindex(conctable.analyte)))
    end
    analyte_map = Table((; analyte = analyte_map.analyte, isd = analyte_map.isd, calibration = analyte_map.calibration), analyte_map)    
    MethodTable{T}(analyte_map, signal, level_map, conctable, signaltable)
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
* `calibration`: `Vector{MultipleCalibration{<: A}}` or `Vector{SingleCalibration{<: A}}`.
* `data`: Data for analysis, `AnalysisTable{A, <: T}` or `Nothing`.

# Properties
* `analyte`: `Vector{A}`, analytes in user-defined types, identical to `method.analyte_map.analyte`.
* `isd`: `Vector{<: A}`, analytes which are internal standards.
* `nonisd`: `Vector{<: A}`, analytes which are not internal standards.
* `point`: `Vector{Symbol}`, calibration points, identical to `method.signaltable.sample`.
* `level`: `Vector{Symbol}`, calibration levels, identical to `method.conctable.sample`.

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
        length(method.conctable.sample) > 1 ? map(method.conctable.analyte) do analyte
            calibration(method, analyte; type, zero, weight)
        end : map(method.conctable.analyte) do analyte
            SingleCalibration((analyte, ), first(getanalyte(method.conctable, analyte)))
        end
        ,
        data
    )
end

function getproperty(batch::Batch, p::Symbol)
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

include("tables.jl")
include("utils.jl")
include("cal.jl")
include("io.jl")

"""
    ColumnDataTable(tbl::RowDataTable{A, T}, sample_name::Symbol) where {A, T}
    ColumnDataTable{S}(tbl::RowDataTable{A, T}, sample_name::Symbol) where {A, S, T}

Convert `RowDataTable` to `ColumnDataTable` with `sample_name` as the column name of sample. `S` must be provided if `T` is not a valid constructor.
"""
ColumnDataTable(tbl::RowDataTable{A, T}, sample_name::Symbol) where {A, T} = ColumnDataTable{T}(tbl, sample_name)
ColumnDataTable{T}(tbl::RowDataTable{A}, sample_name::Symbol) where {A, T} = 
    ColumnDataTable{T}(tbl.analyte, sample_name, T((; (sample_name => tbl.sample, (tbl.analytename .=> getanalyte.(Ref(tbl), tbl.analyte))...)...)))

"""
    RowDataTable(tbl::RowDataTable{A, T}, sample_name::Symbol) where {A, T}
    RowDataTable{S}(tbl::RowDataTable{A, T}, sample_name::Symbol) where {A, S, T}

Convert `RowDataTable` to `ColumnDataTable` with `sample_name` as the column name of sample. `S` must be provided if `T` is not a valid constructor.
"""
RowDataTable(tbl::ColumnDataTable{A, T}, analyte_name::Symbol) where {A, T} = RowDataTable{T}(tbl, analyte_name)
RowDataTable{T}(tbl::ColumnDataTable{A}, analyte_name::Symbol) where {A, T} = 
    RowDataTable{T}(tbl.analyte, tbl.samplename, T((; (analyte_name => string.(tbl.analytename), (tbl.samplename .=> getsample.(Ref(tbl), tbl.samplename))...)...)))

end

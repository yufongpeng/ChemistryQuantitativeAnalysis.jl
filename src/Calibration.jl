module Calibration

using GLM, CSV, TypedTables, LinearAlgebra
export MultipleCalibration, SingleCalibration, 
    ColumnAnalysisTable, RowAnalysisTable, SampleWrapperTable, CalWrapperTable,
    Project, project, calibration,
    read_calibration, read_analysistable, read_project,
    cal_range, lloq, uloq, accuracy, accuracy!,
    inv_predict, inv_predict!, inv_predict_cal!, inv_predict_sample!, inv_predict_accuracy!, inv_predict_cal_accuracy!,
    find_analyte, get_analyte, find_sample, get_sample, find_isd, get_isd, switch_isd!,
    formula_repr, weight_repr, formula_repr_utf8, weight_repr_utf8, format_number

import Base: getproperty, show, write
    
abstract type AbstractCalibration{A, T} end
abstract type AbstractAnalysisTable{A, T} end
abstract type AbstractWrapperTable{A, T} end

"""
    ColumnAnalysisTable{A, T} <: AbstractAnalysisTable{A, T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples.

# Fields
* `sample_name`: `Symbol`, the column name that each element is sample name.
* `config`: `Table` with 3 columns, `analyte_name`, `analytes` and  `isd_map` which are properties of this type. To add new analytes, user should modify this table.
* `table`: Tabular data of type `T`.

# Properties
* `analyte_name`: `Vector{Symbol}`, the column names of `table` that are analytes names.
* `analytes`: `Vector{A}`, analytes in user-defined types. Its length and correponding analytes must matches those of `isd_map`.
* `isd_map`: `Vector{Int}`, each element is the index of the correponding internal standrads in `analytes` of the analyte in this position. `0` means no internal standard, and `-1` means the analyte itself is internal standard. For example, `[2, -1]` means that internal standard of the first analyte is the second analyte, and the second analyte is an internal standrad.
* `samples`: `Vector{Symbol}`, sample names.
* `isds`: `Vector{A}`, internal standards.
* `nonisds`: `Vector{A}`, analytes that are not internal standards.
"""
struct ColumnAnalysisTable{A, T} <: AbstractAnalysisTable{A, T}
    sample_name::Symbol
    config::Table{NamedTuple{(:analyte_name, :analytes, :isd_map), Tuple{Symbol, A, Int}}, 1, NamedTuple{(:analyte_name, :analytes, :isd_map), Tuple{Vector{Symbol}, Vector{A}, Vector{Int}}}}
    table::T
end

"""
    ColumnAnalysisTable(sample_name::Symbol, analyte_name::Vector{Symbol}, analytes::Vector{A}, isd_map::Vector{Int}, table::T) where {A, T}

A user-friendly contructor for `ColumnAnalysisTable`. See the documentation of the type for detail description.
"""
function ColumnAnalysisTable(sample_name::Symbol, analyte_name::Vector{Symbol}, analytes::Vector{A}, isd_map::Vector{Int}, table::T) where {A, T}
    length(analyte_name) == length(analytes) == length(isd_map) || throw(ArgumentError("Arguments `analyte_name`, `analytes` and `isd_map` should have the same length."))
    ColumnAnalysisTable(sample_name, Table(; analyte_name, analytes, isd_map), table)
end

function getproperty(tbl::ColumnAnalysisTable, property::Symbol)
    if property in [:analyte_name, :analytes, :isd_map] 
        getproperty(getfield(tbl, :config), property)
    elseif property == :samples
        Symbol.(getproperty(tbl.table, tbl.sample_name))
    elseif property == :isds
        getfield(tbl, :config).analytes[getfield(tbl, :config).isd_map .< 0]
    elseif property == :nonisds
        getfield(tbl, :config).analytes[getfield(tbl, :config).isd_map .>= 0]
    else
        getfield(tbl, property)
    end
end

"""
    RowAnalysisTable{T} <: AbstractAnalysisTable{T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples.

# Fields
* `sample_name`: `Vector{Symbol}`, the column names that are sample names.
* `analyte_name`: `Symbol`, the column name that each element is analyte name.
* `analytes`: `Vector{A}`, analytes in user-defined types. Its length and correponding analytes must matches those of `isd_map`.
* `isd_map`: `Symbol` or `Nothing`, the column name that contains internal standard information. `nothing` indicates no internal standards.
* `table`: Tabular data of type `T`.

# Properties
* `samples`: `Vector{Symbol}`, sample names.
* `isds`: `Vector{A}`, internal standards.
* `nonisds`: `Vector{A}`, analytes that are not internal standards.
"""
struct RowAnalysisTable{A, T} <: AbstractAnalysisTable{A, T}
    sample_name::Vector{Symbol}
    analyte_name::Symbol
    analytes::Vector{A}
    isd_map::Union{Symbol, Nothing}
    table::T
    function RowAnalysisTable(sample_name::Vector{Symbol}, analyte_name::Symbol, analytes::Vector{A}, isd_map::Union{Symbol, Nothing}, table::T) where {A, T}
        length(analytes) == length(getproperty(table, analyte_name)) || throw(ArgumentError("The length of `analytes` is different from that of `table`."))
        new{A, T}(sample_name, analyte_name, analytes, isd_map, table)
    end
end

function getproperty(tbl::RowAnalysisTable{A}, property::Symbol) where A
    if property == :samples
        getfield(tbl, :sample_name)
    elseif property == :isds
        isd_map = getfield(tbl, :isd_map)
        isnothing(isd_map) && return A[]
        getproperty(getfield(tbl, :table), getfield(tbl, :analyte_name))[getproperty(getfield(tbl, :table), isd_map) .< 0]
    elseif property == :nonisds
        isd_map = getfield(tbl, :isd_map)
        isnothing(isd_map) && return getproperty(getfield(tbl, :table), getfield(tbl, :analyte_name))
        getproperty(getfield(tbl, :table), getfield(tbl, :analyte_name))[getproperty(getfield(tbl, :table), isd_map) .>= 0]
    else
        getfield(tbl, property)
    end
end

"""
    SampleWrapperTable{A, T} <: AbstractWrapperTable{A, T}

Tabular data wrapper for mapping analytes in `analysistable` to calibration curves.

# Fields
* `cal_map`: The length and corresponding analytes matches analyte in `analysistable.analytes`, and each element is the index of another analyte that its calibration curve is used for quantification of this analyte. For example, `[1, -1, 1]` means that the first analyte uses calibration curve of itself, the second analyte is internal standard, and the third analyte also uses calibration curve of the first analyte. 
* `analysistable`: `AbstractAnalysisTable{A, T}`

# Properties
All properties for `AbstractAnalysisTable{A, T}` are available.
"""
struct SampleWrapperTable{A, T} <: AbstractWrapperTable{A, T}
    cal_map::Vector{Int}
    analysistable::AbstractAnalysisTable{A, T}
    function SampleWrapperTable(cal_map::Vector{Int}, table::AbstractAnalysisTable{A, T}) where {A, T}
        length(cal_map) == length(table.analytes) || throw(ArgumentError("The length of `cal_map` is different from that of `table`."))
        new{A, T}(cal_map, table)
    end
end

function getproperty(tbl::SampleWrapperTable, property::Symbol)
    if property == :cal_map
        getfield(tbl, :cal_map)
    elseif property == :analysistable
        getfield(tbl, :analysistable)
    else
        getproperty(getfield(tbl, :analysistable), property)
    end
end

"""
    CalWrapperTable{A, T} <: AbstractWrapperTable{A, T}

Tabular data wrapper for mapping levels to points.

# Fields
* `level_map`: `Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.
* `conctable`: `AbstractAnalysisTable{A, T}` containing concentration data for each level, and ISD information. Sample names must be symbol or string of integers for multiple levels. One level indicates using `SingleCalibration`.
* `signaltable`: `AbstractAnalysisTable{A}` containig signal for each point. It can be `nothing` if signal data is unecessary.
"""
struct CalWrapperTable{A, T} <: AbstractAnalysisTable{A, T}
    level_map::Vector{Int}
    conctable::AbstractAnalysisTable{A}
    signaltable::Union{AbstractAnalysisTable{A}, Nothing}
    function CalWrapperTable(level_map::Vector{Int}, conctable::AbstractAnalysisTable{A, T}, signaltable::Union{AbstractAnalysisTable{A}, Nothing}) where {A, T}
        if length(conctable.samples) > 1
            length(level_map) == length(signaltable.samples) || throw(ArgumentError("The length of `level_map` is different from that of `table.samples`."))
            for a in conctable.analytes
                a in signaltable.analytes || throw(ArgumentError("Analyte `$a` is not in the `signatable`."))
            end
        end
        new{A, T}(level_map, conctable, signaltable)
    end
end

"""
    MultipleCalibration{A, T} <: AbstractCalibration{A, T}

A type holding all data for a calibration curve.

# Fields
* `analyte`: `Int` indicates the index the analyte in `caltable.conctable.analytes`: 
* `isd`: `Int` indicates the index internal standard in `caltable.conctable.analytes`. `0` indicates no internal standard.
* `type`: `Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`).
* `zero`: `Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`).
* `weight`: `Float64` represents the exponential applying to each element of `x` as a weighting vector.
* `formula`: `FormulaTerm`, the formula for fitting calibration curve.
* `caltable`: `CalAnalysisTable{A, T}`, the original calibration data.
* `table`: `TypedTable.Table`, the clean up calibration data, containing 7 columns.
    * `id`: Sample name
    * `level`: The index of concentration level. The larger, the higher concentraion it represents.
    * `y`: Signal
    * `x`: Concentraion
    * `x̂`: Predicted concentration
    * `accuracy`: Accuracy, i.e. `x̂/x`.
    * `include`: Whether this point is included or not
* `model`: `GLM` object.
"""
mutable struct MultipleCalibration{A, T} <: AbstractCalibration{A, T}
    analyte::Int
    isd::Int
    type::Bool
    zero::Bool
    weight::Float64
    formula::FormulaTerm
    caltable::CalWrapperTable{A, T}
    table::Table
    model
end

"""
    SingleCalibration{A, T} <: AbstractCalibration{A, T}

A type holding all data for single point calibration.

# Fields
* `analyte`: `Int` indicates the index the analyte in `caltable.conctable.analytes`.
* `isd`: `Int` indicates the index internal standard in `caltable.conctable.analytes`. `0` indicates no internal standard.
* `caltable`: `CalWrapperTable{A, T}`, the original calibration data.
"""
mutable struct SingleCalibration{A, T} <: AbstractCalibration{A, T}
    analyte::Int
    isd::Int
    caltable::CalWrapperTable{A, T}
end

"""
    Project{A, T}

A type for holding all data.

# Fields
* `calibration`: `Vector{MultipleCalibration{A, T}}` or `Vector{SingleCalibration{A, T}}`.
* `caltable`: `CalWrapperTable{A, T}`, the original calibration data.
* `sampletable`: Sample data, `SampleWrapperTable{A}` or `Nothing`.
* `resulttable`: Result data, `SampleWrapperTable{A}` or `Nothing`.
"""
mutable struct Project{A, T}
    calibration::Vector{<: AbstractCalibration{A, T}}
    caltable::CalWrapperTable{A, T}
    sampletable::Union{SampleWrapperTable{A}, Nothing}
    resulttable::Union{SampleWrapperTable{A}, Nothing}
end

include("utils.jl")
include("cal.jl")
include("io.jl")

end

module Calibration

using GLM, CSV, TypedTables, LinearAlgebra
export AbstractCalibration, MultipleCalibration, SingleCalibration, 
    AbstractAnalysisTable, ColumnAnalysisTable, RowAnalysisTable,
    Project, project, calibration,
    read_calibration, read_analysistable, read_project,
    cal_range, lloq, uloq, accuracy, accuracy!,
    inv_predict, inv_predict_cal!, inv_predict_accuracy!,
    formula_repr, weight_repr, formula_repr_utf8, weight_repr_utf8, format_number

import Base: getproperty, show
    
abstract type AbstractCalibration{T} end
abstract type AbstractAnalysisTable{T} end

mutable struct MultipleCalibration{T} <: AbstractCalibration{T}
    analyte::Int
    isd::Int
    type::Bool
    zero::Bool
    weight::Float64
    formula::FormulaTerm
    source::AbstractAnalysisTable{T}
    table::Table
    model
end

mutable struct SingleCalibration{T} <: AbstractCalibration{T}
    analyte::Int
    isd::Int
    source::AbstractAnalysisTable{T}
    β::Float64
end

struct ColumnAnalysisTable{T} <: AbstractAnalysisTable{T}
    sample_name::Symbol
    analyte_name::Vector{Symbol}
    analytes::Vector{T}
    table
end
getproperty(tbl::ColumnAnalysisTable, property::Symbol) = property == :samples ? Symbol.(getproperty(tbl.table, tbl.sample_name)) : getfield(tbl, property)

struct RowAnalysisTable{T} <: AbstractAnalysisTable{T}
    sample_name::Vector{Symbol}
    analyte_name::Symbol
    analytes::Vector{T}
    table
end
getproperty(tbl::RowAnalysisTable, property::Symbol) = property == :samples ? getfield(tbl, :sample_name) : getfield(tbl, property)

mutable struct Project
    calibration::Vector{<: AbstractCalibration}
    sample::Union{AbstractAnalysisTable, Nothing}
end

include("utils.jl")
include("cal.jl")
include("io.jl")

end

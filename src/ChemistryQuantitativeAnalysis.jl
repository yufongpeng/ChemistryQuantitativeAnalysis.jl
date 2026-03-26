module ChemistryQuantitativeAnalysis

using GLM, LsqFit, CSV, TypedTables, LinearAlgebra, Dictionaries, ThreadsX, Tables, Pkg, Statistics, StatsAPI
using LsqFit: LsqFitResult

export 
    # DataType
    SampleDataTable, AnalyteDataTable, AnalysisTable, analysistable, 
    # Method/Batch
    AnalysisMethod, Batch, 
    # CalibrationModel
    CalibrationModel,
    ProportionalCalibrator,
    LinearCalibrator,
    QuadraticProportionalCalibrator,
    QuadraticCalibrator,
    LogarithmicCalibrator,
    ExponentialCalibrator,
    PowerCalibrator,
    # Calibrator
    ExternalCalibrator, InternalCalibrator, 
    # Weight
    ConstWeight,
    RootXWeight,
    RootYWeight,
    RootXYWeight,
    XWeight,
    YWeight,
    XYWeight,
    SqXWeight,
    SqYWeight,
    SqXYWeight,
    RootLogXWeight,
    RootLogYWeight,
    LogXWeight,
    LogYWeight,
    SqLogXWeight,
    SqLogYWeight,
    RootExpXWeight,
    RootExpYWeight,
    ExpXWeight,
    ExpYWeight,
    SqExpXWeight,
    SqExpYWeight,

    # Workflow1: edit method
    edit_method!, 
    assign_isd!, assign_std!, replace_std!, 
    # Workflow2: calibrate
    calibrate, calibrate!, 
    edit_method_calibrate!, 
    model_calibrator!, load_method!, assign_isd_calibrate!, assign_std_calibrate!, replace_std_calibrate!,
    # Workflow3: quantify
    relative_signal, set_relative_signal, set_relative_signal!, quantify_relative_signal!,
    inv_predict, set_inv_predict, set_inv_predict!, quantify_inv_predict!,
    quantify, set_quantify, set_quantify!, quantify!,
    accuracy, set_accuracy, set_accuracy!, validate!,
    analyze!,

    # table attr
    analyteobj, sampleobj, analytename, samplename, analytecol, samplecol, 
    # analyte related 
    stdof, isdof, isstd, isisd, 
    # finder/getter
    findanalyte, getanalyte, findsample, getsample, findcalibrator, getcalibrator, eachanalyte, eachsample,
    # analytical attr
    dynamic_range, lloq, uloq, signal_range, signal_lloq, signal_uloq, signal_lob, signal_lod, signal_loq, 
    
    # batch construction utils
    mkbatch, 
    # table construction utils
    typedmap

import Base: getproperty, propertynames, show, write, eltype, length, iterate, 
        getindex, setindex!, insert!, get!, delete!, get, 
        pairs, keys, values, haskey, isassigned, 
        copy
import Dictionaries: set!, unset!, isinsertable, issettable
import Tables: istable, rowaccess, rows, columnaccess, columns
import StatsAPI: coef, predict, r2

@deprecate AbstractCalibration AbstractCalibrator
abstract type AbstractCalibrator{A, N} end
abstract type AbstractDataTable{A, S, N, T} end
# abstract type AbstractAnalysisTable{A, T, S} end
const TypeOrFn = Union{<: Type, <: Function}

include("table.jl")
include("method.jl")
include("weight.jl")
include("calibrator.jl")
include("batch.jl")

include("interface.jl")
include("utils.jl")
include("internal_utils.jl")
include("calibration.jl")
include("quant.jl")
include("input.jl")
include("output.jl")
include("additional_constructors.jl")

end
abstract type AbstractCalibrationModel{T} end 

cdoc = """
    CalibrationModel{T <: CurveType} 

    const ProportionalCalibrator = CalibrationModel{Proportional}
    const LinearCalibrator = CalibrationModel{Linear}
    const QuadraticProportionalCalibrator = CalibrationModel{QuadraticProportional}
    const QuadraticCalibrator = CalibrationModel{Quadratic}
    const LogarithmicCalibrator = CalibrationModel{Logarithmic}
    const ExponentialCalibrator = CalibrationModel{Exponential}
    const PowerCalibrator = CalibrationModel{Power}

Calibration model type of curve type `T`. 

# Fields
* `weight::ComposedWeight`: weight object to generate weight function. 
"""
@doc cdoc 
struct CalibrationModel{T <: CurveType} <: AbstractCalibrationModel{T}
    weight::ComposedWeight 
end

"""
    LsqFitMachine

Type wrapping `LsqFitResult` and the target function.

# Fields 
* `fit::LsqFitResult`
* `fn::Function`: target function.
* `totalvariance::Real`: variance of y.
"""
struct LsqFitMachine
    fit::LsqFitResult
    fn::Function
    totalvariance::Real
end

@doc cdoc 
const ProportionalCalibrator = CalibrationModel{Proportional}
@doc cdoc 
const LinearCalibrator = CalibrationModel{Linear}
@doc cdoc 
const QuadraticProportionalCalibrator = CalibrationModel{QuadraticProportional}
@doc cdoc 
const QuadraticCalibrator = CalibrationModel{Quadratic}
@doc cdoc 
const LogarithmicCalibrator = CalibrationModel{Logarithmic}
@doc cdoc 
const ExponentialCalibrator = CalibrationModel{Exponential}
@doc cdoc 
const PowerCalibrator = CalibrationModel{Power}

@deprecate MultipleCalibration ExternalCalibrator
"""
    ExternalCalibrator{A, N, T <: Table} <: AbstractCalibrator{A, N}

A mutable type holding all data and settings for an external calibration curve. `A` determines analyte type, and `N` determines numeric type. 

# Fields
* `analyte::A`: analyte being quantified.
* `isd`: internal standard for which `nothing` indicates no internal standard.
* `table::TypedTable.Table`, the cleaned up calibration data, containing 7 columns.
    * `id`: Point name.
    * `level`: The index of concentration level. The larger, the higher concentraion it represents.
    * `y`: Signal or relative signal
    * `x`: True concentraion
    * `x̂`: Predicted concentration
    * `accuracy`: Accuracy, i.e. `x̂/x`.
    * `include`: Whether this point is included or not.
* `model::AbstractCalibrationModel`: calibration model.
* `machine`: calibration machine. Either a `LsqFitMachine` or fitted `GLM` object.

!!! note 
    For `N`, only `Float64` is supported because of issue of `GLM`.
"""
mutable struct ExternalCalibrator{A, N, T <: Table} <: AbstractCalibrator{A, N}
    analyte::A
    isd
    table::T
    model::AbstractCalibrationModel
    machine
    function ExternalCalibrator(analyte::A, isd, table::T, model, machine) where {A, T}
        N = eltype(table.x)
        N == eltype(table.y) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column y is `", eltype(table.y), "`.")))
        N == eltype(table.x̂) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column x̂ is `", eltype(table.x̂), "`.")))
        N == eltype(table.accuracy) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column accuracy is `", eltype(table.accuracy), "`.")))
        new{A, N, T}(analyte, isd, table, model, machine)
    end
end

@deprecate SingleCalibration InternalCalibrator
"""
    InternalCalibrator{A, N} <: AbstractCalibrator{A, N}

A mutable type holding all data for internal calibration with a single level.

# Fields
* `analyte::A`: the analyte with known concentration (internal standard).
* `conc::Float64`: concentration of analyte.
"""
mutable struct InternalCalibrator{A, N} <: AbstractCalibrator{A, N}
    analyte::A
    conc::N
end
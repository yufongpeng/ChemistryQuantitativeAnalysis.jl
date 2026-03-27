"""
    AnalysisMethod{A, M <: Table, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}

A type containing analytes settings, and source calibration data. `A` determines analyte type.

# Fields
* `analytetable::M` contaning three columns.
    * `analyte::AbstractVector{A}`: analytes in user-defined types.
    * `isd::AbstractVector{Int}: index of internal standard. `0` means no internal standard, and `-1` means the analyte itself is a internal standard.
    * `std`: index of analyte as calibration standard (internal or external). `0` means no standard.
    * `model`: calibration model type.
* `signal::Symbol`: type and key name of experimental acquisition data, e.g. `:area`.
* `rel_sig::Symbol`: key name of relative signal.
* `est_conc::Symbol`: key name of estimated concentration.
* `nom_conc::Symbol`: key name of nominal concentration.
* `acc::Symbol`: key name of accuracy.
* `pointlevel::Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.
* `conctable::C`: containing concentration data for each level. Samples must integers. 
* `signaltable::D` containig signal for each point. It can be `nothing` if signal data is unecessary.

# Properties
* `analyte::AbstractVector{A}`: analytes in user-defined types, identical to `analytetable.analyte`.
* `isd::AbstractVector{A}`: analytes which are internal standards.
* `nonisd::AbstractVector{A}`: analytes which are not internal standards.
* `std::AbstractVector{A}`: analytes which are calibration standards.
* `point::AbstractVector{S}`: calibration points, identical to `sampleobj(signaltable)`. If `signaltable` is `nothing`, this value is `nothing` too.
* `level::AbstractVector{Int}`: calibration levels, identical to `sampleobj(conctable)`.
"""
struct AnalysisMethod{A, M <: Table, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}
    analytetable::M
    signal::Symbol
    rel_sig::Symbol
    est_conc::Symbol
    nom_conc::Symbol
    acc::Symbol
    pointlevel::Vector{Int}
    conctable::C
    signaltable::D
    function AnalysisMethod(
                    analytetable::M,
                    signal::Symbol,
                    rel_sig::Symbol,
                    est_conc::Symbol,
                    nom_conc::Symbol,
                    acc::Symbol,
                    pointlevel::AbstractVector{Int}, 
                    conctable::AbstractDataTable{B, Int}, 
                    signaltable::AbstractDataTable{C, S}
                ) where {B, M <: Table, C, S}
        cp = propertynames(analytetable)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analytetable"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analytetable"))
        :std in cp || throw(ArgumentError("Column `:std` is required in `analytetable"))
        if !in(:model, cp)
            analytetable = Table(analytetable; model = Any[mkcalmodel(CalibrationModel{Linear})  for _ in eachindex(analytetable)])
        end
        A = eltype(analytetable.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable $B should be a subtype of $A"))
        C <: A || throw(ArgumentError("Analyte type of signaltable $C should be a subtype of $A"))
        for a in analyteobj(conctable)
            a in analytetable.analyte || throw(ArgumentError("Analyte `$a` is not in the `analytetable`."))
        end
        if length(sampleobj(conctable)) > 1
            length(pointlevel) == length(sampleobj(signaltable)) || throw(ArgumentError("The length of `pointlevel` is different from that of `sampleobj(signaltable)`."))
            for a in analyteobj(conctable)
                a in analyteobj(signaltable) || throw(ArgumentError("Analyte `$a` is not in the `signatable`."))
            end
        end
        new{A, typeof(analytetable), typeof(conctable), typeof(signaltable)}(analytetable, signal, rel_sig, est_conc, nom_conc, acc, convert(Vector{Int}, pointlevel), conctable, signaltable)
    end
    function AnalysisMethod(
                    analytetable::M,
                    signal::Symbol,
                    rel_sig::Symbol,
                    est_conc::Symbol,
                    nom_conc::Symbol,
                    acc::Symbol,
                    pointlevel::AbstractVector{Int}, 
                    conctable::AbstractDataTable{B, Int}, 
                    signaltable::Nothing
                ) where {B, M <: Table}
        cp = propertynames(analytetable)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analytetable"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analytetable"))
        :std in cp || throw(ArgumentError("Column `:std` is required in `analytetable"))
        if !in(:model, cp)
            analytetable = Table(analytetable; model = Any[mkcalmodel(CalibrationModel{Linear}) for _ in eachindex(analytetable)])
        end
        A = eltype(analytetable.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable should be a subtype of $A"))
        for a in analyteobj(conctable)
            a in analytetable.analyte || throw(ArgumentError("Analyte `$a` is not in the `analytetable`."))
        end
        new{A, typeof(analytetable), typeof(conctable), Nothing}(analytetable, signal, rel_sig, est_conc, nom_conc, acc, convert(Vector{Int}, pointlevel), conctable, signaltable)
    end
end

"""
    AnalysisMethod(
        conctable::AbstractDataTable, 
        signaltable::Union{AbstractDataTable, Nothing}, 
        signal, 
        pointlevel = Int[]; 
        signal = :area, 
        rel_sig = :relative_signal, 
        est_conc = :estimated_concentration, 
        nom_conc = :nominal_concentration, 
        kwargs...
    )

Convenient contructors for `AnalysisMethod` which make constructing a `Table` optional.

`signal` data type for quantification. `levelname` is the column name for `pointlevel`. See `AnalysisMethod` for details of other arguments.

# Keyword arguments
* `rel_sig`: key name of relative signal. 
* `est_conc`: key name of estimated concentration. 
* `nom_conc`: key name of nominal concentration. 
* `acc`: key name of accuracy. 
* Other keyword arguments will be columns in `analytetable`; when `analyte`, `isd` and `std` are not provided, it will use analyte in `conctable`. 
"""
function AnalysisMethod(
        conctable::AbstractDataTable, 
        signaltable::Union{AbstractDataTable, Nothing},
        signal,
        pointlevel = nothing; 
        rel_sig = :relative_signal,
        est_conc = :estimated_concentration,
        nom_conc = :nominal_concentration, 
        acc = :accuracy, 
        kwargs...
    )
    signal = Symbol(signal)
    rel_sig = Symbol(rel_sig)
    est_conc = Symbol(est_conc)
    nom_conc = Symbol(nom_conc)
    acc = Symbol(acc)
    if !isa(pointlevel, AbstractVector{<: Integer})
        pointlevel = getproperty(signaltable, Symbol(pointlevel))
    end
    analyte = get(kwargs, :analyte, analyteobj(conctable))
    isd = get(kwargs, :isd, zeros(Int, length(analyteobj(conctable))))
    std = get(kwargs, :std, collect(eachindex(analyteobj(conctable))))
    model = get(kwargs, :model, 
        Any[mkcalmodel(CalibrationModel{Linear}) for i in eachindex(std)]
    )
    analytetable = Table(; analyte, isd, std, model, kwargs...)    
    AnalysisMethod(analytetable, signal, rel_sig, est_conc, nom_conc, acc, pointlevel, conctable, signaltable)
end
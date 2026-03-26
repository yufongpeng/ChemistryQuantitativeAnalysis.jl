"""
    Batch{A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibrator{<: A}}, D <: Union{AnalysisTable{<: A}, Nothing}}

A type representing a batch for quantitative analysis. `A` determines analyte type, and `T` determines table type.

# Fields
* `method::M`: method.
* `calibrator::C`:, calibrators.
* `data::D`:, data for analysis.

# Properties
* `analyte::AbstractVector{A}`: analytes in user-defined types, identical to `method.analyte`.
* `isd::AbstractVector{<: A}`: analytes which are internal standards, identical to `method.isd`.
* `nonisd::AbstractVector{<: A}`: analytes which are not internal standards, identical to `method.nonisd`.
* `std::AbstractVector{A}: analytes which are calibration standards, identical to `method.std`.
* `point::Union{AbstractVector, Nothing}`: calibration points, identical to `method.point`.
* `level::AbstractVector{Int}`: calibration levels, identical to `method.level`.

# Constructors
* `Batch(method, calibrator, data = nothing)`
"""
struct Batch{A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibrator}, D <: Union{AnalysisTable{<: A}, Nothing}}
    method::M
    calibrator::C
    data::D
    Batch(method::M, calibrator::C, data::D) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibrator}, D <: AnalysisTable{<: A}} = new{A, M, C, D}(method, calibrator, data)
    Batch(method::M, calibrator::C, data::Nothing) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibrator}} = new{A, M, C, Nothing}(method, calibrator, data)
    Batch(method::M, calibrator::C) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibrator}} = new{A, M, C, Nothing}(method, calibrator, nothing)
end

"""
    Batch(method::AnalysisMethod, data = nothing)

Construct a `Batch` from `method`, and optionally `data` with specified calibrator parameters. See "ExternalCalibrator" for detail description of keyword arguments.
"""
function Batch(method::AnalysisMethod{A}, data = nothing) where A
    Batch(
        method,
        (AbstractCalibrator{T} where {T <: A})[],
        # map(analyteobj(method.conctable)) do analyte
        #     calibrate(method, analyte; model, kwargs...)
        # end,
        data
    )
end
"""
    Batch(batch::Batch, at::AnalysisTable)

Construct a new batch from an old batch and a `AnalysisTable` as new data.
"""
Batch(batch::Batch, at::AnalysisTable) = Batch(batch.method, batch.calibrator, at)
"""
    Batch(dt::AbstractDataTable; 
        signal = :area, 
        rel_sig = :relative_signal, 
        est_conc = :estimated_concentration, 
        nom_conc = :nominal_concentration, 
        acc = :accuracy, 
        calid = r"Cal_(\\d)_(\\d*-*\\d*)", 
        order = "LR", 
        ratio = nothing, 
        dilution_factor = nothing, 
        conc_factor = 1,
        parse_decimal = x -> replace(x, "-" => "."))

Construct a batch from data. All analytes are considered as normal analytes, so calibration curves are not contructed immediately; call `calibrate!` after editing and confirming calibration settings in `batch.method.analytetable` using `edit_method!`, `assign_std!` and `assign_isd!`.

# Arguments
* `dt`: data containing both calibration curves and samples.

# Keyword Arguments
* `signal`: type and key name for storing experimental acquisition data, e.g. `:area`.
* `rel_sig`: key name for storing relative signal.
* `est_conc`: key name for storing estimated concentration.
* `nom_conc`: key name for storing nominal concentration.
* `acc`: key name for storing accuracy.
* `calid`: identifier for calibration point. It has two possible values.
    * `Regex`, level and concentration related factors (ratio of concentration or dilution factor) should be captured. The former should be able to be parsed as integer directly; the latter should be able to be parsed as floating number after applying `parse_decimal`. 
* `order::String`: represents the order and identity of captured string; `L` is level, `R` is ratio of concentration, `D` is dilution factor.
* `level`: level of calibration of each sample (`missing` for non calibration sample). If it is a property name of `dt`, the column is used. `nothing` indicates using captured values. 
* `ratio`: ratio of concantrations of each level. If it is a property name of `dt`, the column is used. `nothing` indicates using captured values. 
* `dilution_factor`: dilution factors of each level. If it is a property name of `dt`, the column is used. `nothing` indicates using captured values. 
* `conc_factor::Union{Number, Vector{<: Number}`: concentration equals to conc_factor * ratio or conc_factor / dilution_factor. When a vector is provided, each element represents `conc_factor` value of each analyte. If it is a property name of `dt`, the column is used. 
* `parse_decimal::Function`: converts a string into another string which can be parsed as floating number.
"""
function Batch(dt::AbstractDataTable; 
                signal = :area, 
                rel_sig = :relative_signal, 
                est_conc = :estimated_concentration, 
                nom_conc = :nominal_concentration, 
                acc = :accuracy, 
                calid = r"Cal_(\d)_(\d*-*\d*)", 
                order = "LR", 
                level = nothing,
                ratio = nothing, 
                dilution_factor = nothing, 
                conc_factor = 1,
                parse_decimal = x -> replace(x, "-" => "."))
    signal = Symbol(signal)
    rel_sig = Symbol(rel_sig)
    est_conc = Symbol(est_conc)
    nom_conc = Symbol(nom_conc)
    acc = Symbol(acc)
    if dt isa SampleDataTable 
        if level in propertynames(dt)
            level = getproperty(dt, level)
        end
        if ratio in propertynames(dt)
            ratio = unique!(filter(!ismissing, getproperty(dt, ratio)))
        end
        if dilution_factor in propertynames(dt)
            dilution_factor = unique!(filter(!ismissing, getproperty(dt, dilution_factor)))
        end
    else 
        if conc_factor in propertynames(dt)
            conc_factor = getproperty(dt, conc_factor)
        end
    end
    idc, pointlevel, levels, concs = parse_calibration_level_name(dt, calid, order, level, ratio, dilution_factor, conc_factor, parse_decimal)
    aj = analyteobj(dt)
    id = idcol(dt)
    sj = sampleobj(dt)
    ids = setdiff(eachindex(sj), idc)
    tbl = Table(table(dt))
    if dt isa SampleDataTable
        sampledata = SampleDataTable(aj, sj[ids], id, tbl[ids])
        signaltable = SampleDataTable(aj, sj[idc], id, tbl[idc])
    else
        sampledata = AnalyteDataTable(aj, sj[ids], id, getproperties(tbl, tuple(samplename(dt)[ids]...)))
        signaltable = AnalyteDataTable(aj, sj[idc], id, getproperties(tbl, tuple(samplename(dt)[idc]...)))
    end
    if dt isa SampleDataTable
        conctable = SampleDataTable(Table(; Level = levels, (analytename(dt) .=> (collect(i) for i in zip(concs...)))...), :Level)
    else
        conctable = AnalyteDataTable(Table(; (id => aj, )..., (Symbol.(levels) .=> map(x -> repeat([0], length(aj)) .+ x, concs))...), id, Int)
    end
    analytetable = Table(; analyte = aj, isd = zeros(Int, length(aj)), std = collect(1:length(aj)))
    Batch(AnalysisMethod(analytetable, signal, rel_sig, est_conc, nom_conc, acc, pointlevel, conctable, signaltable), ExternalCalibrator{<: eltype(analyteobj(dt))}[], analysistable((signal => sampledata, )))
end
"""
    quantify_calibrator!(batch::Batch)
    quantify_calibrator!(cal::Vector{<: AbstractCalibrator})
    quantify_calibrator!(cal::InternalCalibrator)
    quantify_calibrator!(cal::ExternalCalibrator)

Inverse predict concentration, update each `cal.table.x̂` with the result(s) and returns `cal` or `batch`.
"""
quantify_calibrator!(batch::Batch) = (quantify_calibrator!(batch.calibrator); batch)
quantify_calibrator!(cal::Vector{<: AbstractCalibrator}) = (foreach(quantify_calibrator!, cal); cal)
quantify_calibrator!(cal::InternalCalibrator) = cal
function quantify_calibrator!(cal::ExternalCalibrator)
    cal.table.x̂ .= inv_predict(cal, cal.table.y)
    cal
end

function fill_result!(dt::SampleDataTable, result)
    for (a, c) in zip(eachanalyte(dt), result)
        a .= c
    end
    dt
end
function fill_result!(dt::AnalyteDataTable, result)
    for (i, p) in enumerate(eachsample(dt))
        p .= getindex.(result, i)
    end
    dt
end 

"""
    validate_calibrator!(batch::Batch)
    validate_calibrator!(cal::Vector{<: AbstractCalibrator})
    validate_calibrator!(cal::ExternalCalibrator)
    validate_calibrator!(cal::InternalCalibrator)

Calculate accuracy for analyte specified by `cal`, update each `cal.table.x̂` with the result(s), and return `cal` or `batch`.
"""
validate_calibrator!(batch::Batch) = (validate_calibrator!(batch.calibrator); batch)
validate_calibrator!(cal::Vector{<: AbstractCalibrator}) = (foreach(validate_calibrator!, cal); cal)
function validate_calibrator!(cal::ExternalCalibrator)
    cal.table.accuracy .= accuracy(cal.table.x̂, cal.table.x)
    cal
end
validate_calibrator!(cal::InternalCalibrator) = cal

"""
    const validate_quantify_calibrator! = validate_calibrator! ∘ quantify_calibrator!

Apply `quantify_calibrator!` and `validate_calibrator!` subsequantly to `Batch` or `AbstractCalibrator`.
"""
const validate_quantify_calibrator! = validate_calibrator! ∘ quantify_calibrator!

"""
    calibration(anisd::Tuple, tbl::Table; analyte = 1, isd = 0, type = true, zero = false, weight = 0)
    calibration(batch::Batch{A}, ana::B; id = sample(batch.method.signaltable), isd = isdof(batch.method, analyte),
                        type = true, zero = false, weight = 0) where {A, B <: A}
    calibration(method::AnalysisMethod{A}, ana::B; id = sample(method.signaltable), isd = isdof(method, analyte),
                type = true, zero = false, weight = 0) where {A, B <: A}

Create `ExternalCalibrator`.

* `anisd`: `Tuple{B <: A, Any}` which will be stored as `analyte`. The first element is the analyte to be quantified, and the second element is its internal standard, which `nothing` means no internal standard.
* `analyte`: `B` which will be stored as first argument of `analyte`.
* `tbl`: `TypedTable.Table` which will be stored as `table`. It is clean-up calibration data for points selection.
* `method`: `AnalysisMethod`, calibration method and data. 
* `batch`: `Batch`.

# Keyword arguments
* `type`: `Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`), which will be stored as `type`.
* `zero`: `Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`), which will be stored as `zero`.
* `weight`: `Float64` represents the exponential applying to each element of `x` as a weighting vector, which will be stored as `weight`.
"""
# recalibration after point select
function externalcalibrate(analyte::A, isd, tbl::Table; model = CalibrationModel{Linear}, kwargs...) where A
    id = findall(x -> isa(x, Number), tbl.y)
    tbl = tbl[id]
    xlevel = filter(>(0), unique(tbl.x))
    table = :id in propertynames(tbl) ? tbl : Table((; id = collect(1:length(tbl)), ), tbl)
    level = map(table.x) do i 
        j = findfirst(x -> i == x, xlevel)
        isnothing(j) ? 0 : j 
    end
    table = :level in propertynames(tbl) ? table : Table(table; level)
    table = :include in propertynames(tbl) ? table : Table(table; include = level .> 0)
    calmodel = mkcalmodel(model; kwargs...)
    calmachine = mkcalmachine(calmodel, table)
    table = Table(; id = table.id, level = table.level, y = table.y, x = table.x, x̂ = zeros(eltype(table.x), length(table)), accuracy = zeros(eltype(table.x), length(table)), include = table.include)
    validate_quantify_calibrator!(ExternalCalibrator(analyte, isd, table, calmodel, calmachine))
end

"""
init_calibration
"""
init_calibrate(batch::Batch{A}, analyte::B; 
            model = CalibrationModel{Linear}, 
            kwargs...
            ) where {A, B} = init_calibrate(batch.method, analyte; model, kwargs...)
function init_calibrate(method::AnalysisMethod{A}, analyte::B; 
                            model = CalibrationModel{Linear}, 
                            kwargs...
                            ) where {A, B <: A}

    isd = isdof(analyte, method)
    isd == analyte && return InternalCalibrator(analyte, first(getanalyte(method.conctable, analyte)))
    ord = sortperm(method.pointlevel)
    level = method.pointlevel[ord]
    conc = getanalyte(method.conctable, analyte)
    ya = getanalyte(method.signaltable, analyte)
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    id = sampleobj(method.signaltable)
    # Use Float64, bug of GLM
    numbertype = Float64 # eltype(conctable)
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        convert(numbertype, conc[findsample(method.conctable, Symbol(l))])
                    end, 
                    y = convert(Vector{numbertype}, y), 
                    x̂ = zeros(numbertype, length(id)), 
                    accuracy = zeros(numbertype, length(id)),
                    include = level .> 0
                    )
    calmodel = mkcalmodel(model; kwargs...)
    calmachine = mkcalmachine(calmodel, table)
    validate_quantify_calibrator!(ExternalCalibrator(analyte, isd, table, calmodel, calmachine))
end
function init_calibrate(method::AnalysisMethod{A}, i::Int; 
                            model = CalibrationModel{Linear}, 
                            kwargs...
                            ) where A
    analyte = analyteobj(method.conctable)[i]
    isd = isdof(analyteobj(method.conctable)[i], method)
    isd == analyte && return InternalCalibrator(analyte, first(getanalyte(method.conctable, analyte)))
    ord = sortperm(method.pointlevel)
    level = method.pointlevel[ord]
    conc = getanalyte(method.conctable, i)
    ya = getanalyte(method.signaltable, i)
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    id = sampleobj(method.signaltable)
    # Use Float64, bug of GLM
    numbertype = Float64 # eltype(conctable)
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        convert(numbertype, conc[findsample(method.conctable, Symbol(l))])
                    end, 
                    y = convert(Vector{numbertype}, y), 
                    x̂ = zeros(numbertype, length(id)), 
                    accuracy = zeros(numbertype, length(id)),
                    include = level .> 0
                    )
    calmodel = mkcalmodel(model; kwargs...)
    calmachine = mkcalmachine(calmodel, table)
    validate_quantify_calibrator!(ExternalCalibrator(analyte, isd, table, calmodel, calmachine))
end

function mkcalmodel(model::Type{CalibrationModel{T}}; wnm = "1", wfn = WFN[wnm]) where T
    f = getformula(T)
    model(f, wnm, wfn)
end
function mkcalmachine(model::CalibrationModel{T}, tbl) where T
    lm1 = lm(model.formula, tbl[tbl.include]; wts = getwts(model.wfn, tbl.x[tbl.include], tbl.y[tbl.include]))
    if T == Quadratic && lm1.model.pp.beta0[1] == 0
        m = hcat(ones(eltype(tbl.x), count(tbl.include)), tbl.x[tbl.include], tbl.x[tbl.include] .^ 2)
        sqrtw = diagm(sqrt.(getwts(model.wfn, tbl.x[tbl.include], tbl.y[tbl.include])))
        y = tbl.y[tbl.include]
        lm1.model.pp.beta0 = (sqrtw * m) \ (sqrtw * y)
        GLM.updateμ!(lm1.model.rr, predict(lm1, tbl[tbl.include]))
    end
    lm1
end
mkcalmachine(cal::ExternalCalibrator) = mkcalmachine(cal.model, cal.table)

# """
#     calibrate(tbl, formula, type, zero, weight)
#     calibrate(cal::ExternalCalibrator)

# Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

# Field `model` will be mutated for mutating version. Calling `calibrate!` on a `Batch` will apply `calibrate!` to all calibration curves.

# It returns `GLM` object, and input object for mutating version.
# """
# function calibrate(tbl, formula, type, zero, weight)
#     model = lm(formula, tbl[tbl.include]; wts = tbl.x[tbl.include] .^ weight)
#     if !type && !zero && model.model.pp.beta0[1] == 0
#         m = hcat(ones(eltype(tbl.x), count(tbl.include)), tbl.x[tbl.include], tbl.x[tbl.include] .^ 2)
#         sqrtw = diagm(sqrt.(tbl.x[tbl.include] .^ weight))
#         y = tbl.y[tbl.include]
#         model.model.pp.beta0 = (sqrtw * m) \ (sqrtw * y)
#         GLM.updateμ!(model.model.rr, predict(model, tbl[tbl.include]))
#     end
#     model
# end
# calibrate(cal::ExternalCalibrator) = calibrate(cal.table, cal.formula, cal.type, cal.zero, cal.weight)

function model_params(n, params; kwargs...)
    params = Table(params)
    m = length(params)
    n == m || throw(ArgumentError("Invalid `params`"))
    pn = propertynames(params)
    add = []
    for (k, v) in kwargs
        if !in(k, pn)
            push!(add, k => repeat([v], n))
        end 
    end
    Table(params; add...)
end
"""
    init_calibrate!(batch::Batch)

Initiate calibration for a batch with empty `batch.calibrator`.
"""
function init_calibrate!(batch::Batch; params = nothing, kwargs...)
    empty!(batch.calibrator)
    analyte = batch.std 
    if isnothing(params)
        for a in analyte
            push!(batch.calibrator, init_calibrate(batch, a; kwargs...))
        end
    else
        params = model_params(length(analyte), params; kwargs...)
        for (a, param) in zip(analyte, params)
            push!(batch.calibrator, init_calibrate(batch, a; param...))
        end
    end
    batch
end

"""
set_calibration!
"""
calibrate!(batch::Batch; kwargs...) = isempty(batch.calibrator) ? init_calibrate!(batch; kwargs...) : (foreach(calibrate!, batch.calibrator); batch)

"""
    calibrate!(batch::Batch)
    calibrate!(cal::Vector{<: AbstractCalibrator})
    calibrate!(cal::ExternalCalibrator)

Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

Field `model` will be mutated for mutating version. Calling `calibrate!` on a `Batch` will apply `calibrate!` to all calibration curves.

It returns `GLM` object input for mutating version.
"""
# calibrate!(cal::Vector{<: AbstractCalibrator}) = (foreach(calibrate!, cal); cal)
calibrate!(cal::InternalCalibrator) = cal
function calibrate!(cal::ExternalCalibrator)
    cal.machine = mkcalmachine(cal.model, cal.table)
    validate_quantify_calibrator!(cal)
end

# set_calibration(batch::Batch; model = CalibrationModel{Linear}, kwargs...) = isempty(batch.calibrator) ? init_calibration!(batch; model, kwargs...) : (foreach(calibrate!, batch.calibrator); batch)
# calibrate(cal::Vector{<: AbstractCalibrator}) = (foreach(calibrate!, cal); cal)
# calibrate(cal::InternalCalibrator) = cal
# function calibrate(cal::ExternalCalibrator)
#     cal.calibrationmachine = mkcalmachine(cal)
#     cal
# end

"""
    recalibrate!(batch::Batch{A}, analyte::B) where {A, B <: A}
    recalibrate!(batch::Batch, cal_id::Int)
    recalibrate!(cal::ExternalCalibrator, method::AnalysisMethod)

Recalibrate `analyte`, `batch.calibrator[cal_id]` or `cal` after modifying any parameters in method or calibrator.
"""
function recalibrate!(batch::Batch; params = nothing, kwargs...)
    if isnothing(params)
        for cal in batch.calibrator
            recalibrate!(cal; kwargs...)
        end
    else
        params = model_params(length(batch.calibrator), params; kwargs...)
        for (cal, param) in zip(batch.calibrator, params)
            recalibrate!(cal; param...)
        end
    end
    batch
end
function recalibrate!(batch::Batch{A}, analyte::B; model = nothing, kwargs...) where {A, B <: A}
    analytes = batch.analyte
    aid = findfirst(==(analyte), analytes)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in this batch"))
    ca = analytes[aid]
    cid = findfirst(x -> ==(x.analyte, ca), batch.calibrator)
    isnothing(cid) && throw(ArgumentError("No fitted calibration data for $analyte"))
    recalibrate!(batch, cid; model, kwargs...)
end
recalibrate!(batch::Batch, cal_id::Int; model = nothing, kwargs...) = recalibrate!(batch.calibrator[cal_id]; model, kwargs...)
# change in only model 
function recalibrate!(cal::ExternalCalibrator; model = nothing, kwargs...)
    if isnothing(model) && isempty(kwargs)
        calibrate!(cal) 
    elseif isnothing(model)
        cal.model = mkcalmodel(typeof(cal.model); kwargs...)
        cal.machine = mkcalmachine(cal.model, cal.table)
    else
        cal.model = mkcalmodel(model; kwargs...)
        cal.machine = mkcalmachine(cal.model, cal.table)
    end
    @warn "Call `quantify!` function on data to get the updated results."
    validate_quantify_calibrator!(cal)
end
function recalibrate!(cal::InternalCalibrator; model = nothing, kwargs...)
    cal
end

"""
"""
function retrieve_method!(batch::Batch; params = nothing, kwargs...)
    if isnothing(params)
        for cal in batch.calibrator
            retrieve_method!(cal, batch.method; kwargs...)
        end
    else
        params = model_params(length(batch.calibrator), params; kwargs...)
        for (cal, param) in zip(batch.calibrator, params)
            retrieve_method!(cal, batch.method; param...)
        end
    end
    batch
end
function retrieve_method!(batch::Batch{A}, analyte::B; model = nothing, kwargs...) where {A, B <: A}
    analytes = batch.analyte
    aid = findfirst(==(analyte), analytes)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in this batch"))
    ca = analytes[aid]
    cid = findfirst(x -> ==(x.analyte, ca), batch.calibrator)
    isnothing(cid) && throw(ArgumentError("No fitted calibration data for $analyte"))
    retrieve_method!(batch, cid; model, kwargs...)
end
retrieve_method!(batch::Batch, cal_id::Int; model = nothing, kwargs...) = retrieve_method!(batch.calibrator[cal_id], batch.method; model, kwargs...)
# change in method and/or model
function retrieve_method!(cal::ExternalCalibrator, method::AnalysisMethod; model = nothing, kwargs...)
    isd = isdof(cal.analyte, method)
    cal.isd = isd
    ord = sortperm(method.pointlevel)
    ya = getanalyte(method.signaltable, cal.analyte)
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    cal.table.y .= y
    # cal.table.include .= true
    recalibrate!(cal; model, kwargs...)
end
function retrieve_method!(cal::InternalCalibrator, method::AnalysisMethod; model = nothing, kwargs...)
    cal.conc = first(getanalyte(method.conctable, cal.analyte))
    cal
end
"""
    set_isd!(batch::Batch{A}, analyte::B, isd::C) where {A, B <: A, C <: A}

Set internal standard of `analyte` to `isd`. `analyte` must have a calibration curve and all analytes use this curve will be affected.
"""
function set_isd!(batch::Batch{A}, analyte::B, isd::Union{C, Nothing} = nothing; model = nothing, kwargs...) where {A, B <: A, C <: A}
    cid = findfirst(x -> ==(x.analyte, analyte), batch.calibrator)
    isnothing(cid) && throw(ArgumentError("Analyte $analyte is not a calibrator"))
    set_isd!(batch.method, analyte, isd)
    retrieve_method!(batch, cid; model, kwargs...)
    batch
end
function set_isd!(method::AnalysisMethod{A}, analyte::B, isd::Union{C, Nothing} = nothing) where {A, B <: A, C <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = isnothing(isd) ? 0 : findfirst(==(isd), method.analytetable.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    cid = findall(==(aid), method.analytetable.calibrator)
    method.analytetable.isd[cid] .= iid
    method
end

function set_isd!(method::AnalysisMethod, aid::Int, iid::Int)
    cid = findall(==(aid), method.analytetable.std)
    method.analytetable.isd[cid] .= iid
    method
end

"""
    set_std!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); type = true, zero = false, weight = 0) where {A, B <: A, C <: A, D <: A}

Set calibration standard for `analyte` to `cal`. `cal` must have a calibration curve. The calibration curve use the new internal standrard `isd`.
"""
function set_std!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); model = nothing, kwargs...) where {A, B <: A, C <: A, D <: A}
    aid = findfirst(==(analyte), batch.method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    cid = findfirst(==(cal), batch.method.analytetable.analyte)
    isnothing(cid) && throw(ArgumentError("Analyte $cal is not in the method"))
    iid = isnothing(isd) ? 0 : findfirst(==(isd), batch.method.analytetable.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    batch.method.analytetable.std[aid] = cid
    calid = findfirst(x -> ==(x.analyte, cal), batch.calibrator)
    if isnothing(calid)
        set_isd!(batch.method, cid, iid)
        push!(batch.calibrator, init_calibrate(batch.method, cal; model, kwargs...))
        @warn "Call `quantify!` function on data to get the updated results."

    else
        batch.calibrator[calid].isd == isd || (set_isd!(batch.method, cid, iid); retrieve_method!(batch, calid; model, kwargs...))
    end
    batch
end

"""
    replace_std!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); type = true, zero = false, weight = 0) where {A, B <: A, C <: A, D <: A}

Delete calibration curve of `analyte` and replace it with calibration curve of `cal`. `analyte` must have a calibration curve. If `cal` has no caliration curve, a new calibration curve is fitted. The calibration curve use the new internal standard `isd`. All analytes use the calibration standard `analyte` are affected.
"""
function replace_std!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); model = nothing, kwargs...) where {A, B <: A, C <: A, D <: A}
    calid1 = findfirst(x -> ==(first(x.analyte), cal), batch.calibrator)
    isnothing(calid1) && throw(ArgumentError("No fitted calibration data for $analyte"))
    aid = findfirst(==(analyte), batch.method.analytetable.analyte)
    isnothing(cid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    cid = findfirst(==(cal), batch.method.analytetable.analyte)
    isnothing(cid) && throw(ArgumentError("Analyte $cal is not in the method"))
    iid = isnothing(isd) ? 0 : findfirst(==(isd), batch.method.analytetable.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    tocal = findall(==(aid), batch.method.analytetable.std)
    batch.method.analytetable.std[tocal] .= cid
    calid = findfirst(x -> ==(first(x.analyte), cal), batch.calibrator)
    if isnothing(calid)
        set_isd!(batch.method, cid, iid)
        push!(batch.calibrator, init_calibrate(batch.method, cal; model, kwargs...))
        method.analytetable.std[aid] = cid
        @warn "Call `quantify!` function on data to get the updated results."
    else
        batch.calibrator[calid].isd == isd || (set_isd!(batch.method, cid, iid); retrieve_method!(batch, calid; model, kwargs...))
    end
    delete!(batch.calibrator, calid1)
    batch
end
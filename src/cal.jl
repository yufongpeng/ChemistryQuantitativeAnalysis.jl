"""
    validate_inv_predict!(batch::Batch)
    validate_inv_predict!(cal::Vector{<: AbstractCalibration})
    validate_inv_predict!(cal::SingleCalibration)
    validate_inv_predict!(cal::MultipleCalibration)

Inverse predict concentration, update each `cal.table.x̂` with the result(s) and returns `cal` or `batch`.
"""
validate_inv_predict!(batch::Batch) = (validate_inv_predict!(batch.calibration); batch)
validate_inv_predict!(cal::Vector{<: AbstractCalibration}) = (foreach(validate_inv_predict!, cal); cal)
validate_inv_predict!(cal::SingleCalibration) = cal
function validate_inv_predict!(cal::MultipleCalibration)
    cal.table.x̂ .= inv_predict(cal, cal.table.y)
    cal
end

function fill_result!(dt::ColumnDataTable, result)
    for (a, c) in zip(eachanalyte(dt), result)
        a .= c
    end
    dt
end
function fill_result!(dt::RowDataTable, result)
    for (i, p) in enumerate(eachsample(dt))
        p .= getindex.(result, i)
    end
    dt
end 

"""
    validate_accuracy!(batch::Batch)
    validate_accuracy!(cal::Vector{<: AbstractCalibration})
    validate_accuracy!(cal::MultipleCalibration)
    validate_accuracy!(cal::SingleCalibration)

Calculate accuracy for analyte specified by `cal`, update each `cal.table.x̂` with the result(s), and return `cal` or `batch`.
"""
validate_accuracy!(batch::Batch) = (validate_accuracy!(batch.calibration); batch)
validate_accuracy!(cal::Vector{<: AbstractCalibration}) = (foreach(validate_accuracy!, cal); cal)
function validate_accuracy!(cal::MultipleCalibration)
    cal.table.accuracy .= accuracy(cal.table.x̂, cal.table.x)
    cal
end
validate_accuracy!(cal::SingleCalibration) = cal

"""
    const inv_predict_accuracy! = accuracy! ∘ inv_predict!

Apply `inv_predict!` and `accuracy!` subsequantly to `Batch` or `AbstractCalibration`.
"""
const validate! = validate_accuracy! ∘ validate_inv_predict!

"""
    calibration(anisd::Tuple, tbl::Table; analyte = 1, isd = 0, type = true, zero = false, weight = 0)
    calibration(batch::Batch{A}, ana::B; id = sample(batch.method.signaltable), isd = isdof(batch.method, analyte),
                        type = true, zero = false, weight = 0) where {A, B <: A}
    calibration(method::AnalysisMethod{A}, ana::B; id = sample(method.signaltable), isd = isdof(method, analyte),
                type = true, zero = false, weight = 0) where {A, B <: A}

Create `MultipleCalibration`.

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
function calibration(anisd::Tuple, tbl::Table;
                    type = true, 
                    zero = false, 
                    weight = 0
                    )
    id = findall(x -> isa(x, Number), tbl.y)
    tbl = tbl[id]
    xlevel = unique(tbl.x)
    table = :id in propertynames(tbl) ? tbl : Table((; id = collect(1:length(tbl)), ), tbl)
    table = :level in propertynames(tbl) ? table : Table(table; level = [findfirst(x -> i == x, xlevel) for i in table.x])
    table = :include in propertynames(tbl) ? table : Table(table; include = trues(length(table)))
    f = getformula(type, zero)
    model = calibrate(table, f, type, zero, weight)
    table = Table(; id = table.id, level = table.level, y = table.y, x = table.x, x̂ = zeros(eltype(table.x), length(table)), accuracy = zeros(eltype(table.x), length(table)), include = table.include)
    validate!(MultipleCalibration(anisd, type, zero, weight, f, table, model))
end
calibration(batch::Batch{A}, analyte::B; 
                    id = sampleobj(batch.method.signaltable),
                    isd = nothing,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B} = calibration(batch.method, analyte; id, isd, type, zero, weight)
function calibration(method::AnalysisMethod{A}, analyte::B; 
                    id = sampleobj(method.signaltable),
                    isd = missing,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B <: A}

    isd = (ismissing(isd) ? isdof(analyte, method) : isd)
    ord = sortperm(method.pointlevel)
    level = method.pointlevel[ord]
    conc = getanalyte(method.conctable, analyte)
    ya = getanalyte(method.signaltable, analyte)
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    # Use Float64, bug of GLM
    numbertype = Float64 # eltype(conctable)
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        convert(Float64, conc[findsample(method.conctable, Symbol(l))])
                    end, 
                    y = convert(Vector{numbertype}, y), 
                    x̂ = zeros(numbertype, length(id)), 
                    accuracy = zeros(numbertype, length(id)),
                    include = trues(length(id))
                    )
    f = getformula(type, zero)
    weight = numbertype(weight)
    model = calibrate(table, f, type, zero, weight)
    validate!(MultipleCalibration((analyte, isd), type, zero, weight, f, table, model))
end
function calibration(method::AnalysisMethod{A}, i::Int; 
                    id = sampleobj(method.signaltable),
                    isd = missing,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where A

    isd = ismissing(isd) ? isdof(analyteobj(method.conctable)[i], method) : isd
    ord = sortperm(method.pointlevel)
    level = method.pointlevel[ord]
    conc = getanalyte(method.conctable, i)
    ya = getanalyte(method.signaltable, i)
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    # Use Float64, bug of GLM
    numbertype = Float64 # eltype(conctable)
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        convert(Float64, conc[findsample(method.conctable, Symbol(l))])
                    end, 
                    y = convert(Vector{numbertype}, y), 
                    x̂ = zeros(numbertype, length(id)), 
                    accuracy = zeros(numbertype, length(id)),
                    include = trues(length(id))
                    )
    f = getformula(type, zero)
    weight = numbertype(weight)
    model = calibrate(table, f, type, zero, weight)
    validate!(MultipleCalibration((analyteobj(method.conctable)[i], isd), type, zero, Float64(weight), f, table, model))
end

"""
    calibrate(tbl, formula, type, zero, weight)
    calibrate(cal::MultipleCalibration)

Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

Field `model` will be mutated for mutating version. Calling `calibrate!` on a `Batch` will apply `calibrate!` to all calibration curves.

It returns `GLM` object, and input object for mutating version.
"""
function calibrate(tbl, formula, type, zero, weight)
    model = lm(formula, tbl[tbl.include]; wts = tbl.x[tbl.include] .^ weight)
    if !type && !zero && model.model.pp.beta0[1] == 0
        m = hcat(ones(eltype(tbl.x), count(tbl.include)), tbl.x[tbl.include], tbl.x[tbl.include] .^ 2)
        sqrtw = diagm(sqrt.(tbl.x[tbl.include] .^ weight))
        y = tbl.y[tbl.include]
        model.model.pp.beta0 = (sqrtw * m) \ (sqrtw * y)
        GLM.updateμ!(model.model.rr, predict(model, tbl[tbl.include]))
    end
    model
end
calibrate(cal::MultipleCalibration) = calibrate(cal.table, cal.formula, cal.type, cal.zero, cal.weight)

"""
    calibrate!(batch::Batch)
    calibrate!(cal::Vector{<: AbstractCalibration})
    calibrate!(cal::MultipleCalibration)

Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

Field `model` will be mutated for mutating version. Calling `calibrate!` on a `Batch` will apply `calibrate!` to all calibration curves.

It returns `GLM` object input for mutating version.
"""
calibrate!(batch::Batch) = (foreach(calibrate!, batch.calibration); batch)
calibrate!(cal::Vector{<: AbstractCalibration}) = (foreach(calibrate!, cal); cal)
function calibrate!(cal::MultipleCalibration)
    cal.model = calibrate(cal.table, cal.formula, cal.type, cal.zero, cal.weight)
    cal
end

"""
    update_calibration!(batch::Batch{A}, analyte::B) where {A, B <: A}
    update_calibration!(batch::Batch, cal_id::Int)
    update_calibration!(cal::MultipleCalibration, method::AnalysisMethod)

Update calibration obeject for `analyte`, `batch.calibration[cal_id]` or `cal` after modifying any parameters in method or calibration object.
"""
function update_calibration!(batch::Batch{A}, analyte::B) where {A, B <: A}
    analytes = batch.analyte
    aid = findfirst(==(analyte), analytes)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in this batch"))
    ca = analytes[aid]
    cid = findfirst(x -> ==(first(x.analyte), ca), batch.calibration)
    isnothing(cid) && throw(ArgumentError("No fitted calibration data for $analyte"))
    update_calibration!(batch, cid)
end
update_calibration!(batch::Batch, cal_id::Int) = update_calibration!(batch.calibration[cal_id], batch.method)
function update_calibration!(cal::MultipleCalibration, method::AnalysisMethod)
    isd = isdof(first(cal.analyte), method)
    cal.analyte = (first(cal.analyte), isd)
    ord = sortperm(method.pointlevel)
    ya = getanalyte(method.signaltable, first(cal.analyte))
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    cal.table.y .= y
    # cal.table.include .= true
    cal.formula = getformula(cal)
    calibrate!(cal)
    validate!(cal)
end
function update_calibration!(cal::SingleCalibration, method::AnalysisMethod)
    cal.conc = first(getanalyte(method.conctable, first(cal.analyte)))
    cal
end

"""
    set_isd!(batch::Batch{A}, analyte::B, isd::C) where {A, B <: A, C <: A}

Set internal standard of `analyte` to `isd`. `analyte` must have a calibration curve and all analytes use this curve will be affected.
"""
function set_isd!(batch::Batch{A}, analyte::B, isd::Union{C, Nothing} = nothing) where {A, B <: A, C <: A}
    cid = findfirst(x -> ==(first(x.analyte), analyte), batch.calibration)
    isnothing(cid) && throw(ArgumentError("No fitted calibration data for $analyte"))
    set_isd!(batch.method, analyte, isd)
    update_calibration!(batch, cid)
    @warn "Call quantification function on data to get the updated results."
    batch
end
function set_isd!(method::AnalysisMethod{A}, analyte::B, isd::Union{C, Nothing} = nothing) where {A, B <: A, C <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = isnothing(isd) ? 0 : findfirst(==(isd), method.analytetable.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    cid = findall(==(aid), method.analytetable.calibration)
    method.analytetable.isd[cid] .= iid
    method
end

function set_isd!(method::AnalysisMethod, aid::Int, iid::Int)
    cid = findall(==(aid), method.analytetable.calibration)
    method.analytetable.isd[cid] .= iid
    method
end

"""
    set_cal!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); type = true, zero = false, weight = 0) where {A, B <: A, C <: A, D <: A}

Set calibration curve for `analyte` to calibration curve of `cal`. `cal` must have a calibration curve. The calibration curve use the new `isd`.
"""
function set_cal!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); type = true, zero = false, weight = 0) where {A, B <: A, C <: A, D <: A}
    aid = findfirst(==(analyte), batch.method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    cid = findfirst(==(cal), batch.method.analytetable.analyte)
    isnothing(cid) && throw(ArgumentError("Analyte $cal is not in the method"))
    iid = isnothing(isd) ? 0 : findfirst(==(isd), batch.method.analytetable.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    batch.method.analytetable.calibration[aid] = cid
    calid = findfirst(x -> ==(first(x.analyte), cal), batch.calibration)
    if isnothing(calid)
        if length(sampleobj(method.conctable)) > 1
            push!(batch.calibration, calibration(batch.method, cal; isd, type, zero, weight))
        else
            push!(batch.calibration, SingleCalibration((cal, ), first(getanalyte(method.conctable, cal))))
        end
    else
        last(batch.calibration[calid].analyte) == isd || (set_isd!(batch.method, cid, iid); update_calibration!(batch, calid))
    end
    @warn "Call quantification function on data to get the updated results."
    batch
end

"""
    replace_cal!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); type = true, zero = false, weight = 0) where {A, B <: A, C <: A, D <: A}

Delete calibration curve of `analyte` and replace it with calibration curve of `cal`. `analyte` must have a calibration curve. If `cal` has no caliration curve, a new calibration curve is fitted. The calibration curve use the new `isd`. All analytes use the calibration curve of `analyte` are affected.
"""
function replace_cal!(batch::Batch{A}, analyte::B, cal::C, isd::Union{D, Nothing} = isdof(cal, batch.method); type = true, zero = false, weight = 0) where {A, B <: A, C <: A, D <: A}
    calid1 = findfirst(x -> ==(first(x.analyte), cal), batch.calibration)
    isnothing(calid1) && throw(ArgumentError("No fitted calibration data for $analyte"))
    aid = findfirst(==(analyte), batch.method.analytetable.analyte)
    isnothing(cid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    cid = findfirst(==(cal), batch.method.analytetable.analyte)
    isnothing(cid) && throw(ArgumentError("Analyte $cal is not in the method"))
    iid = isnothing(isd) ? 0 : findfirst(==(isd), batch.method.analytetable.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    tocal = findall(==(aid), batch.method.analytetable.calibration)
    batch.method.analytetable.calibration[tocal] .= cid
    calid = findfirst(x -> ==(first(x.analyte), cal), batch.calibration)
    if isnothing(calid)
        if length(sampleobj(method.conctable)) > 1
            push!(batch.calibration, calibration(batch.method, cal; isd, type, zero, weight))
        else
            push!(batch.calibration, SingleCalibration((cal, ), first(getanalyte(method.conctable, cal))))
        end
        method.analytetable.calibration[aid] = cid

    else
        last(batch.calibration[calid].analyte) == isd || (set_isd!(batch.method, cid, iid); update_calibration!(batch, calid))
    end
    delete!(batch.calibration, calid1)
    @warn "Call quantification function on data to get the updated results."
    batch
end
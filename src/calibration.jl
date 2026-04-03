"""
    edit_method!(batch::Batch; kwargs...)
    edit_method!(batch::Batch, params; kwargs...)
    edit_method!(batch::Batch, analyte_std, isd, std_isd, params = nothing; kwargs...) -> Batch
    edit_method!(method::AnalysisMethod; kwargs...)
    edit_method!(method::AnalysisMethod, params; kwargs...)
    edit_method!(method::AnalysisMethod, analyte_std, isd, std_isd, params = nothing; kwargs...) -> AnalysisMethod

Edit analysis method.

* `params`: object compatible with `Table.jl` interfaces. Columns may include
    * `:analyte`: analytes to be assined parameters. Both analyte objects or ids are accepted. If this is not provided, the first n analytes are used, where n is the length of `params`.
    * `:isd`: assigin internal standards for calibration standards (do not assign if it is not a calibration standard). Both analyte objects or ids are accepted. Use `0` to not assign isd and `-1` if it is an internal standard.
    * `:std`: assigin calibration standards. Both analyte objects or ids are accepted. Use `0` to not assign std and `-1` if it is an internal standard.
    * `:model`: calibration model types (`AbstractCalibrationModel`). Use `Nothing` to not assign meodel.
    * Keyword arguments of function `mkcalmodel`.

Calibration standards and internal standards can be provided separately
* `analyte_std::Vector{<: Pair}`: analyte-std pairs. Assign calibration standard to analyte. Both analyte objects or ids are accepted. Use `0` to not assign std.
* `isd::Vector`: analyte objects or ids of internal standards.  
* `std_isd::Vector{<: Pair}`: std_isd pairs. Assign internal stadards to calibration standards. Both analyte objects or ids are accepted. Use `0` to not assign isd.

Other keyword arguments are applied to all included models in function `mkcalmodel`.
"""
function edit_method!(batch::Batch, args...; kwargs...)
    edit_method!(batch.method, args...; kwargs...) 
    batch
end
function edit_method!(method::AnalysisMethod; kwargs...)
    _edit_method!(method; kwargs...)
end
function edit_method!(method::AnalysisMethod, params; kwargs...)
    analyte_std, isd, std_isd, params = cal_params(params)
    edit_method!(method, analyte_std, isd, std_isd, params; kwargs...) 
end
function edit_method!(method::AnalysisMethod, analyte_std::Vector, isd::Vector, std_isd::Vector, params = nothing; kwargs...)
    if !isempty(analyte_std) || !isempty(isd) || !isempty(std_isd)
        _edit_method!(
            method, 
            [getanalyteid_check(method, i) => getanalyteid_check(method, j) for (i, j) in analyte_std],
            [getanalyteid_check(method, i) for i in isd],
            [getanalyteid_check(method, i) => getanalyteid_check(method, j) for (i, j) in std_isd],
        )
    end
    _edit_method!(method, params; kwargs...)
end

function _edit_method!(method::AnalysisMethod, analyte_std::Vector{<: Pair}, isd::Vector, std_isd::Vector{<: Pair})
    visd = deepcopy(method.analytetable.isd)
    visd[isd] .= -1
    vstd = deepcopy(method.analytetable.std)
    vstd[isd] .= -1
    # input check
    for i in eachindex(method.analyte)
        if visd[i] < 0
            # assign as isd, check if it is used in external calibration
            # it is not assigned as std for other isd except itself (internal calibration)
            isnothing(findfirst(x -> first(x) == i, analyte_std)) || throw(ArgumentError("Analyte $(method.analyte[i]) ($i) is interanl standard and cannot be external calibrated."))
            # throw(ArgumentError("Analyte $(method.analyte[i]) cannot be both calibration standard and internal standard."))
            j = findfirst(x -> first(x) == i, std_isd)
            isnothing(j) || ==(std_isd[j]...) || throw(ArgumentError("Analyte $(method.analyte[i]) ($i) cannot be isd and futher internal calibrated."))
        else
            # find std if it is external calibrated 
            # default std otherwise
            j = findfirst(x -> first(x) == i, analyte_std)
            if isnothing(j)
                j = vstd[i]
            else
                j = last(analyte_std[j])
            end
            # find isd for that std
            k = findfirst(x -> first(x) == i, std_isd)
            # if calibrated by itself (it is a std) 
            # elseif it is not assigned isd 
            i == j || isnothing(k) || throw(ArgumentError("Analyte $(method.analyte[i]) ($i) is not a calibration standard, so it can be internal calibrated but cannot be assgined internal standard (assign isd to calibration standard only)."))
            vstd[i] = j 
        end
    end
    for i in eachindex(method.analyte)
        if visd[i] < 0
            vstd[i] = -1 
        elseif vstd[i] < 0
            visd[i] = -1 
        end
        j = findfirst(x -> first(x) == vstd[i], std_isd)
        if !isnothing(j)
            visd[i] = last(std_isd[j])
        end
    end
    # final check
    analytetable_check(method, vstd, visd)
    method.analytetable.std .= vstd
    method.analytetable.isd .= visd
    method
end

"""
    analytetable_check(batch::Batch)
    analytetable_check(batch::Batch, vstd, visd)    
    analytetable_check(method::AnalysisMethod)
    analytetable_check(method::AnalysisMethod, vstd, visd)

Check if analytetable is valid.
"""
analytetable_check(batch::Batch, args...) = analytetable_check(batch.method, args...)
analytetable_check(method::AnalysisMethod, vstd, visd) = analytetable_check(method.analyte, vstd, visd)
analytetable_check(method::AnalysisMethod) = analytetable_check(method.analyte, method.analytetable.std, method.analytetable.isd)
function analytetable_check(analyte, vstd, visd)
    for i in eachindex(analyte)
        if visd[i] < 0
            # it is not assigned as std for other isd except itself (internal calibration)
            j = findall(==(i), vstd)
            isempty(j) || all(k -> visd[k] == i, j) || throw(ArgumentError("Analyte $(analyte[i]) ($i) cannot be both calibration standard and internal standard."))
        elseif vstd[i] == i 
            # external calibrated, check isd is isd or no isd
            j = visd[i]
            j == 0 && continue
            j == i && throw(ArgumentError("Analyte $(analyte[i]) ($i) cannot be internal calibrated by itself."))
            j == 0 || visd[j] < 0 || throw(ArgumentError("Analyte $(analyte[i]) ($i) is internal calibrated by non internal standard."))
        else
            j = visd[i]
            j == i && throw(ArgumentError("Analyte $(analyte[i]) ($i) cannot be internal calibrated by itself."))
            k = vstd[i]
            k == 0 && continue
            if j == k && visd[j] < 0 
                continue 
            elseif j == k 
                throw(ArgumentError("Analyte $(analyte[i]) ($i) is internal calibrated by non internal standard."))
            elseif vstd[k] != k 
                throw(ArgumentError("Analyte $(analyte[i]) ($i) is external calibrated by non calibration standard."))
            elseif visd[k] != j
                throw(ArgumentError("Analyte $(analyte[i]) ($i) is internal calibrated by incorrect internal standard."))
            end 
        end
    end
end

_edit_method!(method::AnalysisMethod; kwargs...) = _edit_method!(method, nothing; kwargs...)
function _edit_method!(method::AnalysisMethod, params::Nothing; kwargs...)
    isempty(kwargs) && return method
    model, kwargs = model_params(params; kwargs...)
    method.analytetable.model .= [modifycalmodel(m, model; kwargs...) for m in method.analytetable.model]
    method
end
function _edit_method!(method::AnalysisMethod, params; kwargs...)
    analyte, models, param = model_params(params; kwargs...)
    aid = getanalyteid_check.(Ref(method), analyte)
    for (a, model, kwargs) in zip(aid, models, param)
        method.analytetable.model[a] = modifycalmodel(method.analytetable.model[a], model; kwargs...)
    end
    method
end

"""
    assign_std!(batch::Batch, analyte, cal, isd = isdof(cal, batch.method); calibrate = false, model = Nothing, kwargs...) -> Batch

Assign calibration standard for `analyte` to `cal`. 
    
If `analyte` is a calibration standard, all other analytes calibrated by `analyte` are affected.

If `cal` is not a calibration standard, a new calibration curve is fitted. 

The calibration curve use the new internal standard `isd` (`nothing` means no isd). All analytes use the calibration standard `analyte` are affected.
If `cal` and `isd` are the same, `InternalCalibrator` is used.

# Keyword Arguments
* `calibrate`: apply `calibrate!` subsequently or not. This function do not affect other unrelated analytes.
* `model`: calibration model type.
* Other keyword arguments for function `mkcalmodel`.
"""
function assign_std!(batch::Batch{A}, 
                analyte::B, 
                cal::C, 
                isd::Union{D, Nothing} = isdof(cal, batch.method); 
                calibrate = false, 
                model = Nothing, 
                kwargs...) where {A, B <: A, C <: A, D <: A}
    aid = getanalyteid_check(batch, analyte)
    cid = getanalyteid_check(batch, cal)
    iid = getisdid_check(batch, isd)
    assign_std!(batch, aid, cid, iid; calibrate, model, kwargs...)
end
function assign_std!(batch::Batch, 
                aid::Int, 
                cid::Int, 
                iid::Int = batch.method.analytetable.isd[cid]; 
                calibrate = false, 
                model = Nothing, 
                kwargs...) 
    assign_std!(batch.method, aid, cid, iid; model, kwargs...)
    calibrate || return batch 
    if aid == cid 
        # No deletion
        calibrate!(batch, cid) 
    else
        calid = findfirst(x -> ==(x.analyte, batch.analyte[aid]), batch.calibrator)
        if !isnothing(calid) 
            # delete replaced cal
            deleteat!(batch.calibrator, calid)
        end
        calibrate!(batch, cid)
    end
end
function assign_std!(method::AnalysisMethod{A}, 
                analyte::B, 
                cal::C, 
                isd::Union{D, Nothing} = isdof(cal, method); 
                model = Nothing, 
                kwargs...) where {A, B <: A, C <: A, D <: A}
    aid = getanalyteid_check(batch, analyte)
    cid = getanalyteid_check(batch, cal)
    iid = getisdid_check(batch, isd)
    assign_std!(method, aid, cid, iid; model, kwargs...)
end
function assign_std!(method::AnalysisMethod, 
                aid::Int, 
                cid::Int, 
                iid::Int = method.analytetable.isd[cid]; 
                model = Nothing, 
                kwargs...)
    method.analytetable.isd[aid] < 0 && throw(ArgumentError("Analyte $(method.analytetable.analyte[aid]) ($aid) cannot be an internal standard."))
    tocal = findall(==(aid), method.analytetable.std)
    push!(tocal, aid)
    if method.analytetable.isd[cid] < 0
        # check if cid is really others isd
        j = findall(==(cid), method.analytetable.isd)
        isempty(j) || throw(ArgumentError("Analyte $(method.analytetable.analyte[cid]) ($cid) cannot be an internal standard."))
        push!(tocal, cid)
    else
        push!(tocal, cid)
    end
    unique!(tocal)
    if iid <= 0 
        # no isd
        method.analytetable.isd[tocal] .= 0
    elseif method.analytetable.isd[iid] < 0
        method.analytetable.isd[tocal] .= iid 
    else 
        # check if iid is others std except interanl calibrated
        j = findall(==(iid), method.analytetable.std)
        isempty(j) || all(k -> method.analytetable.isd[k] == iid, j) || throw(ArgumentError("Analyte $(method.analytetable.analyte[iid]) ($iid) cannot be a calibration standard."))
        method.analytetable.isd[tocal] .= iid 
    end
    method.analytetable.std[tocal] .= cid
    method.analytetable.model[cid] = modifycalmodel(method.analytetable.model[cid], model; kwargs...)
    method
end

"""
    assign_isd!(batch::Batch, analyte, isd = nothing; calibrate = false, model = Nothing, kwargs...) -> Batch

Assign internal standard of `analyte` to `isd` (`nothing` means no isd). 

`analyte` must be a calibration standard and all analytes use this curve will be affected.

# Keyword Arguments
* `calibrate`: apply `calibrate!` subsequently or not. This function do not affect other unrelated analytes.
* `model`: calibration model type.
* Other keyword arguments for function `mkcalmodel`.
"""
function assign_isd!(batch::Batch{A}, 
                analyte::B, 
                isd::Union{C, Nothing} = nothing; 
                calibrate = false, 
                model = Nothing, 
                kwargs...) where {A, B <: A, C <: A}
    aid = getanalyteid_check(batch, analyte)
    iid = getisdid_check(batch, isd)
    assign_isd!(batch, aid, iid; calibrate, model, kwargs...)
end
function assign_isd!(batch::Batch, 
                aid::Int, 
                iid::Int = 0; 
                calibrate = false, 
                model = Nothing, 
                kwargs...) 
    assign_isd!(batch.method, aid, iid; model, kwargs...)
    calibrate ? calibrate!(batch, aid) : batch 
end
function assign_isd!(method::AnalysisMethod{A}, 
                analyte::B, 
                isd::Union{C, Nothing} = nothing; 
                model = Nothing, 
                kwargs...) where {A, B <: A, C <: A}
    aid = getanalyteid_check(batch, analyte)
    iid = getisdid_check(batch, isd)
    assign_isd!(method, aid, iid; model, kwargs...)
end
function assign_isd!(method::AnalysisMethod, 
                aid::Int, 
                iid::Int; 
                model = Nothing, 
                kwargs...)
    method.analytetable.std[aid] == aid || throw(ArgumentError("Analyte $(method.analytetable.analyte[aid]) ($aid) must be a calibration standard."))
    tocal = findall(==(aid), method.analytetable.std)
    if iid <= 0 
        method.analytetable.isd[tocal] .= 0
    elseif method.analytetable.isd[iid] < 0
        method.analytetable.isd[tocal] .= iid 
    else 
        j = findall(==(iid), method.analytetable.std)
        isempty(j) || all(k -> method.analytetable.isd[k] == iid, j) || throw(ArgumentError("Analyte $(method.analytetable.analyte[iid]) ($iid) cannot be a calibration standard."))
        method.analytetable.isd[tocal] .= iid 
    end
    method.analytetable.model[aid] = modifycalmodel(method.analytetable.model[aid], model; kwargs...)
    method
end

@deprecate replace_std! assign_std!

"""
    calibrate!(batch::Batch; check = true) -> Batch

Initiate calibration for a batch with empty `batch.calibrator` with `check` doing `analytetable_check` first.
"""
function calibrate!(batch::Batch; check = true)
    check && analytetable_check(batch)
    empty!(batch.calibrator)
    analyte = batch.std 
    for a in analyte
        push!(batch.calibrator, calibrate(batch, a))
    end
    batch
end

"""
    calibrate!(batch::Batch, analyte)
    calibrate!(batch::Batch, analyte_id::Int) 
  
Initiate calibration for `analyte` or `batch.analyte[analyte_id]`. 
"""
function calibrate!(batch::Batch{A}, analyte::B) where {A, B <: A}
    cal_id = findfirst(x -> ==(x.analyte, analyte), batch.calibrator)
    if isnothing(cal_id) 
        push!(batch.calibrator, calibrate(batch, analyte))
    else
        batch.calibrator[cal_id] = calibrate(batch, analyte)
    end
end
calibrate!(batch::Batch, cid::Int) = 
    calibrate!(batch, batch.analyte[cid])

# recalibration after point select
"""
    externalcalibrate(analyte, isd, tbl::Table, model; modeltype = Nothing, kwargs...) -> ExternalCalibrator

Construct `ExternalCalibrator` from existing data.

* `tbl::TypedTable.Table`: this will be stored as `table`. It is clean-up calibration data for points selection.
* `model`: calibration model.
* `modeltype`: new calibration model type.
* `kwargs`: new parameters; keyword arguments for function `mkcalmodel`.  
"""
function externalcalibrate(analyte::A, isd, tbl::Table, model; modeltype = Nothing, kwargs...) where A
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
    calmodel = modifycalmodel(model, modeltype; kwargs...)
    calmachine = mkcalmachine(calmodel, table)
    table = Table(; id = table.id, level = table.level, y = table.y, x = table.x, x̂ = zeros(eltype(table.x), length(table)), accuracy = zeros(eltype(table.x), length(table)), include = table.include)
    analyze_calibrator!(ExternalCalibrator(analyte, isd, table, calmodel, calmachine))
end

"""
    calibrate(batch::Batch, analyte)
    calibrate(batch::Batch, analyte_id::Int) -> Batch
    calibrate(method::AnalysisMethod{A}, analyte::B) where {A, B <: A}
    calibrate(method::AnalysisMethod{A}, analyte_id::Int) -> AnalysisMethod

Construct `InternalCalibrator` or `ExternalCalibrator` for `analyte` or `analyte` id.
"""
calibrate(batch::Batch, analyte) = calibrate(batch.method, analyte)
function calibrate(method::AnalysisMethod{A}, analyte::B) where {A, B <: A}
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
                    y = convert(Vector{numbertype}, y), 
                    x = map(level) do l
                        convert(numbertype, conc[findsample(method.conctable, Symbol(l))])
                    end, 
                    x̂ = zeros(numbertype, length(id)), 
                    accuracy = zeros(numbertype, length(id)),
                    include = collect(level .> 0)
                    )
    calmodel = method.analytetable.model[getanalyteid(method, analyte)]
    isnothing(calmodel) && throw(ArgumentError(""))
    calmachine = mkcalmachine(calmodel, table)
    analyze_calibrator!(ExternalCalibrator(analyte, isd, table, calmodel, calmachine))
end
calibrate(method::AnalysisMethod{A}, i::Int) where A = calibrate(method, method.analyte[i])

"""
    model_calibrator!(batch::Batch, params = nothing; kwargs...)
    model_calibrator!(batch::Batch, analyte_id::Int; kwargs...)
    model_calibrator!(batch::Batch{A}, analyte::B; kwargs...) where {A, B <: A} -> Batch
    model_calibrator!(batch::Batch, cal::AbstractCalibrator; kwargs...)

Modify model-related paramters and recalibrate `analyte`, `batch.analyte[analyte_id]`, or `cal.analyte`.

* `params`: object compatible with `Table.jl` interfaces. Columns may include
    * `:analyte`: analytes to be assined parameters. Both analyte objects or ids are accepted. If this is not provided, the first n analytes are used, where n is the length of `models`.
    * `:model`: calibration model types (`AbstractCalibrationModel`). Use `Nothing` to not assign meodel.
    * Keyword arguments of function `mkcalmodel`.
* Keyword arguments are applied to all included models in function `mkcalmodel`.
"""
function model_calibrator!(batch::Batch, params = nothing; kwargs...)
    _edit_method!(batch.method, params; kwargs...)
    for cal in batch.calibrator
        aid = getanalyteid_check(batch, cal.analyte)
        model_calibrator!(cal, batch.method.analytetable.model[aid])
    end
    batch
end
function model_calibrator!(batch::Batch{A}, analyte::B; kwargs...) where {A, B <: A}
    # cal_id = findfirst(x -> ==(x.analyte, analyte), batch.calibrator)
    # isnothing(cal_id) && throw(ArgumentError("No fitted calibration data for $analyte"))
    model_calibrator!(batch, getanalyteid_check(batch, analyte); kwargs...)
end
function model_calibrator!(batch::Batch, cid::Int; kwargs...) 
    cal_id = getcalibratorid_check(batch, batch.analyte[cid])
    model_calibrator!(batch, cid, batch.calibrator[cal_id]; kwargs...)
end
function model_calibrator!(batch::Batch, cal::AbstractCalibrator; kwargs...) 
    cid = getanalyteid_check(batch, cal.analyte)
    model_calibrator!(batch, cid, cal)
end

function model_calibrator!(batch::Batch, cid::Int, cal::AbstractCalibrator; kwargs...) 
    kwargs = Dict{Symbol, Any}(kwargs)
    model = get!(kwargs, :model, Nothing)
    delete!(kwargs, :model)
    batch.method.analytetable.model[cid] = modifycalmodel(batch.method.analytetable.model[cid], model; kwargs...)
    model_calibrator!(cal, batch.method.analytetable.model[cid])
end

"""
    model_calibrator!(cal::InternalCalibrator, model) -> InternalCalibrator
    model_calibrator!(cal::ExternalCalibrator, model) -> ExternalCalibrator

Update calibration model with `model`, construct calibration machine and calibrate. 
"""
function model_calibrator!(cal::ExternalCalibrator, model)
    if isnothing(model) 
        cal.machine = mkcalmachine(cal)
    else
        cal.model = model
        cal.machine = mkcalmachine(cal)
    end
    # @warn "Call `quantify!` function on data to get the updated results."
    analyze_calibrator!(cal)
end
function model_calibrator!(cal::InternalCalibrator, model)
    cal
end

"""
    edit_method_calibrate!(batch::Batch{A}, args...; kwargs...) -> Batch

Apply `edit_method!` and `calibrate!` to entire batch. See documentation of both funtions for details.
"""
edit_method_calibrate!(batch::Batch, args...; kwargs...) = calibrate!(edit_method!(batch, args...; kwargs...))
"""
    assign_isd_calibrate!(batch::Batch, analyte, isd = nothing; kwargs...) -> Batch

Apply `assign_isd!` and `calibrate!`. See documentation of `assign_isd!` for details.
"""
assign_isd_calibrate!(batch::Batch, args...; kwargs...) = assign_isd!(batch, args...; calibrate = true, kwargs...)

"""
    assign_std_calibrate!(batch::Batch, analyte, cal, isd = isdof(cal, batch.method); kwargs...) -> Batch 

Apply `assign_std!` and `calibrate!`. See documentation of `assign_std!`for details.
"""
assign_std_calibrate!(batch::Batch, args...; kwargs...) = assign_std!(batch, args...; calibrate = true, kwargs...)

@deprecate replace_std_calibrate! assign_std_calibrate!

"""
    quantify_calibrator!(batch::Batch) -> Batch
    quantify_calibrator!(cal::Vector{<: AbstractCalibrator}) -> Vector{<: AbstractCalibrator}
    quantify_calibrator!(cal::InternalCalibrator) -> InternalCalibrator
    quantify_calibrator!(cal::ExternalCalibrator) -> ExternalCalibrator

Inverse predict concentration, update each `cal.table.x̂` with the result(s) and returns `cal` or `batch`.
"""
quantify_calibrator!(batch::Batch) = (quantify_calibrator!(batch.calibrator); batch)
quantify_calibrator!(cal::Vector{<: AbstractCalibrator}) = (foreach(quantify_calibrator!, cal); cal)
quantify_calibrator!(cal::InternalCalibrator) = cal
function quantify_calibrator!(cal::ExternalCalibrator)
    cal.table.x̂ .= inv_predict(cal, cal.table.y)
    cal
end

quantify_calibrator(cal::InternalCalibrator, tbl = nothing) = isnothing(tbl) ? [cal.conc] : repeat([cal.conc], length(tbl))
quantify_calibrator(cal::ExternalCalibrator, tbl = cal.table) = inv_predict(cal, tbl.y)

"""
    validate_calibrator!(batch::Batch) -> Batch
    validate_calibrator!(cal::Vector{<: AbstractCalibrator}) -> Vector{<: AbstractCalibrator}
    validate_calibrator!(cal::InternalCalibrator) -> InternalCalibrator
    validate_calibrator!(cal::ExternalCalibrator) -> ExternalCalibrator

Calculate accuracy for analyte specified by `cal`, update each `cal.table.x̂` with the result(s), and return `cal` or `batch`.
"""
validate_calibrator!(batch::Batch) = (validate_calibrator!(batch.calibrator); batch)
validate_calibrator!(cal::Vector{<: AbstractCalibrator}) = (foreach(validate_calibrator!, cal); cal)
function validate_calibrator!(cal::ExternalCalibrator)
    cal.table.accuracy .= accuracy(cal.table.x̂, cal.table.x)
    cal
end
validate_calibrator!(cal::InternalCalibrator) = cal

validate_calibrator(cal::ExternalCalibrator, tbl = cal.table) = accuracy(tbl.x̂, tbl.x)
validate_calibrator(cal::InternalCalibrator, tbl) = repeat([1.0], length(tbl))

"""
    const analyze_calibrator! = validate_calibrator! ∘ quantify_calibrator!

Apply `quantify_calibrator!` and `validate_calibrator!` subsequantly to `Batch` or `AbstractCalibrator`.
"""
const analyze_calibrator! = validate_calibrator! ∘ quantify_calibrator!

analyze_calibrator(cal::ExternalCalibrator, tbl = cal.table) = accuracy(inv_predict(cal, tbl.y), tbl.x)
analyze_calibrator(cal::InternalCalibrator, tbl) = repeat([1.0], length(tbl))

# Interfacing AbstractCalibrator
"""
    getformula(::Type{<: CurveType}) -> FormulaTerm

Get correspond `FormulaTerm` of `CurveType` for GLM fitting.
"""
getformula(::Type{Proportional}) = @formula(y ~ 0 + x)
getformula(::Type{Linear}) = @formula(y ~ x)
getformula(::Type{QuadraticOrigin}) = @formula(y ~ 0 + x + x ^ 2)
getformula(::Type{Quadratic}) = @formula(y ~ x + x ^ 2)
getformula(::Type{Logarithmic}) = @formula(y ~ log(x))
getformula(::Type{Exponential}) = @formula(log(y) ~ x)
getformula(::Type{Power}) = @formula(log(y) ~ log(x))

"""
    mkcalmodel(model::Type{CalibrationModel}; weight = ConstWeight()) -> CalibrationModel
    mkcalmodel(::Type{Nothing}; kwargs...) -> Nothing

Construct a calibration model. 
"""
mkmodel(model::Type{T}) where {T <: AbstractCalibrationModel} = throw(ArgumentError("`mkmodel` not defined for `$T`"))
function mkcalmodel(model::Type{CalibrationModel{T}}; weight = ConstWeight()) where T
    model(weight)
end
mkcalmodel(::Type{Nothing}; kwargs...) = nothing

"""
    modifycalmodel(old, model; kwargs...) -> CalibrationModel

Modify old model with new model type `model` and new parameters `kwargs`.

If `model` is type `Nothing`, model type remains the same.
"""
function modifycalmodel(old, model::Type{Nothing}; kwargs...)
    isempty(kwargs) ? old : modifycalmodel(old, typeof(old); kwargs...)
end
function modifycalmodel(old::S, model::Type{T}; kwargs...) where {S, T}
    kwargs = collect(kwargs)
    for p in propertynames(old) 
        i = findfirst(x -> first(x) == p, kwargs)
        if isnothing(i)
            push!(kwargs)
        end
    end
    mkcalmodel(model; kwargs...)
end

"""
    mkcalmachine(model::CalibrationModel, tbl)
    mkcalmachine(cal::ExternalCalibrator)

Construct calibration machine from `tbl` (`cal.table`) of type `model` (`cal.model`). 
"""
mkcalmachine(model::Type{T}, tbl) where {T <: AbstractCalibrationModel} = throw(ArgumentError("`mkcalmachine` not defined for `$T`"))
function mkcalmachine(model::CalibrationModel{T}, tbl) where T
    lm1 = lm(getformula(T), tbl[tbl.include]; weights = getweights(model.weight, tbl.x[tbl.include], tbl.y[tbl.include]))
    # if T == Quadratic && lm1.model.pp.beta0[1] == 0
    #     m = hcat(ones(eltype(tbl.x), count(tbl.include)), tbl.x[tbl.include], tbl.x[tbl.include] .^ 2)
    #     sqrtw = diagm(sqrt.(getweights(model.wfn, tbl.x[tbl.include], tbl.y[tbl.include])))
    #     y = tbl.y[tbl.include]
    #     lm1.model.pp.beta0 = (sqrtw * m) \ (sqrtw * y)
    #     GLM.updateμ!(lm1.model.rr, predict(lm1, tbl[tbl.include]))
    # end
    lm1
end
function mkcalmachine(model::CalibrationModel{T}, tbl) where {T <: Union{Exponential, Power}}
    lm1 = lm(getformula(T), tbl[tbl.include])
    wts = getweights(model.weight, tbl.x[tbl.include], tbl.y[tbl.include])
    beta0 = lm1.model.pp.beta0
    fn = getlsqfn(T)
    fit = curve_fit(fn, tbl.x[tbl.include], tbl.y[tbl.include], wts, [exp(beta0[1]), beta0[2]])
    LsqFitMachine(fit, fn, var(tbl.y[tbl.include]))
end
mkcalmachine(cal::ExternalCalibrator) = mkcalmachine(cal.model, cal.table)

"""
    getlsqfn(::Type{Exponential})
    getlsqfn(::Type{Power}) -> Function

Get correspond target function for LsqFit fitting.
"""
function getlsqfn(::Type{Exponential})
    function fn(x, p)
        @. p[1] * exp(p[2] * x)
    end
end
function getlsqfn(::Type{Power})
    function fn(x, p)
        @. p[1] * x ^ p[2]
    end
end
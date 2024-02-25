"""
    relative_signal(batch::Batch, at::AnalysisTable; signal = batch.method.signal)
    relative_signal(method::AnalysisMethod, at::AnalysisTable; signal = method.signal)
    relative_signal(batch::Batch, dt::AbstractDataTable)
    relative_signal(method::AnalysisMethod, dt::AbstractDataTable)

Calculate relative signal using `getproperty(at, signal)` or `dt` as signal data, and return the result as `AbstractDataTable`.
"""
relative_signal(method::AnalysisMethod, at::AnalysisTable; signal = method.signal) = relative_signal(method, getproperty(at, signal))
relative_signal(batch::Batch, at::AnalysisTable; signal = batch.method.signal) = relative_signal(batch.method, at; signal)
function relative_signal(method::AnalysisMethod, dt::AbstractDataTable{A, S, N}) where {A, S, N}
    analytetable = method.analytetable
    aid = [findfirst(==(analyte), analytetable.analyte) for analyte in analyteobj(dt)]
    cs = length(aid) > 10 ? ThreadsX.map(eachindex(aid)) do i
        convert(Vector{N}, analytetable.isd[aid[i]] < 0 ? repeat([N(NaN)], length(sampleobj(dt))) :
                analytetable.isd[aid[i]] == 0 ? getanalyte(dt, i) :
                getanalyte(dt, i) ./ getanalyte(dt, analytetable.analyte[analytetable.isd[aid[i]]]))::Vector{N}
    end : map(eachindex(aid)) do i
        convert(Vector{N}, analytetable.isd[aid[i]] < 0 ? repeat([N(NaN)], length(sampleobj(dt))) :
                analytetable.isd[aid[i]] == 0 ? getanalyte(dt, i) :
                getanalyte(dt, i) ./ getanalyte(dt, analytetable.analyte[analytetable.isd[aid[i]]]))::Vector{N}
    end
    fill_result!(deepcopy(dt), cs::Vector{Vector{N}})
end
relative_signal(batch::Batch, dt::AbstractDataTable) = relative_signal(batch.method, dt)

"""
    set_relative_signal(at::AbstractDataTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal)
    set_relative_signal(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = :relative_signal)

Calculate relative signal using `getproperty(at, signal)` as signal data, update or insert the value into a copy of `at` at index `rel_sig`, and return the copy.
"""
set_relative_signal(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal) = 
    set_relative_signal(at, batch.method; signal, rel_sig)
function set_relative_signal(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = :relative_signal)
    result = relative_signal(method, at; signal)
    new = copy(at)
    set!(new, rel_sig, result)
    new
end

"""
    set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal)
    set_relative_signal!(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = :relative_signal)

Calculate relative signal using `getproperty(at, signal)` as signal data, update and insert the value into `at` at index `rel_sig`, and return `at`.
"""
set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal) = 
    set_relative_signal!(at, batch.method; signal, rel_sig)
function set_relative_signal!(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = :relative_signal)
    result = relative_signal(method, at; signal)
    set!(at, rel_sig, result)
    at
end

"""
    update_relative_signal!(batch::Batch; signal = batch.method.signal)

Calculate relative signal using `getproperty(batch.data, signal)` as signal data, update and insert the value into `batch.data` at index `rel_sig`, and return the updated `batch`.
"""
function update_relative_signal!(batch::Batch{A, M, C, D}; signal = batch.method.signal, rel_sig = :relative_signal) where {A, M, C, D}
    set_relative_signal!(batch.data, batch.method; signal, rel_sig)
    batch
end
function update_relative_signal!(batch::Batch{A, M, C, Nothing}; signal = batch.method.signal, rel_sig = :relative_signal) where {A, M, C}
    @warn "There is no data!"
    batch
end
"""
    inv_predict(batch::Batch, at::AnalysisTable; rel_sig = :relative_signal)
    inv_predict(batch::Batch, dt::AbstractDataTable)

Inversely predict concentration based on relative signal data, `getproperty(at, rel_sig)` or `dt`, and return the result as `AbstractDataTable`.
"""
inv_predict(batch::Batch, at::AnalysisTable; rel_sig = :relative_signal) = 
    inv_predict(batch, getproperty(at, rel_sig))

function inv_predict(batch::Batch, dt::AbstractDataTable{A, S, N}) where {A, S, N}
    analytetable = batch.method.analytetable
    cid = [analytetable.calibration[findfirst(==(analyte), analytetable.analyte)] for analyte in analyteobj(dt)]
    cal_id = [id > 0 ? findfirst(cal -> first(cal.analyte) == analytetable.analyte[id], batch.calibration) : nothing for id in cid]
    cs = length(cal_id) > 10 ? ThreadsX.map(eachindex(cal_id)) do i
        convert(Vector{N}, isnothing(cal_id[i]) ? repeat([N(NaN)], length(sampleobj(dt))) : 
            inv_predict(batch.calibration[cal_id[i]], getanalyte(dt, i)))::Vector{N}
    end : map(eachindex(cal_id)) do i
        convert(Vector{N}, isnothing(cal_id[i]) ? repeat([N(NaN)], length(sampleobj(dt))) : 
            inv_predict(batch.calibration[cal_id[i]], getanalyte(dt, i)))::Vector{N}
    end
    fill_result!(deepcopy(dt), cs::Vector{Vector{N}})
end

"""
    inv_predict(cal::AbstractCalibration, dt::AbstractDataTable; analyte = first(cal.analyte))
    inv_predict(cal::SingleCalibration, y::AbstractArray)
    inv_predict(cal::MultipleCalibration, y::AbstractArray)
    inv_predict(cal::SingleCalibration)
    inv_predict(cal::MultipleCalibration)

Inversely predict concentration of `analyte` or analyte specified in `cal`, and return the result as a vector using data in `dt`, `y` or `cal.table.y` as inverse predictors.
"""
inv_predict(cal::SingleCalibration) = [cal.conc]
inv_predict(cal::MultipleCalibration) = inv_predict(cal, cal.table.y)
inv_predict(cal::AbstractCalibration, dt::AbstractDataTable; analyte = first(cal.analyte)) = inv_predict(cal, getanalyte(dt, analyte))
function inv_predict(cal::MultipleCalibration, y::AbstractArray{T}) where T
    β = convert(Vector{T}, cal.model.model.pp.beta0)::Vector{T}
    if cal.type && cal.zero
        y ./ β[1]
    elseif cal.type
        @. (y - β[1]) / β[2]
    else
        c, b, a = cal.zero ? (zero(T), β...) : β
        d = @. max(b ^ 2 + 4 * a * (y - c), 0)
        @. (-b + sqrt(d)) / 2a
    end
end
inv_predict(cal::SingleCalibration, y::AbstractArray) = y .* cal.conc

"""
    set_inv_predict(at::AnalysisTable, batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)

Inversely predict concantration using `getproperty(at, rel_sig)` as relstive signal data, update or insert the value into a copy of `at` at index `est_conc`, and return the copy.
"""
function set_inv_predict(at::AnalysisTable, batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)
    result = inv_predict(batch, at; rel_sig)
    new = copy(at)
    set!(new, est_conc, result)
    new
end

"""
    set_inv_predict!(at::AnalysisTable, batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)

Inversely predict concantration using `getproperty(at, rel_sig)` as relstive signal data, update or insert the value into `at` at index `est_conc`, and return `at`.
"""
function set_inv_predict!(at::AnalysisTable, batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)
    set!(at, est_conc, inv_predict(batch, at; rel_sig))
    at
end

"""
    update_inv_predict!(batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)

Inversely predict concentration using `getproperty(batch.data, rel_sig)` or `dt` as relstive signal data, update or insert the value into `batch.data` at index `est_conc` and returns the updated `batch`.
"""
function update_inv_predict!(batch::Batch{A, M, C, D}; rel_sig = :relative_signal, est_conc = :estimated_concentration) where {A, M, C, D}
    set!(batch.data, est_conc, inv_predict(batch, batch.data; rel_sig))
    batch
end
function update_inv_predict!(batch::Batch{A, M, C, Nothing}; rel_sig = :relative_signal, est_conc = :estimated_concentration) where {A, M, C}
    @warn "There is no data!"
    batch
end

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
    quantification(batch::Batch, at::AnalysisTable; signal = batch.method.signal)
    quantification(batch::Batch, dt::AbstractDataTable)

Quantify all analyes based on relative signal data, and return the result as `AbstractDataTable` using `getproperty(at, signal)` or `dt` as signal data.
"""
quantification(batch::Batch, at::AnalysisTable; signal = batch.method.signal) = quantification(batch, getproperty(at, signal))
function quantification(batch::Batch, dt::AbstractDataTable{A, S, N}) where {A, S, N}
    analytetable = batch.method.analytetable
    cid = [analytetable.calibration[findfirst(==(analyte), analytetable.analyte)] for analyte in analyteobj(dt)]
    cal_id = [id > 0 ? findfirst(cal -> first(cal.analyte) == analytetable.analyte[id], batch.calibration) : nothing for id in cid]
    cs = length(cal_id) > 10 ? ThreadsX.map(eachindex(cal_id)) do i
        isnothing(cal_id[i]) && return repeat([N(NaN)], length(sampleobj(dt)))
        cal = batch.calibration[cal_id[i]]
        convert(Vector{N}, quantification(cal, dt, i, last(cal.analyte)))::Vector{N}
    end : map(eachindex(cal_id)) do i
        isnothing(cal_id[i]) && return repeat([N(NaN)], length(sampleobj(dt)))
        cal = batch.calibration[cal_id[i]]
        convert(Vector{N}, quantification(cal, dt, i, last(cal.analyte)))::Vector{N}
    end
    fill_result!(deepcopy(dt), cs::Vector{Vector{N}})
end

"""
    quantification(cal::AbstractCalibration, dt::AbstractDataTable; analyte = cal.analyte)

Quantify `analyte` using data in `dt` as signals, and return the result as a vector.
"""
quantification(cal::SingleCalibration) = [cal.conc]
quantification(cal::MultipleCalibration) = inv_predict(cal, cal.table.y)
function quantification(cal::AbstractCalibration, dt::AbstractDataTable{A, S, N}; analyte = cal.analyte) where {A, S, N}
    isnothing(last(analyte)) && return inv_predict(cal, getanalyte(dt, first(analyte))::AbstractVector{N})
    inv_predict(cal, (getanalyte(dt, first(analyte)) ./ getanalyte(dt, last(analyte)))::AbstractVector{N})
end

function quantification(cal::AbstractCalibration, dt::AbstractDataTable{A, S, N}, analyte, isd) where {A, S, N}
    isnothing(isd) && return inv_predict(cal, getanalyte(dt, analyte)::AbstractVector{N})
    inv_predict(cal, (getanalyte(dt, analyte) ./ getanalyte(dt, isd))::AbstractVector{N})
end

"""
    set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)

Quantify all analytes using `getproperty(at, signal)` as signal data., update or insert the values into a copy of `at` at index `rel_sig` for relative signal and `est_conc` for concentration, and return the copy.
"""
function set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)
    new = set_relative_signal(at, batch; signal, rel_sig)
    result = inv_predict(batch, new,; rel_sig)
    set!(new, est_conc, result)
    new
end

"""
    set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)

Quantify all analytes, update or insert the values into `at` at index `rel_sig` for relative signal and `est_conc` for concentration, and return `at`.
"""
function set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)
    set_relative_signal!(at, batch; signal, rel_sig)
    set_inv_predict!(at, batch; rel_sig, est_conc)
end

"""
    update_quantification!(batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)

Quantify all analytes using `getproperty(at, signal)` as signal data, update or insert the values into `batch.data` at index `rel_sig` for relative signal and `est_conc` for concentration, and returns the updated `batch`.
"""
function update_quantification!(batch::Batch{A, M, C, D}; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration) where {A, M, C, D}
    update_relative_signal!(batch; signal, rel_sig)
    update_inv_predict!(batch; rel_sig, est_conc)
end
function update_quantification!(batch::Batch{A, M, C, Nothing}; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration) where {A, M, C}
    @warn "There is no data!"
    batch
end
"""
    accuracy(at::AnalysisTable; true_conc = :true_concentration, est_conc = :estimated_concentration) 
    accuracy(dtp::AbstractDataTable, dtt::AbstractDataTable)
    accuracy(cal::MultipleCalibration, tbl = cal.table)
    accuracy(cal::SingleCalibration, tbl)
    accuracy(x̂::AbstractVector, x::AbstractVector)

Calculate accuracy and return the values as `AbstractDataTable` or `Vector`. `tbl` must contain two properties, `y` and `x`, as same as `cal.table`.
"""
accuracy(at::AnalysisTable; true_conc = :true_concentration, est_conc = :estimated_concentration) = 
    accuracy(getproperty(at, est_conc), getproperty(at, true_conc))

function accuracy(dtp::AbstractDataTable{A, S, N}, dtt::AbstractDataTable) where {A, S, N}
    cs = length(analyteobj(dtp)) > 10 ? ThreadsX.map(analyteobj(dtp)) do analyte
        convert(Vector{N}, accuracy(getanalyte(dtp, analyte), getanalyte(dtt, analyte)))::Vector{N}
    end : map(analyteobj(dtp)) do analyte
        convert(Vector{N}, accuracy(getanalyte(dtp, analyte), getanalyte(dtt, analyte)))::Vector{N}
    end
    fill_result!(deepcopy(dtp), cs::Vector{Vector{N}})
end
accuracy(cal::MultipleCalibration, tbl = cal.table) = accuracy(inv_predict(cal, tbl.y), tbl.x)
accuracy(cal::SingleCalibration, tbl) = accuracy(inv_predict(cal, tbl.y), tbl.x)
accuracy(x̂::AbstractVector, x::AbstractVector) = @. x̂ / x

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
    set_accuracy(at::AnalysisTable; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)

Calculate accuracy, update or insert the values into a copy of `at` at index `acc`, and return the copy. 
using `getproperty(at, true_conc)` as true concentration and `getproperty(at, est_conc)` as estimated concentration.
"""
function set_accuracy(at::AnalysisTable; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)
    result = accuracy(at; true_conc, est_conc)
    new = AnalysisTable(analyteobj(at), sampleobj(at), Dictionary(tables(at)))
    set!(new, acc, result)
    new
end

"""
    set_accuracy!(at::AnalysisTable; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)

Calculate accuracy, update or insert the values into `at` at index `acc`, and return the object 
using `getproperty(at, true_conc)` as true concentration and `getproperty(at, est_conc)` as estimated concentration.
"""
function set_accuracy!(at::AnalysisTable; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)
    set!(at, acc, accuracy(at; true_conc, est_conc))
    at
end

"""
    update_accuracy!(batch::Batch; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)

Calculate accuracy, and update or insert the values into `batch.data` or a copy of `batch.data` at index `acc`, and returns the updated `batch`.
Use `getproperty(batch.data, true_conc)` as true concentration and `getproperty(batch.data, est_conc)` as estimated concentration. 
This function assigns `at` to `batch.data` 
"""
function update_accuracy!(batch::Batch{A, M, C, D}; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy) where {A, M, C, D}
    set!(batch.data, acc, accuracy(batch.data; true_conc, est_conc))
    batch
end
function update_accuracy!(batch::Batch{A, M, C, Nothing}; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy) where {A, M, C}
    @warn "There is no data!"
    batch
end
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
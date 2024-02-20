"""
    relative_signal(batch::Batch, at::AnalysisTable; signal = batch.method.signal)
    relative_signal(method::MethodTable, at::AnalysisTable; signal = method.signal)
    relative_signal(batch::Batch, dt::AbstractDataTable)
    relative_signal(method::MethodTable, dt::AbstractDataTable)

Calculate relative signal, and return the result as `AbstractDataTable` using `getproperty(at, signal)` or `dt` as signal data.
"""
relative_signal(method::MethodTable, at::AnalysisTable; signal = method.signal) = relative_signal(method, getproperty(at, signal))
relative_signal(batch::Batch, at::AnalysisTable; signal = batch.method.signal) = relative_signal(batch.method, at; signal)
function relative_signal(method::MethodTable, dt::AbstractDataTable)
    analytetable = method.analytetable
    aid = [findfirst(==(analyte), analytetable.analyte) for analyte in analyteobj(dt)]
    cs = length(aid) > 10 ? ThreadsX.map(eachindex(aid)) do i
        convert(Vector{Float64}, analytetable.isd[aid[i]] < 0 ? repeat([NaN], length(sampleobj(dt))) :
                analytetable.isd[aid[i]] == 0 ? getanalyte(dt, i) :
                getanalyte(dt, i) ./ getanalyte(dt, analytetable.analyte[analytetable.isd[aid[i]]]))::Vector{Float64}
    end :  map(eachindex(aid)) do i
        convert(Vector{Float64}, analytetable.isd[aid[i]] < 0 ? repeat([NaN], length(sampleobj(dt))) :
                analytetable.isd[aid[i]] == 0 ? getanalyte(dt, i) :
                getanalyte(dt, i) ./ getanalyte(dt, analytetable.analyte[analytetable.isd[aid[i]]]))::Vector{Float64}
    end
    fill_result!(deepcopy(dt), cs::Vector{Vector{Float64}})
end
relative_signal(batch::Batch, dt::AbstractDataTable) = relative_signal(batch.method, dt)

"""
    set_relative_signal(at::AbstractDataTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal)
    set_relative_signal(at::AnalysisTable, method::MethodTable; signal = method.signal, rel_sig = :relative_signal)

Calculate relative signal, update or insert the value into a copy of `at` at index `rel_sig`, and return the object using `getproperty(at, signal)` as signal data.
"""
set_relative_signal(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal) = 
    set_relative_signal(at, batch.method; signal, rel_sig)
function set_relative_signal(at::AnalysisTable, method::MethodTable; signal = method.signal, rel_sig = :relative_signal)
    result = relative_signal(method, at; signal)
    new = AnalysisTable(analyteobj(at), sampleobj(at), deepcopy(tables(at)))
    set!(new, rel_sig, result)
    new
end

"""
    set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal)
    set_relative_signal!(at::AnalysisTable, method::MethodTable; signal = method.signal, rel_sig = :relative_signal)

Calculate relative signal, update and insert the value into `at` at index `rel_sig`, and return the object using `getproperty(at, signal)` as signal data.
"""
set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal) = 
    set_relative_signal!(at, batch.method; signal, rel_sig)
function set_relative_signal!(at::AnalysisTable, method::MethodTable; signal = method.signal, rel_sig = :relative_signal)
    result = relative_signal(method, at; signal)
    set!(at, rel_sig, result)
    at
end

"""
    update_relative_signal!(batch::Batch, at::AnalysisTable = batch.data; signal = batch.method.signal)

Calculate relative signal and update or insert the value into `at` at index `rel_sig` using `getproperty(at, signal)` or `dt` as signal data.
This function assigns `at` to `batch.data` and returns the updated `batch`.
"""
function update_relative_signal!(batch::Batch, at::Union{AnalysisTable, Nothing} = batch.data; signal = batch.method.signal, rel_sig = :relative_signal)
    isnothing(at) && throw(ArgumentError("Data can't be nothing."))
    batch.data = at
    set_relative_signal!(at, batch.method; signal, rel_sig)
    batch
end
"""
    inv_predict(batch::Batch, at::AnalysisTable; rel_sig = :relative_signal)
    inv_predict(batch::Batch, dt::AbstractDataTable)

Inversely predict concentration based on relative signal data, and return the result as `AbstractDataTable` using `getproperty(at, signal)` or `dt` as signal data.
"""
inv_predict(batch::Batch, at::AnalysisTable; rel_sig = :relative_signal) = 
    inv_predict(batch, getproperty(at, rel_sig))

function inv_predict(batch::Batch, dt::AbstractDataTable)
    analytetable = batch.method.analytetable
    cid = [analytetable.calibration[findfirst(==(analyte), analytetable.analyte)] for analyte in analyteobj(dt)]
    cal_id = [id > 0 ? findfirst(cal -> first(cal.analyte) == analytetable.analyte[id], batch.calibration) : nothing for id in cid]
    cs = length(cal_id) > 10 ? ThreadsX.map(eachindex(cal_id)) do i
        convert(Vector{Float64}, isnothing(cal_id[i]) ? repeat([NaN], length(sampleobj(dt))) : 
            inv_predict(batch.calibration[cal_id[i]], getanalyte(dt, i)))::Vector{Float64}
    end : map(eachindex(cal_id)) do i
        convert(Vector{Float64}, isnothing(cal_id[i]) ? repeat([NaN], length(sampleobj(dt))) : 
            inv_predict(batch.calibration[cal_id[i]], getanalyte(dt, i)))::Vector{Float64}
    end
    fill_result!(deepcopy(dt), cs::Vector{Vector{Float64}})
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

Inversely predict concantration, update or insert the value into a copy of `at` at index `est_conc`, and return the object using `getproperty(at, rel_sig)` as relstive signal data.
"""
function set_inv_predict(at::AnalysisTable, batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)
    result = inv_predict(batch, at; rel_sig)
    new = AnalysisTable(analyteobj(at), sampleobj(at), Dictionary(tables(at)))
    set!(new, est_conc, result)
    new
end

"""
    set_inv_predict!(at::AnalysisTable, batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)

Inversely predict concantration, update or insert the value into `at` at index `est_conc`, and return the object using `getproperty(at, rel_sig)` as relstive signal data.
"""
function set_inv_predict!(at::AnalysisTable, batch::Batch; rel_sig = :relative_signal, est_conc = :estimated_concentration)
    set!(at, est_conc, inv_predict(batch, at; rel_sig))
    at
end

"""
    update_inv_predict!(batch::Batch, at::AnalysisTable = batch.data; rel_sig = :relative_signal, est_conc = :estimated_concentration)

Inversely predict concentration and update or insert the value into `at` at index `est_conc` using `getproperty(at, rel_sig)` or `dt` as relstive signal data.
This function assigns `at` to `batch.data` and returns the updated `batch`.
"""
function update_inv_predict!(batch::Batch, at::Union{AnalysisTable, Nothing} = batch.data; rel_sig = :relative_signal, est_conc = :estimated_concentration)
    isnothing(at) && throw(ArgumentError("Data can't be nothing."))
    batch.data = at
    set!(at, est_conc, inv_predict(batch, at; rel_sig))
    batch
end

"""
    inv_predict!(batch::Batch) = (inv_predict!.(batch.calibration); batch)
    inv_predict!(cal::SingleCalibration)
    inv_predict!(cal::MultipleCalibration)

Inverse predict concentration, update each `cal.table.x̂` with the result(s) and returns `cal` or `batch`.
"""
inv_predict!(batch::Batch) = (foreach(inv_predict!, batch.calibration); batch)
inv_predict!(cal::SingleCalibration) = cal
function inv_predict!(cal::MultipleCalibration)
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
function quantification(batch::Batch, dt::AbstractDataTable)
    analytetable = batch.method.analytetable
    cid = [analytetable.calibration[findfirst(==(analyte), analytetable.analyte)] for analyte in analyteobj(dt)]
    cal_id = [id > 0 ? findfirst(cal -> first(cal.analyte) == analytetable.analyte[id], batch.calibration) : nothing for id in cid]
    cs = length(cal_id) > 10 ? ThreadsX.map(eachindex(cal_id)) do i
        isnothing(cal_id[i]) && return repeat([NaN], length(sampleobj(dt)))
        cal = batch.calibration[cal_id[i]]
        convert(Vector{Float64}, quantification(cal, dt, i, last(cal.analyte)))::Vector{Float64}
    end : map(eachindex(cal_id)) do i
        isnothing(cal_id[i]) && return repeat([NaN], length(sampleobj(dt)))
        cal = batch.calibration[cal_id[i]]
        convert(Vector{Float64}, quantification(cal, dt, i, last(cal.analyte)))::Vector{Float64}
    end
    fill_result!(deepcopy(dt), cs::Vector{Vector{Float64}})
end

"""
    quantification(cal::AbstractCalibration, dt::AbstractDataTable; analyte = cal.analyte)

Quantify `analyte`, and return the result as a vector using data in `dt` as signals.
"""
quantification(cal::SingleCalibration) = [cal.conc]
quantification(cal::MultipleCalibration) = inv_predict(cal, cal.table.y)
function quantification(cal::AbstractCalibration, dt::AbstractDataTable; analyte = cal.analyte)
    isnothing(last(analyte)) && return inv_predict(cal, getanalyte(dt, first(analyte)))
    inv_predict(cal, getanalyte(dt, first(analyte)) ./ getanalyte(dt, last(analyte)))
end

function quantification(cal::AbstractCalibration, dt::AbstractDataTable, analyte, isd)
    isnothing(isd) && return inv_predict(cal, getanalyte(dt, analyte))
    inv_predict(cal, getanalyte(dt, analyte) ./ getanalyte(dt, isd))
end

"""
    set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)

Quantify all analytes, update or insert the values into a copy of `at` at index `rel_sig` for relative signal and `est_conc` for concentration, and return the object using `getproperty(at, signal)` as signal data.
"""
function set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)
    new = set_relative_signal(at, batch; signal, rel_sig)
    result = inv_predict(batch, new,; rel_sig)
    set!(new, est_conc, result)
    new
end

"""
    set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)

Quantify all analytes, update or insert the values into `at` at index `rel_sig` for relative signal and `est_conc` for concentration, and return the object using `getproperty(at, signal)` as signal data.
"""
function set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)
    set_relative_signal!(at, batch; signal, rel_sig)
    set_inv_predict!(at, batch; rel_sig, est_conc)
end

"""
    update_quantification!(batch::Batch, at::AnalysisTable = batch.data; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)

Quantify all analytes, and update or insert the values into `at` or a copy of `at` at index `rel_sig` for relative signal and `est_conc` for concentration using `getproperty(at, signal)` as signal data. 
This function assigns `at` to `batch.data` and returns the updated `batch`.

"""
function update_quantification!(batch::Batch, at::AnalysisTable = batch.data; signal = batch.method.signal, rel_sig = :relative_signal, est_conc = :estimated_concentration)
    set_relative_signal!(at, batch; signal, rel_sig)
    update_inv_predict!(batch, at; rel_sig, est_conc)
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

function accuracy(dtp::AbstractDataTable, dtt::AbstractDataTable)
    cs = length(analyteobj(dtp)) > 10 ? ThreadsX.map(analyteobj(dtp)) do analyte
        accuracy(getanalyte(dtp, analyte), getanalyte(dtt, analyte))
    end : map(analyteobj(dtp)) do analyte
        accuracy(getanalyte(dtp, analyte), getanalyte(dtt, analyte))
    end
    fill_result!(deepcopy(dtp), cs::Vector{Vector{Float64}})
end
accuracy(cal::MultipleCalibration, tbl = cal.table) = accuracy(inv_predict(cal, tbl.y), tbl.x)
accuracy(cal::SingleCalibration, tbl) = accuracy(inv_predict(cal, tbl.y), tbl.x)
accuracy(x̂::AbstractVector, x::AbstractVector) = @. x̂ / x

"""
    accuracy!(batch::Batch) = (accuracy!.(batch.calibration); batch)
    accuracy!(cal::MultipleCalibration)
    accuracy!(cal::SingleCalibration)

Calculate accuracy for analyte specified by `cal`, update each `cal.table.x̂` with the result(s), and return `cal` or `batch`.
"""
accuracy!(batch::Batch) = (foreach(accuracy!, batch.calibration); batch)
function accuracy!(cal::MultipleCalibration)
    cal.table.accuracy .= accuracy(cal.table.x̂, cal.table.x)
    cal
end
accuracy!(cal::SingleCalibration) = cal

"""
    set_accuracy(at::AnalysisTable; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)

Calculate accuracy, update or insert the values into a copy of `at` at index `acc`, and return the object 
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
    update_accuracy!(batch::Batch, at::AnalysisTable = batch.data; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)

Calculate accuracy, and update or insert the values into `at` or a copy of `at` at index `acc` 
using `getproperty(at, true_conc)` as true concentration and `getproperty(at, est_conc)` as estimated concentration. 
This function assigns `at` to `batch.data` and returns the updated `batch`.
"""
function update_accuracy!(batch::Batch, at::Union{AnalysisTable, Nothing} = batch.data; true_conc = :true_concentration, est_conc = :estimated_concentration, acc = :accuracy)
    isnothing(at) && throw(ArgumentError("Data can't be nothing."))
    batch.data = at
    set!(at, acc, accuracy(at; true_conc, est_conc))
    batch
end

"""
    const inv_predict_accuracy! = accuracy! ∘ inv_predict!

Apply `inv_predict!` and `accuracy!` subsequantly to `Batch` or `AbstractCalibration`.
"""
const inv_predict_accuracy! = accuracy! ∘ inv_predict!

"""
    calibration(anisd::Tuple, tbl::Table; analyte = 1, isd = 0, type = true, zero = false, weight = 0)
    calibration(batch::Batch{A}, ana::B; id = sample(batch.method.signaltable), isd = isd_of(batch.method, analyte),
                        type = true, zero = false, weight = 0) where {A, B <: A}
    calibration(method::MethodTable{A}, ana::B; id = sample(method.signaltable), isd = isd_of(method, analyte),
                type = true, zero = false, weight = 0) where {A, B <: A}

Create `MultipleCalibration`.

* `anisd`: `Tuple{B <: A, Any}` which will be stored as `analyte`. The first element is the analyte to be quantified, and the second element is its internal standard, which `nothing` means no internal standard.
* `analyte`: `B` which will be stored as first argument of `analyte`.
* `tbl`: `TypedTable.Table` which will be stored as `table`. It is clean-up calibration data for points selection.
* `method`: `MethodTable`, calibration method and data. 
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
    model = calfit(table, f, type, zero, weight)
    table = Table(; id = table.id, level = table.level, y = table.y, x = table.x, x̂ = zeros(Float64, length(table)), accuracy = zeros(Float64, length(table)), include = table.include)
    inv_predict_accuracy!(MultipleCalibration(anisd, type, zero, weight, f, table, model))
end
calibration(batch::Batch{A}, analyte::B; 
                    id = sampleobj(batch.method.signaltable),
                    isd = nothing,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B} = calibration(batch.method, analyte; id, isd, type, zero, weight)
function calibration(method::MethodTable{A}, analyte::B; 
                    id = sampleobj(method.signaltable),
                    isd = missing,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B <: A}

    isd = (ismissing(isd) ? isd_of(method, analyte) : isd)
    ord = sortperm(method.pointlevel)
    level = method.pointlevel[ord]
    conc = getanalyte(method.conctable, analyte)
    ya = getanalyte(method.signaltable, analyte)
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        conc[findsample(method.conctable, Symbol(l))]
                    end, 
                    y, 
                    x̂ = zeros(Float64, length(id)), 
                    accuracy = zeros(Float64, length(id)),
                    include = trues(length(id))
                    )
    f = getformula(type, zero)
    model = calfit(table, f, type, zero, weight)
    inv_predict_accuracy!(MultipleCalibration((analyte, isd), type, zero, Float64(weight), f, table, model))
end
function calibration(method::MethodTable{A}, i::Int; 
                    id = sampleobj(method.signaltable),
                    isd = missing,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where A

    isd = ismissing(isd) ? isd_of(method, analyteobj(method.conctable)[i]) : isd
    ord = sortperm(method.pointlevel)
    level = method.pointlevel[ord]
    conc = getanalyte(method.conctable, i)
    ya = getanalyte(method.signaltable, i)
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        conc[findsample(method.conctable, Symbol(l))]
                    end, 
                    y, 
                    x̂ = zeros(Float64, length(id)), 
                    accuracy = zeros(Float64, length(id)),
                    include = trues(length(id))
                    )
    f = getformula(type, zero)
    model = calfit(table, f, type, zero, weight)
    inv_predict_accuracy!(MultipleCalibration((analyteobj(method.conctable)[i], isd), type, zero, Float64(weight), f, table, model))
end

"""
    calfit(tbl, formula, type, zero, weight)
    calfit(cal::MultipleCalibration)
    calfit!(batch::Batch)
    calfit!(cal::MultipleCalibration)

Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

Field `model` will be mutated for mutating version. Calling `calfit!` on a `Batch` will apply `calfit!` to all calibration curves.

It returns `GLM` object for non-mutating version and input object for mutating version.
"""
function calfit(tbl, formula, type, zero, weight)
    model = lm(formula, tbl[tbl.include]; wts = tbl.x[tbl.include] .^ weight)
    if !type && !zero && model.model.pp.beta0[1] == 0
        m = hcat(ones(Float64, count(tbl.include)), tbl.x[tbl.include], tbl.x[tbl.include] .^ 2)
        sqrtw = diagm(sqrt.(tbl.x[tbl.include] .^ weight))
        y = tbl.y[tbl.include]
        model.model.pp.beta0 = (sqrtw * m) \ (sqrtw * y)
        GLM.updateμ!(model.model.rr, predict(model, tbl[tbl.include]))
    end
    model
end
calfit(cal::MultipleCalibration) = calfit(cal.table, cal.formula, cal.type, cal.zero, cal.weight)

"""
    calfit(tbl, formula, type, zero, weight)
    calfit(cal::MultipleCalibration)
    calfit!(batch::Batch)
    calfit!(cal::MultipleCalibration)

Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

Field `model` will be mutated for mutating version. Calling `calfit!` on a `Batch` will apply `calfit!` to all calibration curves.

It returns `GLM` object for non-mutating version and input object for mutating version.
"""
calfit!(batch::Batch) = (calfit!.(batch.calibration); batch)
function calfit!(cal::MultipleCalibration)
    cal.model = calfit(cal.table, cal.formula, cal.type, cal.zero, cal.weight)
    cal
end

"""
    update_calibration!(batch::Batch{A}, analyte::B) where {A, B <: A}
    update_calibration!(batch::Batch, cal_id::Int)
    update_calibration!(cal::MultipleCalibration, method::MethodTable)

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
function update_calibration!(cal::MultipleCalibration, method::MethodTable)
    isd = isd_of(method, first(cal.analyte))
    cal.analyte = (first(cal.analyte), isd)
    ord = sortperm(method.pointlevel)
    ya = getanalyte(method.signaltable, first(cal.analyte))
    yi = isnothing(isd) ? [1.0] : getanalyte(method.signaltable, isd)
    y = @. (ya / yi)[ord]
    cal.table.y .= y
    # cal.table.include .= true
    cal.formula = getformula(cal)
    calfit!(cal)
    inv_predict_accuracy!(cal)
end
function update_calibration!(cal::SingleCalibration, method::MethodTable)
    cal.conc = first(getanalyte(method.conctable, first(cal.analyte)))
    cal
end

"""
    set_isd!(batch::Batch{A}, analyte::B, isd::C) where {A, B <: A, C <: A}

Set internal standard of `analyte` to `isd`.
"""
function set_isd!(batch::Batch{A}, analyte::B, isd::Union{C, Nothing} = nothing) where {A, B <: A, C <: A}
    aid = findfirst(==(analyte), batch.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in this batch"))
    ca = batch.analyte[aid]
    cid = findfirst(x -> ==(first(x.analyte), ca), batch.calibration)
    isnothing(cid) && throw(ArgumentError("No fitted calibration data for $analyte"))
    set_isd!(batch.method, analyte, isd)
    update_calibration!(batch, cid)
    @warn "Call quantification function on data to get the updated results."
    batch
end
function set_isd!(method::MethodTable{A}, analyte::B, isd::Union{C, Nothing} = nothing) where {A, B <: A, C <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    isnothing(isd) && (method.analytetable.isd[aid] = 0; return method)
    iid = findfirst(==(isd), method.analytetable.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    method.analytetable.isd[aid] = iid
    method
end
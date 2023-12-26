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
    aid = [findfirst(==(analyte), analytetable.analyte) for analyte in dt.analyte]
    cs = ThreadsX.map(eachindex(aid)) do i
        analytetable.isd[aid[i]] < 0 ? repeat([NaN], length(dt.sample)) :
                analytetable.isd[aid[i]] == 0 ? getanalyte(dt, dt.analyte[i]) :
                getanalyte(dt, dt.analyte[i]) ./ getanalyte(dt, analytetable.analyte[analytetable.isd[aid[i]]])
    end
    fill_result!(deepcopy(dt), cs)
end
relative_signal(batch::Batch, dt::AbstractDataTable) = relative_signal(batch.method, dt)

"""
    set_relative_signal(at::AbstractDataTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal)
    set_relative_signal(at::AnalysisTable, method::MethodTable; signal = method.signal, relsig = :relative_signal)
    set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal)
    set_relative_signal!(at::AnalysisTable, method::MethodTable; signal = method.signal, relsig = :relative_signal)

Calculate relative signal, update or insert the value into `at` or a copy of `at` at index `relsig`, and return the object using `getproperty(at, signal)` as signal data.
"""
set_relative_signal(at::AnalysisTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal) = 
    set_relative_signal(at, batch.method; signal, relsig)
function set_relative_signal(at::AnalysisTable, method::MethodTable; signal = method.signal, relsig = :relative_signal)
    result = relative_signal(method, at; signal)
    new = AnalysisTable(at.analyte, at.sample, Dictionary(at.tables))
    set!(new.tables, relsig, result)
    new
end
set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal) = 
    set_relative_signal!(at, batch.method; signal, relsig)
function set_relative_signal!(at::AnalysisTable, method::MethodTable; signal = method.signal, relsig = :relative_signal)
    result = relative_signal(method, at; signal)
    set!(at.tables, relsig, result)
    at
end

"""
    update_relative_signal!(batch::Batch, at::AnalysisTable = batch.data; signal = batch.method.signal)

Calculate relative signal and update or insert the value into `at` at index `relsig` using `getproperty(at, signal)` or `dt` as signal data.
This function assigns `at` to `batch.data` and returns the updated `batch`.
"""
function update_relative_signal!(batch::Batch, at::AnalysisTable = batch.data; signal = batch.method.signal, relsig = :relative_signal)
    batch.data = at
    set_relative_signal!(at, batch.method; signal, relsig)
    batch
end
"""
    inv_predict(batch::Batch, at::AnalysisTable; relsig = :relative_signal)
    inv_predict(batch::Batch, dt::AbstractDataTable)

Inversely predict concentration based on relative signal data, and return the result as `AbstractDataTable` using `getproperty(at, signal)` or `dt` as signal data.
"""
inv_predict(batch::Batch, at::AnalysisTable; relsig = :relative_signal) = 
    inv_predict(batch, getproperty(at, relsig))

function inv_predict(batch::Batch, dt::AbstractDataTable)
    analytetable = batch.method.analytetable
    cid = [analytetable.calibration[findfirst(==(analyte), analytetable.analyte)] for analyte in dt.analyte]
    cal_id = [id > 0 ? findfirst(cal -> first(cal.analyte) == analytetable.analyte[id], batch.calibration) : nothing for id in cid]
    cs = ThreadsX.map(eachindex(cal_id)) do i
        isnothing(cal_id[i]) ? repeat([NaN], length(dt.sample)) : 
            inv_predict(batch.calibration[cal_id[i]], getanalyte(dt, dt.analyte[i]))
    end
    fill_result!(deepcopy(dt), cs)
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
function inv_predict(cal::MultipleCalibration, y::AbstractArray)
    β = cal.model.model.pp.beta0
    if cal.type && cal.zero
        y ./ β[1]
    elseif cal.type
        (y .- β[1]) ./ β[2]
    else
        c, b, a = cal.zero ? (0, β...) : β
        d = @. max(b ^ 2 + 4 * a * (y - c), 0)
        @. (-b + sqrt(d)) / 2a
    end
end
inv_predict(cal::SingleCalibration, y::AbstractArray) = y .* cal.conc

"""
    set_inv_predict(at::AnalysisTable, batch::Batch; relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    set_inv_predict!(at::AnalysisTable, batch::Batch; relsig = :relative_signal, estimated_concentration = :estimated_concentration)

Inversely predict concantration, update or insert the value into `at` or a copy of `at` at index `estimated_concentration`, and return the object using `getproperty(at, relsig)` as relstive signal data.
"""
function set_inv_predict(at::AnalysisTable, batch::Batch; relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    result = inv_predict(batch, at; relsig)
    new = AnalysisTable(at.analyte, at.sample, Dictionary(at.tables))
    set!(new.tables, estimated_concentration, result)
    new
end
function set_inv_predict!(at::AnalysisTable, batch::Batch; relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    set!(at.tables, estimated_concentration, inv_predict(batch, at; relsig))
    at
end

"""
    update_inv_predict!(batch::Batch, at::AnalysisTable = batch.data; relsig = :relative_signal, estimated_concentration = :estimated_concentration)

Inversely predict concentration and update or insert the value into `at` at index `estimated_concentration` using `getproperty(at, relsig)` or `dt` as relstive signal data.
This function assigns `at` to `batch.data` and returns the updated `batch`.
"""
function update_inv_predict!(batch::Batch, at::AnalysisTable = batch.data; relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    batch.data = at
    set!(at.tables, estimated_concentration, inv_predict(batch, at; relsig))
    batch
end

"""
    inv_predict!(batch::Batch) = (inv_predict!.(batch.calibration); batch)
    inv_predict!(cal::SingleCalibration)
    inv_predict!(cal::MultipleCalibration)

Inverse predict concentration, update each `cal.table.x̂` with the result(s) and returns `cal` or `batch`.
"""
inv_predict!(batch::Batch) = (inv_predict!.(batch.calibration); batch)
inv_predict!(cal::SingleCalibration) = cal
function inv_predict!(cal::MultipleCalibration)
    cal.table.x̂ .= inv_predict(cal, cal.table.y)
    cal
end

function fill_result!(dt::ColumnDataTable, result::Vector{Vector{Float64}})
    for (a, c) in zip(eachanalyte(dt), result)
        a .= c
    end
    dt
end
function fill_result!(dt::RowDataTable, result::Vector{Vector{Float64}})
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
    cid = [analytetable.calibration[findfirst(==(analyte), analytetable.analyte)] for analyte in dt.analyte]
    cal_id = [id > 0 ? findfirst(cal -> first(cal.analyte) == analytetable.analyte[id], batch.calibration) : nothing for id in cid]
    cs = ThreadsX.map(eachindex(cal_id)) do i
        isnothing(cal_id[i]) && return repeat([NaN], length(dt.sample))
        cal = batch.calibration[cal_id[i]]
        quantification(cal, dt; analyte = (dt.analyte[i], last(cal.analyte)))
    end
    fill_result!(deepcopy(dt), cs)
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

"""
    set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal, estimated_concentration = :estimated_concentration)

Quantify all analytes, update or insert the values into `at` or a copy of `at` at index `relsig` for relative signal and `estimated_concentration` for concentration, and return the object using `getproperty(at, signal)` as signal data.
"""
function set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    new = set_relative_signal(at, batch; signal, relsig)
    result = inv_predict(batch, new,; relsig)
    set!(new.tables, estimated_concentration, result)
    new
end
function set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    set_relative_signal!(at, batch; signal, relsig)
    set_inv_predict!(at, batch; relsig, estimated_concentration)
end

"""
    update_quantification!(batch::Batch, at::AnalysisTable = batch.data; signal = batch.method.signal, relsig = :relative_signal, estimated_concentration = :estimated_concentration)

Quantify all analytes, and update or insert the values into `at` or a copy of `at` at index `relsig` for relative signal and `estimated_concentration` for concentration using `getproperty(at, signal)` as signal data. 
This function assigns `at` to `batch.data` and returns the updated `batch`.

"""
function update_quantification!(batch::Batch, at::AnalysisTable = batch.data; signal = batch.method.signal, relsig = :relative_signal, estimated_concentration = :estimated_concentration)
    set_relative_signal!(at, batch; signal, relsig)
    update_inv_predict!(batch, at; relsig, estimated_concentration)
end

"""
    accuracy(at::AnalysisTable; true_concentration = :true_concentration, estimated_concentration = :estimated_concentration) 
    accuracy(dtp::AbstractDataTable, dtt::AbstractDataTable)
    accuracy(cal::MultipleCalibration, tbl = cal.table)
    accuracy(cal::SingleCalibration, tbl)
    accuracy(x̂::AbstractVector, x::AbstractVector)

Calculate accuracy and return the values as `AbstractDataTable` or `Vector`. `tbl` must contain two properties, `y` and `x`, as same as `cal.table`.
"""
accuracy(at::AnalysisTable; true_concentration = :true_concentration, estimated_concentration = :estimated_concentration) = 
    accuracy(getproperty(at, estimated_concentration), getproperty(at, true_concentration))

function accuracy(dtp::AbstractDataTable, dtt::AbstractDataTable)
    cs = ThreadsX.map(dtp.analyte) do analyte
        accuracy(getanalyte(dtp, analyte), getanalyte(dtt, analyte))
    end
    fill_result!(deepcopy(dtp), cs)
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
accuracy!(batch::Batch) = (accuracy!.(batch.calibration); batch)
function accuracy!(cal::MultipleCalibration)
    cal.table.accuracy .= accuracy(cal.table.x̂, cal.table.x)
    cal
end
accuracy!(cal::SingleCalibration) = cal

"""
    set_accuracy(at::AnalysisTable; true_concentration = :true_concentration, estimated_concentration = :estimated_concentration, acc = :accuracy)

Calculate accuracy, update or insert the values into `at` or a copy of `at` at index `acc`, and return the object 
using `getproperty(at, true_concentration)` as true concentration and `getproperty(at, estimated_concentration)` as estimated concentration.
"""
function set_accuracy(at::AnalysisTable; true_concentration = :true_concentration, estimated_concentration = :estimated_concentration, acc = :accuracy)
    result = accuracy(at; true_concentration, estimated_concentration)
    new = AnalysisTable(at.analyte, at.sample, Dictionary(at.tables))
    set!(new.tables, acc, result)
    new
end
function set_accuracy!(at::AnalysisTable; true_concentration = :true_concentration, estimated_concentration = :estimated_concentration, acc = :accuracy)
    set!(at.tables, acc, accuracy(at; true_concentration, estimated_concentration))
    at
end

"""
    update_accuracy!(batch::Batch, at::AnalysisTable = batch.data; true_concentration = :true_concentration, estimated_concentration = :estimated_concentration, acc = :accuracy)

Calculate accuracy, and update or insert the values into `at` or a copy of `at` at index `acc` 
using `getproperty(at, true_concentration)` as true concentration and `getproperty(at, estimated_concentration)` as estimated concentration. 
This function assigns `at` to `batch.data` and returns the updated `batch`.
"""
function update_accuracy!(batch::Batch, at::AnalysisTable = batch.data; true_concentration = :true_concentration, estimated_concentration = :estimated_concentration, acc = :accuracy)
    batch.data = at
    set!(at.tables, acc, accuracy(at; true_concentration, estimated_concentration))
    batch
end

"""
    const inv_predict_accuracy! = accuracy! ∘ inv_predict!

Apply `inv_predict!` and `accuracy!` subsequantly to `Batch` or `AbstractCalibration`.
"""
const inv_predict_accuracy! = accuracy! ∘ inv_predict!

"""
    calibration(analyte::Tuple, tbl::Table; analyte = 1, isd = 0, type = true, zero = false, weight = 0)
    calibration(batch::Batch{A}, analyte::B; id = batch.method.signaltable.sample, isd = isd_of(batch.method, analyte),
                        type = true, zero = false, weight = 0) where {A, B <: A}
    calibration(method::MethodTable{A}, analyte::B; id = method.signaltable.sample, isd = isd_of(method, analyte),
                type = true, zero = false, weight = 0) where {A, B <: A}

Create `MultipleCalibration`.

* `tbl`: `TypedTable.Table` which will be stored as `table`. It is clean-up calibration data for points selection.
* `method`: `MethodTable`, calibration method and data. 
* `batch`: `Batch`.

# Keyword arguments
* `analyte`: `Tuple{B <: A, Any}` which will be stored as `analyte`. The first element is the analyte to be quantified, and the second element is its internal standard, which `nothing` means no internal standard.
* `type`: `Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`), which will be stored as `type`.
* `zero`: `Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`), which will be stored as `zero`.
* `weight`: `Float64` represents the exponential applying to each element of `x` as a weighting vector, which will be stored as `weight`.
"""
function calibration(analyte::Tuple, tbl::Table;
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
    inv_predict_accuracy!(MultipleCalibration(analyte, type, zero, weight, f, table, model))
end
calibration(batch::Batch{A}, analyte::B; 
                    id = batch.method.signaltable.sample,
                    isd = isd_of(batch.method, analyte),
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B <: A} = calibration(batch.method, analyte; id, isd, type, zero, weight)
function calibration(method::MethodTable{A}, analyte::B; 
                    id = method.signaltable.sample,
                    isd = isd_of(method, analyte),
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B <: A}

    ord = sortperm(method.pointlevel)
    level = method.pointlevel[ord]
    conc = getanalyte(method.conctable, analyte)
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        conc[findsample(method.conctable, Symbol(l))]
                    end, 
                    y = isnothing(isd) ? getanalyte(method.signaltable, analyte)[ord] : (getanalyte(method.signaltable, analyte) ./ getanalyte(method.signaltable, isd))[ord], 
                    x̂ = zeros(Float64, length(id)), 
                    accuracy = zeros(Float64, length(id)),
                    include = trues(length(id)))
    f = getformula(type, zero)
    model = calfit(table, f, type, zero, weight)
    inv_predict_accuracy!(MultipleCalibration((analyte, isd), type, zero, Float64(weight), f, table, model))
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
    aid = findfirst(==(analyte), batch.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in this batch"))
    ca = batch.analyte[aid]
    cid = findfirst(x -> ==(first(x.analyte), ca), batch.calibration)
    isnothing(cid) && throw(ArgumentError("No fitted calibration data for $analyte"))
    update_calibration!(batch, cid)
end
update_calibration!(batch::Batch, cal_id::Int) = update_calibration!(batch.calibration[cal_id], batch.method)
function update_calibration!(cal::MultipleCalibration, method::MethodTable)
    isd = isd_of(method, first(cal.analyte))
    cal.analyte = (first(cal.analyte), isd)
    ord = sortperm(method.pointlevel)
    cal.table.y .= isnothing(isd) ? getanalyte(method.signaltable, first(cal.analyte))[ord] : (getanalyte(method.signaltable, first(cal.analyte)) ./ getanalyte(method.signaltable, isd))[ord]
    cal.table.include .= true
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
    aid = findfirst(==(analyte), method.analytemap.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    isnothing(isd) && (method.analytetable.isd[aid] = 0; return method)
    iid = findfirst(==(isd), method.analytemap.analyte)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    method.analytetable.isd[aid] = iid
    method
end
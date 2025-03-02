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
    set_relative_signal(at::AbstractDataTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig)
    set_relative_signal(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = method.rel_sig)

Calculate relative signal using `getproperty(at, signal)` as signal data, update or insert the value into a copy of `at` at index `rel_sig`, and return the copy.
"""
set_relative_signal(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig) = 
    set_relative_signal(at, batch.method; signal, rel_sig)
function set_relative_signal(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = method.rel_sig)
    signal == rel_sig && return copy(at)
    result = relative_signal(method, at; signal)
    new = copy(at)
    set!(new, rel_sig, result)
    new
end

"""
    set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig)
    set_relative_signal!(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = method.rel_sig)

Calculate relative signal using `getproperty(at, signal)` as signal data, update and insert the value into `at` at index `rel_sig`, and return `at`.
"""
set_relative_signal!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig) = 
    set_relative_signal!(at, batch.method; signal, rel_sig)
function set_relative_signal!(at::AnalysisTable, method::AnalysisMethod; signal = method.signal, rel_sig = method.rel_sig)
    signal == rel_sig && return at
    result = relative_signal(method, at; signal)
    set!(at, rel_sig, result)
    at
end

"""
    update_relative_signal!(batch::Batch; signal = batch.method.signal)

Calculate relative signal using `getproperty(batch.data, signal)` as signal data, update and insert the value into `batch.data` at index `rel_sig`, and return the updated `batch`.
"""
function update_relative_signal!(batch::Batch{A, M, C, D}; signal = batch.method.signal, rel_sig = batch.method.rel_sig) where {A, M, C, D}
    signal == rel_sig || set_relative_signal!(batch.data, batch.method; signal, rel_sig)
    batch
end
function update_relative_signal!(batch::Batch{A, M, C, Nothing}; signal = batch.method.signal, rel_sig = batch.method.rel_sig) where {A, M, C}
    @warn "There is no data!"
    batch
end
"""
    inv_predict(batch::Batch, at::AnalysisTable; rel_sig = batch.method.rel_sig)
    inv_predict(batch::Batch, dt::AbstractDataTable)

Inversely predict concentration based on relative signal data, `getproperty(at, rel_sig)` or `dt`, and return the result as `AbstractDataTable`.
"""
inv_predict(batch::Batch, at::AnalysisTable; rel_sig = batch.method.rel_sig) = 
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
    set_inv_predict(at::AnalysisTable, batch::Batch; rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)

Inversely predict concantration using `getproperty(at, rel_sig)` as relstive signal data, update or insert the value into a copy of `at` at index `est_conc`, and return the copy.
"""
function set_inv_predict(at::AnalysisTable, batch::Batch; rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)
    result = inv_predict(batch, at; rel_sig)
    new = copy(at)
    set!(new, est_conc, result)
    new
end

"""
    set_inv_predict!(at::AnalysisTable, batch::Batch; rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)

Inversely predict concantration using `getproperty(at, rel_sig)` as relstive signal data, update or insert the value into `at` at index `est_conc`, and return `at`.
"""
function set_inv_predict!(at::AnalysisTable, batch::Batch; rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)
    set!(at, est_conc, inv_predict(batch, at; rel_sig))
    at
end

"""
    update_inv_predict!(batch::Batch; rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)

Inversely predict concentration using `getproperty(batch.data, rel_sig)` or `dt` as relstive signal data, update or insert the value into `batch.data` at index `est_conc` and returns the updated `batch`.
"""
function update_inv_predict!(batch::Batch{A, M, C, D}; rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc) where {A, M, C, D}
    set!(batch.data, est_conc, inv_predict(batch, batch.data; rel_sig))
    batch
end
function update_inv_predict!(batch::Batch{A, M, C, Nothing}; rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc) where {A, M, C}
    @warn "There is no data!"
    batch
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
    set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)

Quantify all analytes using `getproperty(at, signal)` as signal data., update or insert the values into a copy of `at` at index `rel_sig` for relative signal and `est_conc` for concentration, and return the copy.
"""
function set_quantification(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)
    signal == rel_sig && return set_inv_predict(at, batch; rel_sig, est_conc)
    new = set_relative_signal(at, batch; signal, rel_sig)
    result = inv_predict(batch, new,; rel_sig)
    set!(new, est_conc, result)
    new
end

"""
    set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)

Quantify all analytes, update or insert the values into `at` at index `rel_sig` for relative signal and `est_conc` for concentration, and return `at`.
"""
function set_quantification!(at::AnalysisTable, batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)
    signal == rel_sig || set_relative_signal!(at, batch; signal, rel_sig)
    set_inv_predict!(at, batch; rel_sig, est_conc)
end

"""
    update_quantification!(batch::Batch; signal = batch.method.signal, rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc)

Quantify all analytes using `getproperty(at, signal)` as signal data, update or insert the values into `batch.data` at index `rel_sig` for relative signal and `est_conc` for concentration, and returns the updated `batch`.
"""
function update_quantification!(batch::Batch{A, M, C, D}; signal = batch.method.signal, rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc) where {A, M, C, D}
    signal == rel_sig || update_relative_signal!(batch; signal, rel_sig)
    update_inv_predict!(batch; rel_sig, est_conc)
end
function update_quantification!(batch::Batch{A, M, C, Nothing}; signal = batch.method.signal, rel_sig = batch.method.rel_sig, est_conc = batch.method.est_conc) where {A, M, C}
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

accuracy(at::AnalysisTable, method::AnalysisMethod; true_conc = method.true_conc, est_conc = method.est_conc) = 
    accuracy(getproperty(at, est_conc), getproperty(at, true_conc))

accuracy(at::AnalysisTable, batch::Batch; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc) = 
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
    set_accuracy(at::AnalysisTable, method::AnalysisMethod; true_conc = method.true_conc, est_conc = method.est_conc, acc = method.acc)
    set_accuracy(at::AnalysisTable, batch::Batch; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc, acc = batch.method.acc)

Calculate accuracy, update or insert the values into a copy of `at` at index `acc`, and return the copy. 
using `getproperty(at, true_conc)` as true concentration and `getproperty(at, est_conc)` as estimated concentration.
"""
function set_accuracy(at::AnalysisTable, method::AnalysisMethod; true_conc = method.true_conc, est_conc = method.est_conc, acc = method.acc)
    result = accuracy(at; true_conc, est_conc)
    new = AnalysisTable(analyteobj(at), sampleobj(at), Dictionary(tables(at)))
    set!(new, acc, result)
    new
end
set_accuracy(at::AnalysisTable, batch::Batch; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc, acc = batch.method.acc) = 
    set_accuracy(at, batch.method; true_conc, est_conc, acc)

"""
    set_accuracy!(at::AnalysisTable, method::AnalysisMethod; true_conc = method.true_conc, est_conc = method.est_conc, acc = method.acc)
    set_accuracy!(at::AnalysisTable, batch::Batch; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc, acc = batch.method.acc)

Calculate accuracy, update or insert the values into `at` at index `acc`, and return the object 
using `getproperty(at, true_conc)` as true concentration and `getproperty(at, est_conc)` as estimated concentration.
"""
function set_accuracy!(at::AnalysisTable, method::AnalysisMethod; true_conc = method.true_conc, est_conc = method.est_conc, acc = method.acc)
    set!(at, acc, accuracy(at; true_conc, est_conc))
    at
end

set_accuracy!(at::AnalysisTable, batch::Batch; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc, acc = batch.method.acc) = 
    set_accuracy!(at, batch.method; true_conc, est_conc, acc)

"""
    update_accuracy!(batch::Batch; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc, acc = batch.method.acc)

Calculate accuracy, and update or insert the values into `batch.data` or a copy of `batch.data` at index `acc`, and returns the updated `batch`.
Use `getproperty(batch.data, true_conc)` as true concentration and `getproperty(batch.data, est_conc)` as estimated concentration. 
This function assigns `at` to `batch.data` 
"""
function update_accuracy!(batch::Batch{A, M, C, D}; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc, acc = batch.method.acc) where {A, M, C, D}
    set_accuracy!(batch.data, batch; true_conc, est_conc, acc)
    batch
end
function update_accuracy!(batch::Batch{A, M, C, Nothing}; true_conc = batch.method.true_conc, est_conc = batch.method.est_conc, acc = batch.method.acc) where {A, M, C}
    @warn "There is no data!"
    batch
end
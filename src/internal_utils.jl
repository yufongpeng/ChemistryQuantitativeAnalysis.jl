"""
    cal_params(params)

Extract `analyte_std`, `isd`, `std_isd` from `params`.
"""
function cal_params(params)
    params = Table(params)
    pn = propertynames(params)
    analyte = :analyte in pn ? params.analyte : collect(eachindex(params))
    visd = []
    if :std in pn 
        std = params.std
        id = findall(x -> !(x isa Int && x < 0), std)
        analyte_std = [a => s for (a, s) in zip(analyte[id], std[id])]
        params = length(pn) == 1 ? nothing : Table(params; std = nothing)
        visd = vcat(visd, setdiff(eachindex(std), id))
    else
        analyte_std = Pair[]
    end
    if :isd in pn
        isd = params.isd
        id = findall(x -> !(x isa Int && x < 0), isd)
        std_isd = [a => i for (a, i) in zip(analyte[id], isd[id])]
        params = length(pn) == 1 ? nothing : Table(params; isd = nothing)
        visd = vcat(visd, setdiff(eachindex(isd), id))
    else
        std_isd = Pair[]
    end
    analyte_std, unique(visd), std_isd, params
end

"""
    model_params(params; kwargs...)

Extract `model` and vectorize `kwargs`.
"""
function model_params(params, analytetable; kwargs...)
    params = Table(params)
    pn = propertynames(params)
    analyte = :analyte in pn ? params.analyte : collect(eachindex(params))
    add = []
    n = length(params)
    for (k, v) in kwargs
        if !in(k, pn)
            push!(add, k => repeat([v], n))
        end 
    end
    params = Table(params; analyte, add...)
    pt = propertynames(params)
    tid = findall(x -> x in propertynames(analytetable), pt)
    did = findall(x -> x in [:analyte, :model], pt)
    tid = setdiff(tid, did)
    mid = setdiff(eachindex(pt), tid, did)
    (
        analyte, 
        :model in pt ? params.model : [Nothing for _ in eachindex(params)], 
        isempty(tid) ? [Pair[] for _ in 1:n] : getproperties(params, pt[tid]),
        isempty(mid) ? [Pair[] for _ in 1:n] : getproperties(params, pt[mid])
    )
end
function model_params(params::Nothing, analytetable; kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    model = get!(kwargs, :model, Nothing)
    delete!(kwargs, :model)
    param = Dict{Symbol, Any}()
    del = Symbol[]
    for (k, v) in kwargs
        if k in propertynames(analytetable)
            param[k] = v 
            push!(del, k)
        end
    end
    for k in del 
        delete!(kwargs, k)
    end
    model, param, kwargs
end

"""
    getanalyteid(batch::Batch, id::Int)
    getanalyteid(method::AnalysisMethod, id::Int)
    getanalyteid(batch::Batch, analyte)
    getanalyteid(method::AnalysisMethod, analyte)

Return the index (`id`) of `analyte` in `method.analyte` or `batch.analyte`.
"""
getanalyteid(batch::Batch, id::Int) = id
getanalyteid(method::AnalysisMethod, id::Int) = id
getanalyteid(batch::Batch, analyte) = getanalyteid(batch.method, analyte)
getanalyteid(method::AnalysisMethod, analyte) = findfirst(==(analyte), method.analytetable.analyte)

"""
    getcalibratorid(batch::Batch, id::Int)
    getcalibratorid(batch::Batch, analyte)

Return the index of the first element of `batch.calibrator` for which the element is the calibrator of `analyte` or `batch.analyte[id]`.
"""
getcalibratorid(batch::Batch, id::Int) = getcalibratorid(batch, batch.analyte[id])
getcalibratorid(batch::Batch, analyte) = findfirst(x -> ==(x.analyte, analyte), batch.calibrator)

"""
    getanalyteid_check(x, analyte)

`getanalyteid` but throw error if nothing is found.
"""
function getanalyteid_check(x, analyte)
    aid = getanalyteid(x, analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    aid
end
"""
    getisdid_check(x, isd)

`getanalyteid` but throw error if nothing is found and return 0 when `isd` is `nothing`.  
"""
function getisdid_check(x, isd)
    iid = isnothing(isd) ? 0 : getanalyteid(x, isd)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the method"))
    iid
end
"""
    getcalibratorid_check(x, analyte)

`getcalibratorid` but throw error if nothing is found.
"""
function getcalibratorid_check(x, analyte)
    cid = getcalibratorid(x, analyte)
    isnothing(cid) && throw(ArgumentError("Analyte $analyte does not have calibration data"))
    cid
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

# function critical_point(cal::ExternalCalibrator{A, N}) where {A, N}
#     β = cal.model.model.pp.beta0::Vector{N}
#     c, b, a = cal.zero ? (0, β...) : β
#     -b / 2a
# end

function parse_calibration_level_name(dt::AbstractDataTable, calid::Regex, order, clevel, ratio, df, f2c, parse_decimal)
    so = split(order, "")
    dilutionid = findfirst(==("D"), so)
    ratioid = findfirst(==("R"), so)
    levelid = findfirst(==("L"), so)
    cs = map(samplename(dt)) do s
            m = match(calid, string(s))
            isnothing(m) ? missing : collect(String, m)
    end
    if isnothing(clevel)
        isnothing(levelid) && throw(ArgumentError("No valid level id."))
        clevel = map(s -> ismissing(s) ? missing : parse(Int, s[levelid]), cs)
    end
    id = findall(!ismissing, clevel)
    pointlevel = convert(Vector{Int}, clevel[id])
    levels = unique(pointlevel)
    idx = [findfirst(==(x), clevel) for x in levels]
    if !isnothing(ratio)
        concs = map(s -> s .* f2c, ratio)
    elseif !isnothing(df)
        concs = map(s -> f2c ./ s, df)
    elseif !isnothing(ratioid)
        concs = map(s -> parse(Float64, parse_decimal(s[ratioid])) .* f2c, cs[idx])
    elseif !isnothing(dilutionid)
        concs = map(s -> f2c ./ parse(Float64, parse_decimal(s[dilutionid])), cs[idx])
    else
        throw(ArgumentError("No valid ratios or dilution factors are obtained."))
    end
    id, pointlevel, levels, concs
end

# function parse_calibration_level_name(dt::AbstractDataTable, calid, order, clevel, ratio, df, f2c, parse_decimal)
#     clevel = isnothing(clevel) ? replace(calid, missing => nothing) : clevel
#     id = findall(!isnothing, clevel)
#     pointlevel = convert(Vector{Int}, clevel[id])
#     levels = unique(pointlevel)
#     if !isnothing(ratio)
#         concs = map(s -> s .* f2c, ratio)
#     elseif !isnothing(df)
#         concs = map(s -> f2c ./ s, df)
#     else
#         throw(ArgumentError("No valid ratios or dilution factors are obtained."))
#     end
#     id, pointlevel, levels, concs
# end

vectorize(x) = [x]
vectorize(x::AbstractVector) = x
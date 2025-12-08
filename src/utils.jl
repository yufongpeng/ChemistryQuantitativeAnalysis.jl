"""
    cqaconvert(::Type{T}, x::S)
    cqaconvert(fn::Function, x::S)

Call constructor, `parse`, or `fn`. Extend this function if defining `T(x)` is not safe.

Return
* `parse(T, x)` if `T <: Number` and `S <: Union{AbstractString, AbstractChar}`.
* `fn(x)` if `fn` is a `Function`.
* `T(x)` otherwise.

Note that for any `x` of type `T`, `x == cqaconvert(T, string(x)))`.
"""
cqaconvert(::Type{T}, x) where T = T(x)::T
cqaconvert(::Type{T}, x::T) where T = x
cqaconvert(::Type{T}, x::Union{AbstractString, AbstractChar}) where {T <: Number} = parse(T, x)::T
cqaconvert(fn::T, x) where {T <: Function} = fn(x)

"""
    cqamap(::Type{T}, x::AbstractVector{S})
    cqamap(fn::Function, x::AbstractVector)

Convert `x` to `AbstractVector{T}`. If direct construction is not possible, i.e. neither `T <: S` nor `S <: T`, it applys `cqaconvert` on every elements.

For a function, `T` is inferred by `cqatype(fn)` first and then `cqatype(fn, v)` for the returned vector `v` to avoid abstract element type. 
"""
cqamap(::Type{T}, x::AbstractVector{T}) where T = x
cqamap(::Type{T}, x::AbstractVector{<: T}) where T = Vector{T}(x)
cqamap(::Type{T}, x::AbstractVector{S}) where {S, T <: S} = Vector{T}(x)
cqamap(::Type{T}, x::AbstractVector) where T = typedmap(e -> cqaconvert(T, e), T, x)
function cqamap(fn::T, x::AbstractVector) where {T <: Function}
    R = cqatype(fn)
    if isabstracttype(R)
        v = map(e -> cqaconvert(fn, e), x)
        cqamap(cqatype(fn, v), v)
    else
        typedmap(e -> cqaconvert(fn, e), R, x)
    end
end

"""
    cqatype(fn::Function)
    cqatype(fn::Function, x::AbstractVector)

Default returned type of `cqaconvert(fn, x)`. By default, the union of types of each element is used to avoid abstract element type. Extend this function if neccessary.
"""
cqatype(fn::T) where {T <: Function} = Any
cqatype(fn::T, v::AbstractVector) where {T <: Function} = Union{typeof.(v)...}

"""
    typedmap(::Type{T}, c...) 
    typedmap(f, ::Type{T}, c...)

Transform collection `c` by applying `f` to each element. For multiple collection arguments, apply `f` elementwise, and stop when any of them  
is exhausted. The element type of returned collection will be forced to `T`. In addtion, `typedmap(T, c...)` is equivalent to `typedmap(T, T, c...)`.
"""
typedmap(::Type{T}, iters...) where T = collect(T, Base.Generator(T, iters...))
typedmap(fn::F, ::Type{T}, iters...) where {F, T} = collect(T, Base.Generator(fn, iters...))

"""
    getformula(cal::ExternalCalibrator)
    getformula(type::Bool, zero::Bool)

Get a `FormulaTerm` based on `type` and `zero` or parameters from `cal`.

See `ExternalCalibrator` for detail description of `type` and `zero`.
"""
getformula(::Type{Proportional}) = @formula(y ~ 0 + x)
getformula(::Type{Linear}) = @formula(y ~ x)
getformula(::Type{QuadraticProportional}) = @formula(y ~ 0 + x + x ^ 2)
getformula(::Type{Quadratic}) = @formula(y ~ x + x ^ 2)
getformula(::Type{Logarithmic}) = @formula(y ~ log(x))
getformula(::Type{Exponential}) = @formula(y ~ exp(x))
getformula(::Type{Power}) = @formula(log(y) ~ log(x))

"""
    const_weight(x, y) = 1

Constant weight function
"""
const_weight(x, y) = 1.0

"""
    WFN::Dict{String, Function}

The default dictionary that maps names of weight function to the actual functions
"""
const WFN = Dict{String, Function}(
    "1"     => const_weight,
    "1/√x"  => (x, y) -> 1/sqrt(x),
    "1/x^0.5"  => (x, y) -> 1/sqrt(x),
    "1/x^(1/2)"  => (x, y) -> 1/sqrt(x),
    "1/x"   => (x, y) -> 1/x,
    "1/x²"  => (x, y) -> 1/x^2,
    "1/x^2"  => (x, y) -> 1/x^2,
    "1/√y"  => (x, y) -> 1/sqrt(y),
    "1/y^0.5"  => (x, y) -> 1/sqrt(y),
    "1/y^(1/2)"  => (x, y) -> 1/sqrt(y),
    "1/y"   => (x, y) -> 1/y,
    "1/y²"  => (x, y) -> 1/y^2,
    "1/y^2"  => (x, y) -> 1/y^2,
    "1/√(x+y)"  => (x, y) -> 1/sqrt(x+y),
    "1/(x+y)^0.5"  => (x, y) -> 1/sqrt(x+y),
    "1/(x+y)^(1/2)"  => (x, y) -> 1/sqrt(x+y),
    "1/(x+y)"   => (x, y) -> 1/(x+y),
    "1/(x+y)²"  => (x, y) -> 1/(x+y)^2,
    "1/(x+y)^2"  => (x, y) -> 1/(x+y)^2,
    "1/√|log(x)|"  => (x, y) -> 1/sqrt(abs(log(x))),
    "1/|log(x)|^0.5"  => (x, y) -> 1/sqrt(abs(log(x))),
    "1/|log(x)|^(1/2)"  => (x, y) -> 1/sqrt(abs(log(x))),
    "1/|log(x)|"   => (x, y) -> 1/abs(log(x)),
    "1/log(x)²"  => (x, y) -> 1/abs(log(x))^2,
    "1/log(x)^2"  => (x, y) -> 1/log(x)^2,
    "1/√|log(y)|"  => (x, y) -> 1/sqrt(abs(log(y))),
    "1/|log(y)|^0.5"  => (x, y) -> 1/sqrt(abs(log(y))),
    "1/|log(y)|^(1/2)"  => (x, y) -> 1/sqrt(abs(log(y))),
    "1/|log(y)|"   => (x, y) -> 1/abs(log(y)),
    "1/log(y)²"  => (x, y) -> 1/log(y)^2,
    "1/log(y)^2"  => (x, y) -> 1/log(y)^2
)

const WNM = Dict{String, Tuple{String, String}}(
    "1"     => ("1", "1"),
    "1/√x"  => ("1/√x", "1/x^(1/2)"),
    "1/x^0.5" => ("1/√x", "1/x^(1/2)"),
    "1/x^(1/2)" => ("1/√x", "1/x^(1/2)"),
    "1/x"   => ("1/x", "1/x"),
    "1/x²"  => ("1/x²", "1/x^2"),
    "1/x^2" => ("1/x²", "1/x^2"),
    "1/√y"  => ("1/√y", "1/y^(1/2)"),
    "1/y^0.5" => ("1/√y", "1/y^(1/2)"),
    "1/y^(1/2)" => ("1/√y", "1/y^(1/2)"),
    "1/y"   => ("1/y", "1/y"),
    "1/y²"  => ("1/y²", "1/y^2"),
    "1/y^2" => ("1/y²", "1/y^2"),
    "1/√(x+y)"  => ("1/√(x+y)", "1/(x+y)^(1/2)"),
    "1/(x+y)^0.5" => ("1/√(x+y)", "1/(x+y)^(1/2)"),
    "1/(x+y)^(1/2)" => ("1/√(x+y)", "1/(x+y)^(1/2)"),
    "1/(x+y)"   => ("1/(x+y)", "1/(x+y)"),
    "1/(x+y)²"  => ("1/(x+y)²", "1/(x+y)^2"),
    "1/(x+y)^2" => ("1/(x+y)²", "1/(x+y)^2"),
    "1/√|log(x)|"  => ("1/√|log(x)|", "1/|log(x)|^(1/2)"),
    "1/|log(x)|^0.5"  => ("1/√|log(x)|", "1/|log(x)|^(1/2)"),
    "1/|log(x)|^(1/2)"  => ("1/√|log(x)|", "1/|log(x)|^(1/2)"),
    "1/|log(x)|"   => ("1/|log(x)|", "1/|log(x)|"),
    "1/log(x)²"  => ("1/log(x)²", "1/log(x)^2"),
    "1/log(x)^2"  => ("1/log(x)²", "1/log(x)^2"),
    "1/√|log(y)|"  => ("1/|log(y)|", "1/|log(y)|^(1/2)"),
    "1/|log(y)|^0.5"  => ("1/√|log(y)|", "1/|log(y)|^(1/2)"),
    "1/|log(y)|^(1/2)"  => ("1/√|log(y)|", "1/|log(y)|^(1/2)"),
    "1/|log(y)|"   => ("1/|log(y)|", "1/|log(y)|"),
    "1/log(y)²"  => ("1/log(y)²", "1/log(y)^2"),
    "1/log(y)^2"  => ("1/log(y)²", "1/log(y)^2"),
)

getwts(wfn, x, y) = map(wfn, x, y)
# getwts(wfn::Nothing, x, y) = similar(y, 0)

"""
    stdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return calibratio standard of `analyte` based on `method`.
"""
function stdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = method.analytetable.std[aid]
    (iid > 0 ? method.analytetable.analyte[iid] : nothing)::Union{A, Nothing}
end
"""
    isdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return internal standard of `analyte` based on `method`.
"""
function isdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = method.analytetable.isd[aid]
    (iid > 0 ? method.analytetable.analyte[iid] : iid == 0 ? nothing : analyte)::Union{A, Nothing}
end
"""
    isstd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return if `analyte` is a internal standard based on `method`.
"""
function isstd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    !isnothing(findfirst(==(aid), method.analytetable.std))
end
"""
    isisd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return if `analyte` is a internal standard based on `method`.
"""
function isisd(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    method.analytetable.isd[aid] < 0
end

# getter method
table(dt::AbstractDataTable) = getfield(dt, :table)
tables(at::AnalysisTable) = getfield(at, :tables)
"""
    analyteobj(dt::AbstractDataTable{A}) -> Vector{A}
    analyteobj(at::AnalysisTable{A}) -> Vector{A}

Equivalent to `getfield(dt, :analyte)` or `getfield(at, :analyte)`.
"""
analyteobj(dt::AbstractDataTable) = getfield(dt, :analyte)
analyteobj(at::AnalysisTable) = getfield(at, :analyte)
"""
    sampleobj(dt::AbstractDataTable{A, S}) -> Vector{S}
    sampleobj(at::AnalysisTable{A, S, T}) -> Vector{S}

Equivalent to `getfield(dt, :sample)` or `getfield(at, :sample)`.
"""
sampleobj(dt::AbstractDataTable) = getfield(dt, :sample)
sampleobj(at::AnalysisTable) = getfield(at, :sample)
"""
    analytename(dt::AbstractDataTable) -> Vector{Symbol}
    analytename(at::AnalysisTable) -> Vector{Symbol}

Equivalent to `Symbol.(analyteobj(dt))` or `Symbol.(analyteobj(at))`.
"""
analytename(dt::AbstractDataTable) = Symbol.(analyteobj(dt))
analytename(at::AnalysisTable) = Symbol.(analyteobj(at))
"""
    samplename(dt::AbstractDataTable) -> Vector{Symbol}
    samplename(at::AnalysisTable) -> Vector{Symbol}

Equivalent to `Symbol.(sampleobj(dt))` or `Symbol.(sampleobj(at))`.
"""
samplename(dt::AbstractDataTable) = Symbol.(sampleobj(dt))
samplename(at::AnalysisTable) = Symbol.(sampleobj(at))
"""
    samplecol(dt::SampleDataTable) -> Symbol

Equivalent to `getfield(dt, :samplecol)`.
"""
samplecol(dt::SampleDataTable) = getfield(dt, :samplecol)
"""
    analytecol(dt::SampleDataTable) -> Symbol

Equivalent to `getfield(dt, :analytecol)`.
"""
analytecol(dt::AnalyteDataTable) = getfield(dt, :analytecol)
"""
    idcol(dt::AbstractDataTable) -> Symbol

Unified interface for identity column. Equivalent to `samplecol(dt::SampleDataTable)` or `analytecol(dt::SampleDataTable)`.
"""
idcol(dt::SampleDataTable) = samplecol(dt)
idcol(dt::AnalyteDataTable) = analytecol(dt)

"""
    findanalyte(dt::AbstractDataTable{A}, analyte::B) where {A, B <: A}
    findanalyte(dt::AbstractDataTable, analyte::Symbol)

Return the index of the first element of `analyteobj(dt)` or `analytename(dt)` for which the element equals to `analyte`.
"""
findanalyte(dt::AbstractDataTable{A}, analyte::B) where {A, B <: A} = findfirst(==(analyte), analyteobj(dt))
findanalyte(dt::AbstractDataTable, analyte::Symbol) = findfirst(==(analyte), analytename(dt))

"""
    findsample(dt::AbstractDataTable{A, S}, sample::T) where {A, S, T <: S}
    findsample(dt::AbstractDataTable, sample::Symbol)

Return the index of the first element of `sampleobj(dt)` or `samplename(dt)` for which the element equals to `sample`.
"""
findsample(dt::AbstractDataTable{A, S}, sample::T) where {A, S, T <: S} = findfirst(==(sample), sampleobj(dt))
findsample(dt::AbstractDataTable, sample::Symbol) = findfirst(==(sample), samplename(dt))

"""
    getanalyte(dt::AnalyteDataTable, id::Int)
    getanalyte(dt::SampleDataTable, id::Int)
    getanalyte(dt::AnalyteDataTable, analyte)
    getanalyte(dt::SampleDataTable, analyte)

Get data belonging to `analyte` or `analyteobj(dt)[id]` as a `Vector`. For `AnalyteDataTable`, a new vector is created; mutating this vector will not change the value in `dt`.
"""
getanalyte(dt::AnalyteDataTable{A, S, N}, id::Int) where {A, S, N} = [getproperty(dt, p)[id] for p in samplename(dt)]::Vector{N}
function getanalyte(dt::AnalyteDataTable{A, S, N}, analyte) where {A, S, N}
    id = findanalyte(dt, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    [getproperty(dt, p)[id] for p in samplename(dt)]::Vector{N}
end
getanalyte(dt::SampleDataTable{A, S, N}, id::Int) where {A, S, N} = getproperty(dt, analytename(dt)[id])::AbstractVector{N}
function getanalyte(dt::SampleDataTable{A, S, N}, analyte) where {A, S, N}
    id = findanalyte(dt, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    getproperty(dt, Symbol(analyteobj(dt)[id]))::AbstractVector{N}
end

"""
    getsample(dt::AnalyteDataTable, id::Int)
    getsample(dt::SampleDataTable, id::Int)
    getsample(dt::AnalyteDataTable, sample)
    getsample(dt::SampleDataTable, sample)

Get data belonging to `sample` or `sampleobj(dt)[id]` as a `Vector`. For `SampleDataTable`, a new vector is created; mutating this vector will not change the value in `dt`.
"""
getsample(dt::AnalyteDataTable{A, S, N}, id::Int) where {A, S, N} = getproperty(dt, Symbol(sampleobj(dt)[id]))::AbstractVector{N}
function getsample(dt::AnalyteDataTable{A, S, N}, sample) where {A, S, N}
    id = findsample(dt, sample)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    getproperty(dt, Symbol(sampleobj(dt)[id]))::AbstractVector{N}
end
getsample(dt::SampleDataTable{A, S, N}, id::Int) where {A, S, N} = [getproperty(dt, p)[id] for p in analytename(dt)]::Vector{N}
function getsample(dt::SampleDataTable{A, S, N}, sample) where {A, S, N}
    id = findsample(dt, sample)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    [getproperty(dt, p)[id] for p in analytename(dt)]::Vector{N}
end

findcalibrator(batch::Batch, analyte) = findcalibrator(batch.calibrator, analyte)
function findcalibrator(calibrator::AbstractVector{A}, analyte) where {A <: AbstractCalibrator}
    findfirst(x -> x.analyte == analyte, calibrator)
end
function findcalibrator(calibrator::AbstractVector{A}, analyte::Symbol) where {A <: AbstractCalibrator}
    findfirst(x -> Symbol(x.analyte) == analyte, calibrator)
end

getcalibrator(batch::Batch, analyte) = getcalibrator(batch.calibrator, analyte)
function getcalibrator(calibrator::AbstractVector{A}, analyte) where {A <: AbstractCalibrator}
    id = findcalibrator(calibrator, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte does not have a calibration curve."))
    calibrator[id]
end
getcalibrator(batch::Batch, id::Int) = getcalibrator(batch.calibrator, id)
getcalibrator(calibrator::AbstractVector{A}, id::Int) where {A <: AbstractCalibrator} = calibrator[id]

function critical_point(cal::ExternalCalibrator{A, N}) where {A, N}
    β = cal.model.model.pp.beta0::Vector{N}
    c, b, a = cal.zero ? (0, β...) : β
    -b / 2a
end

function parse_calibration_level_name(dt::AbstractDataTable, calid::Regex, order, ratio, df, f2c, parse_decimal)
    so = split(order, "")
    dilutionid = findfirst(==("D"), so)
    ratioid = findfirst(==("R"), so)
    levelid = findfirst(==("L"), so)
    cs = map(samplename(dt)) do s
            m = match(calid, string(s))
            isnothing(m) ? nothing : collect(String, m)
    end
    isnothing(levelid) && throw(ArgumentError("No valid level id."))
    clevel = map(s -> isnothing(s) ? nothing : parse(Int, s[levelid]), cs)
    id = findall(!isnothing, clevel)
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


function parse_calibration_level_name(dt::AbstractDataTable, calid, order, ratio, df, f2c, parse_decimal)
    clevel = replace(calid, missing => nothing)
    id = findall(!isnothing, clevel)
    pointlevel = convert(Vector{Int}, clevel[id])
    levels = unique(pointlevel)
    if !isnothing(ratio)
        concs = map(s -> s .* f2c, ratio)
    elseif !isnothing(df)
        concs = map(s -> f2c ./ s, df)
    else
        throw(ArgumentError("No valid ratios or dilution factors are obtained."))
    end
    id, pointlevel, levels, concs
end

# function parse_calibration_name(s, r, d, calid, levelid, ratioid, dilutionid, f2c, parse_decimal)
#     m = match(calid, string(s))
#     isnothing(m) && return missing
#     m = collect(String, m)
#     level = parse(Int, m[levelid])
#     conc = parse(Float64, parse_decimal(m[only(setdiff(eachindex(m), levelid))]))
#     (level, dilution ? f2c ./ conc : f2c .* conc)
# end

table_convert(::Type{D}, data::AbstractDataTable{A, S, N, T}) where {D, A, S, N, T <: D} = data
table_convert(::Type{D}, data::SampleDataTable) where D = 
    SampleDataTable(analyteobj(data), sampleobj(data), idcol(data), D(table(data)))
table_convert(::Type{D}, data::AnalyteDataTable) where D = 
    AnalyteDataTable(analyteobj(data), sampleobj(data), idcol(data), D(table(data)))
"""
    dynamic_range(cal::ExternalCalibrator)

Return dynamic range as a `Tuple` (lloq, uloq).
"""
dynamic_range(cal::ExternalCalibrator) = (lloq(cal), uloq(cal))
"""
    signal_range(cal::ExternalCalibrator)

Return theoretical signal of dynamic range as a `Tuple` (lloq, uloq).
"""
signal_range(cal::ExternalCalibrator) = (predict(cal.machine, Table(; x = collect(dynamic_range(cal))))..., )

"""
    lloq(cal::ExternalCalibrator)

Return lower limit of quantification.
"""
lloq(cal::ExternalCalibrator) = cal.table.x[findfirst(cal.table.include)]
"""
    uloq(cal::ExternalCalibrator)

Return upper limit of quantification.
"""
uloq(cal::ExternalCalibrator) = cal.table.x[findlast(cal.table.include)]
"""
    signal_lloq(cal::ExternalCalibrator)

Return theoretical signal of lower limit of quantification.
"""
signal_lloq(cal::ExternalCalibrator) = only(predict(cal.machine, Table(; x = [lloq(cal)])))
"""
    uloq(cal::ExternalCalibrator)

Return theoretical signal of upper limit of quantification.
"""
signal_uloq(cal::ExternalCalibrator) = only(predict(cal.machine, Table(; x = [uloq(cal)])))

# function blank(cal::ExternalCalibrator{A, N}) where {A, N}
#     β = cal.model.model.pp.beta0::Vector{N}
#     if cal.type
#         b = length(β) == 1 ? 0 : first(β)
#     else
#         length(β) == 2 && return 0
#         c, b, a = β
#         a < 0 && return c
#         max(0, b ^ 2 / 4a - b ^ 2 / 2a + c)
#     end
# end

blank(cal::ExternalCalibrator) = blank(cal.model)
blank(cal::CalibrationModel{Proportional}) = 0
blank(cal::CalibrationModel{Linear}) = 0
blank(cal::CalibrationModel{QuadraticProportional}) = 0
blank(cal::CalibrationModel{Quadratic}) = 0
blank(cal::CalibrationModel{Logarithmic}) = -Inf
blank(cal::CalibrationModel{Exponential}) = 0
blank(cal::CalibrationModel{Power}) = 0

function signal_lob(cal::ExternalCalibrator)
    zero_id = findall(==(0), cal.table.level)
    isempty(zero_id) ? blank(cal) + 1.645 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]]) : mean(cal.table.y[zero_id]) + 1.645 * std(cal.table.y[zero_id])
end

"""
    signal_loq(cal::ExternalCalibrator)

Theoretical limit of detection (LOD) in signal space.
"""
function signal_lod(cal::ExternalCalibrator) 
    zero_id = findall(==(0), cal.table.level)
    isempty(zero_id) ? blank(cal) + 3.3 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]]) : signal_lob(cal) + 1.645 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]])
end
"""
    loq(cal::ExternalCalibrator)

Theoretical limit of quantification (LOQ) in signal space.
"""
signal_loq(cal::ExternalCalibrator) = blank(cal) + 10 * std(cal.table.y[cal.table.x .== cal.table.x[findfirst(cal.table.include)]])

"""
    weight_repr(cal::ExternalCalibrator)
    weight_repr(weight::Number)

Return string representation of `cal.weight` or `weight`. 

"none" for 0, "1/√x" for -0.5, "1/x" for -1, "1/x²" for -2, 
"x^`\$weight`" for other positive `weight`, and "1/x^`\$(abs(weight))`" for other negative `weight`
"""
weight_repr(cal::InternalCalibrator) = "1"
weight_repr(cal::ExternalCalibrator) = weight_repr(cal.model)
weight_repr(model::CalibrationModel) = first(WNM[model.wnm])

"""
    formula_repr(cal::InternalCalibrator; digits = nothing, sigdigits = 4)
    formula_repr(cal::ExternalCalibrator; digits = nothing, sigdigits = 4)

Return string representation of formula of `cal` with specified `digits` and `sigdigits`. See `format_number`.
"""
formula_repr(cal::InternalCalibrator; digits = nothing, sigdigits = 4) = "y = $(format_number(1/cal.conc)); digits, sigdigits))x"
formula_repr(cal::ExternalCalibrator; digits = nothing, sigdigits = 4) = formula_repr(cal.model, cal.machine; digits, sigdigits)
function formula_repr(model::CalibrationModel, machine; digits = nothing, sigdigits = 4)
    β = machine.model.pp.beta0
    formula_repr(model, β; digits, sigdigits)
end
formula_repr(::CalibrationModel{Proportional}, β::AbstractVector; digits = nothing, sigdigits = 4) = "y = $(format_number(β[1]; digits, sigdigits))x"
function formula_repr(::CalibrationModel{Linear}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "x")
end
function formula_repr(::CalibrationModel{QuadraticProportional}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), "x", op[1], abs(format_number(β[2]; digits, sigdigits)), "x²")
end
function formula_repr(::CalibrationModel{Quadratic}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "x", op[2], abs(format_number(β[3]; digits, sigdigits)), "x²")
end
function formula_repr(::CalibrationModel{Logarithmic}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "log(x)")
end
function formula_repr(::CalibrationModel{Exponential}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("log(y) = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "x")
end
function formula_repr(::CalibrationModel{Power}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("log(y) = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "log(x)")
end
"""
    formula_repr_ascii(cal::AbstractCalibrator; digits = nothing, sigdigits = 4)

Return string representation of formula of `cal` for text file output.
"""
formula_repr_ascii(cal::AbstractCalibrator; digits = nothing, sigdigits = 4) = replace(formula_repr(cal; digits, sigdigits), "x²" => "x^2")
"""
    weight_repr_ascii(cal::AbstractCalibrator)
    weight_repr_ascii(weight::Number)

Return string representation of `cal.weight` or `weight` for text file output.
"""
weight_repr_ascii(cal::InternalCalibrator) = "1"
weight_repr_ascii(cal::ExternalCalibrator) = weight_repr_ascii(cal.model)
weight_repr_ascii(model::CalibrationModel) = last(WNM[model.wnm])

"""
    format_number(x; digits = nothing, sigdigits = 4)

Return string representation of number with specified `digits`.

If `digits` is `nothing`, the function uses `sigdigits` instead.
"""
format_number(x; digits = nothing, sigdigits = 4) = isnothing(digits) ? format_number2int(round(x; sigdigits)) : format_number2int(round(x; digits))
format_number2int(x) = 
    x == round(x) ? round(Int, x) : x

vectorize(x) = [x]
vectorize(x::AbstractVector) = x
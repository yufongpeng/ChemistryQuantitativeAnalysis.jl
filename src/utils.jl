"""
    cqaconvert(::Type{T}, x::S)
    cqaconvert(fn::Function, x::S)

Call constructor, `parse`, or `fn`. Extend this function if defining `T(x)` is not safe.

Return
* `parse(T, x)` if `T <: Number` and `S <: Union{AbstractString, AbstractChar}`.
* `fn(x)` if `fn` is a `Function`.
* `T(x)` otherwise.

Note that `string(cqaconvert(T, x))` should equal `x` for a valid string `x`, i.e. `string` and `x -> cqaconvert(T, x)` are inverse functions.
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
    getformula(cal::MultipleCalibration)
    getformula(type::Bool, zero::Bool)

Get a `FormulaTerm` based on `type` and `zero` or parameters from `cal`.

See `MultipleCalibration` for detail description of `type` and `zero`.
"""
getformula(cal::MultipleCalibration) = getformula(cal.type, cal.zero)
getformula(type::Bool, zero::Bool) = if type 
    zero ? @formula(y ~ 0 + x) : @formula(y ~ x)
else
    zero ? @formula(y ~ 0 + x + x ^ 2) : @formula(y ~ x + x ^ 2)
end

"""
    isdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}

Return internal standard of `analyte` based on `method`.
"""
function isdof(analyte::B, method::AnalysisMethod{A}) where {A, B <: A}
    aid = findfirst(==(analyte), method.analytetable.analyte)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = method.analytetable.isd[aid]
    (iid > 0 ? method.analytetable.analyte[iid] : nothing)::Union{A, Nothing}
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
    samplecol(dt::ColumnDataTable) -> Symbol

Equivalent to `getfield(dt, :samplecol)`.
"""
samplecol(dt::ColumnDataTable) = getfield(dt, :samplecol)
"""
    analytecol(dt::ColumnDataTable) -> Symbol

Equivalent to `getfield(dt, :analytecol)`.
"""
analytecol(dt::RowDataTable) = getfield(dt, :analytecol)

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
    getanalyte(dt::RowDataTable, id::Int)
    getanalyte(dt::ColumnDataTable, id::Int)
    getanalyte(dt::RowDataTable, analyte)
    getanalyte(dt::ColumnDataTable, analyte)

Get data belonging to `analyte` or `analyteobj(dt)[id]` as a `Vector`. For `RowDataTable`, a new vector is created; mutating this vector will not change the value in `dt`.
"""
getanalyte(dt::RowDataTable{A, S, N}, id::Int) where {A, S, N} = [getproperty(dt, p)[id] for p in samplename(dt)]::Vector{N}
function getanalyte(dt::RowDataTable{A, S, N}, analyte) where {A, S, N}
    id = findanalyte(dt, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    [getproperty(dt, p)[id] for p in samplename(dt)]::Vector{N}
end
getanalyte(dt::ColumnDataTable{A, S, N}, id::Int) where {A, S, N} = getproperty(dt, analytename(dt)[id])::AbstractVector{N}
function getanalyte(dt::ColumnDataTable{A, S, N}, analyte) where {A, S, N}
    id = findanalyte(dt, analyte)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    getproperty(dt, Symbol(analyteobj(dt)[id]))::AbstractVector{N}
end

"""
    getsample(dt::RowDataTable, id::Int)
    getsample(dt::ColumnDataTable, id::Int)
    getsample(dt::RowDataTable, sample)
    getsample(dt::ColumnDataTable, sample)

Get data belonging to `sample` or `sampleobj(dt)[id]` as a `Vector`. For `ColumnDataTable`, a new vector is created; mutating this vector will not change the value in `dt`.
"""
getsample(dt::RowDataTable{A, S, N}, id::Int) where {A, S, N} = getproperty(dt, Symbol(sampleobj(dt)[id]))::AbstractVector{N}
function getsample(dt::RowDataTable{A, S, N}, sample) where {A, S, N}
    id = findsample(dt, sample)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    getproperty(dt, Symbol(sampleobj(dt)[id]))::AbstractVector{N}
end
getsample(dt::ColumnDataTable{A, S, N}, id::Int) where {A, S, N} = [getproperty(dt, p)[id] for p in analytename(dt)]::Vector{N}
function getsample(dt::ColumnDataTable{A, S, N}, sample) where {A, S, N}
    id = findsample(dt, sample)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    [getproperty(dt, p)[id] for p in analytename(dt)]::Vector{N}
end

function critical_point(cal::MultipleCalibration{A, N}) where {A, N}
    β = cal.model.model.pp.beta0::Vector{N}
    c, b, a = cal.zero ? (0, β...) : β
    -b / 2a
end

"""
    dynamic_range(cal::MultipleCalibration)

Return dynamic range as a `Tuple` (lloq, uloq).
"""
dynamic_range(cal::MultipleCalibration) = (lloq(cal), uloq(cal))
"""
    signal_range(cal::MultipleCalibration)

Return theoretical signal of dynamic range as a `Tuple` (lloq, uloq).
"""
signal_range(cal::MultipleCalibration) = Tuple(predict(cal.model, Table(; x = collect(dynamic_range(cal)))))

"""
    lloq(cal::MultipleCalibration)

Return lower limit of quantification.
"""
lloq(cal::MultipleCalibration) = (cal.type || last(cal.model.model.pp.beta0) < 0) ? cal.table.x[findfirst(cal.table.include)] : max(cal.table.x[findfirst(cal.table.include)], critical_point(cal))
"""
    uloq(cal::MultipleCalibration)

Return upper limit of quantification.
"""
uloq(cal::MultipleCalibration) = (cal.type || last(cal.model.model.pp.beta0) > 0) ? cal.table.x[findlast(cal.table.include)] : min(cal.table.x[findlast(cal.table.include)], critical_point(cal))
"""
    signal_lloq(cal::MultipleCalibration)

Return theoretical signal of lower limit of quantification.
"""
signal_lloq(cal::MultipleCalibration) = only(predict(cal.model, Table(; x = [lloq(cal)])))
"""
    uloq(cal::MultipleCalibration)

Return theoretical signal of upper limit of quantification.
"""
signal_uloq(cal::MultipleCalibration) = only(predict(cal.model, Table(; x = [uloq(cal)])))

"""
    weight_repr(cal::MultipleCalibration)
    weight_repr(weight::Number)

Return string representation of `cal.weight` or `weight`. 

"none" for 0, "1/√x" for -0.5, "1/x" for -1, "1/x²" for -2, 
"x^`\$weight`" for other positive `weight`, and "1/x^`\$(abs(weight))`" for other negative `weight`
"""
function weight_repr(cal::MultipleCalibration)
    weight_repr(cal.weight)
end
weight_repr(weight::Number) = if weight == -0.5
    "1/√x"
elseif weight == -1
    "1/x"
elseif weight == -2
    "1/x²"
elseif weight == 0
    "none"
elseif weight > 0
    "x^$weight"
else
    "1/x^$(abs(weight))"
end
"""
    weight_value(weight)

Return value of weight from string representation. See "weight_repr".
"""
weight_value(weight) = if weight == "1/√x"
    -0.5
elseif weight == "1/x"
    -1
elseif weight == "1/x²"
    -2
else
    neg = match(r"1/x\^(\d*)", weight)
    isnothing(neg) ? parse(Float64, first(match(r"x\^(\d*)", weight))) : -parse(Float64, first(neg))
end

"""
    formula_repr(cal::SingleCalibration; digits = nothing, sigdigits = 4)
    formula_repr(cal::MultipleCalibration; digits = nothing, sigdigits = 4)

Return string representation of formula of `cal` with specified `digits` and `sigdigits`. See `format_number`.
"""
formula_repr(cal::SingleCalibration; digits = nothing, sigdigits = 4) = "y = $(format_number(1/first(getanalyte(cal.caltable.conctable, cal.isd)); digits, sigdigits))x"
function formula_repr(cal::MultipleCalibration; digits = nothing, sigdigits = 4)
    β = cal.model.model.pp.beta0
    cal.type && cal.zero && return "y = $(format_number(β[1]; digits, sigdigits))x"
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    if cal.type
        string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "x")
    elseif cal.zero
        string("y = ", format_number(β[1]; digits, sigdigits), "x", op[1], abs(format_number(β[2]; digits, sigdigits)), "x²")
    else
        string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "x", op[2], abs(format_number(β[3]; digits, sigdigits)), "x²")
    end
end

"""
    formula_repr_utf8(cal::AbstractCalibration)

Return string representation of formula of `cal` for text file output.
"""
formula_repr_utf8(cal::AbstractCalibration) = replace(formula_repr(cal), "x²" => "x^2")
"""
    weight_repr_utf8(cal::AbstractCalibration)
    weight_repr_utf8(weight::Number)

Return string representation of `cal.weight` or `weight` for text file output.
"""
weight_repr_utf8(cal::AbstractCalibration) = replace(weight_repr(cal), "x²" => "x^2", "√x" => "x^0.5")
weight_repr_utf8(weight::Number) = replace(weight_repr(weight), "x²" => "x^2", "√x" => "x^0.5")
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
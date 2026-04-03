abstract type CurveType end 
struct LinearOrigin <: CurveType end
const Proportional = LinearOrigin
struct Linear <: CurveType end
struct QuadraticOrigin <: CurveType end
@deprecate QuadraticProportional QuadraticOrigin
struct Quadratic <: CurveType end
struct Logarithmic <: CurveType end
struct Exponential <: CurveType end
struct Power <: CurveType end

@deprecate ComposedWeight CalibrationWeight
wdoc = """
    struct CalibrationWeight{X, Y, S, N} end 

    const ConstWeight = CalibrationWeight{Nothing, Nothing, WSum, 1}
    const RootXWeight = CalibrationWeight{Proportional, Nothing, WSum, 1//2}
    const RootYWeight = CalibrationWeight{Nothing, Proportional, WSum, 1//2}
    const RootXYWeight = CalibrationWeight{Proportional, Proportional, WSum, 1//2}
    const XWeight = CalibrationWeight{Proportional, Nothing, WSum, 1}
    const YWeight = CalibrationWeight{Nothing, Proportional, WSum, 1}
    const XYWeight = CalibrationWeight{Proportional, Proportional, WSum, 1}
    const SqXWeight = CalibrationWeight{Proportional, Nothing, WSum, 2}
    const SqYWeight = CalibrationWeight{Nothing, Proportional, WSum, 2}
    const SqXYWeight = CalibrationWeight{Proportional, Proportional, WSum, 2}

    const RootLogXWeight = CalibrationWeight{Logarithmic, Nothing, WSum, 1//2}
    const RootLogYWeight = CalibrationWeight{Nothing, Logarithmic, WSum, 1//2}
    const LogXWeight = CalibrationWeight{Logarithmic, Nothing, WSum, 1}
    const LogYWeight = CalibrationWeight{Nothing, Logarithmic, WSum, 1}
    const SqLogXWeight = CalibrationWeight{Logarithmic, Nothing, WSum, 2}
    const SqLogYWeight = CalibrationWeight{Nothing, Logarithmic, WSum, 2}

    const RootExpXWeight = CalibrationWeight{Exponential, Nothing, WSum, 1//2}
    const RootExpYWeight = CalibrationWeight{Nothing, Exponential, WSum, 1//2}
    const ExpXWeight = CalibrationWeight{Exponential, Nothing, WSum, 1}
    const ExpYWeight = CalibrationWeight{Nothing, Exponential, WSum, 1}
    const SqExpXWeight = CalibrationWeight{Exponential, Nothing, WSum, 2}
    const SqExpYWeight = CalibrationWeight{Nothing, Exponential, WSum, 2}

Weight object to generate weight function. 

* `X::Type{<: CurveType}`: operation for x. `Nothing` indicates no x weight.
* `Y::Type{<: CurveType}`: operation for y. `Nothing` indicates no y weight.
* `S::Type`: operation for combining x and y.
* `N::Real`: degree of final polynomial or rational function after combing x and y. 

`CurveType` is restructed to `Proportional`, `Exponential`, `Logarithmic`.

Weight values can be computed using `getweights`.
"""
@doc wdoc
struct CalibrationWeight{X, Y, S, N} end 

"""
    struct WSum end 

Operation for summing x and y.
"""
struct WSum end 

@doc wdoc
const ConstWeight = CalibrationWeight{Nothing, Nothing, WSum, 1}
@doc wdoc
const RootXWeight = CalibrationWeight{Proportional, Nothing, WSum, 1//2}
@doc wdoc
const RootYWeight = CalibrationWeight{Nothing, Proportional, WSum, 1//2}
@doc wdoc
const RootXYWeight = CalibrationWeight{Proportional, Proportional, WSum, 1//2}
@doc wdoc
const XWeight = CalibrationWeight{Proportional, Nothing, WSum, 1}
@doc wdoc
const YWeight = CalibrationWeight{Nothing, Proportional, WSum, 1}
@doc wdoc
const XYWeight = CalibrationWeight{Proportional, Proportional, WSum, 1}
@doc wdoc
const SqXWeight = CalibrationWeight{Proportional, Nothing, WSum, 2}
@doc wdoc
const SqYWeight = CalibrationWeight{Nothing, Proportional, WSum, 2}
@doc wdoc
const SqXYWeight = CalibrationWeight{Proportional, Proportional, WSum, 2}

@doc wdoc
const RootLogXWeight = CalibrationWeight{Logarithmic, Nothing, WSum, 1//2}
@doc wdoc
const RootLogYWeight = CalibrationWeight{Nothing, Logarithmic, WSum, 1//2}
@doc wdoc
const LogXWeight = CalibrationWeight{Logarithmic, Nothing, WSum, 1}
@doc wdoc
const LogYWeight = CalibrationWeight{Nothing, Logarithmic, WSum, 1}
@doc wdoc
const SqLogXWeight = CalibrationWeight{Logarithmic, Nothing, WSum, 2}
@doc wdoc
const SqLogYWeight = CalibrationWeight{Nothing, Logarithmic, WSum, 2}

@doc wdoc
const RootExpXWeight = CalibrationWeight{Exponential, Nothing, WSum, 1//2}
@doc wdoc
const RootExpYWeight = CalibrationWeight{Nothing, Exponential, WSum, 1//2}
@doc wdoc
const ExpXWeight = CalibrationWeight{Exponential, Nothing, WSum, 1}
@doc wdoc
const ExpYWeight = CalibrationWeight{Nothing, Exponential, WSum, 1}
@doc wdoc
const SqExpXWeight = CalibrationWeight{Exponential, Nothing, WSum, 2}
@doc wdoc
const SqExpYWeight = CalibrationWeight{Nothing, Exponential, WSum, 2}

"""
    const_weight(x, y) = 1

Constant weight function
"""
const_weight(x::S, y::T) where {S, T} = promote_type(S, T)(1.0)

"""
    getweights(::CalibrationWeight, x, y)

Get weight values from weight object.
"""
function getweights(::CalibrationWeight{X, Y, S, N}, x, y) where {X, Y, S, N}
    if X == Nothing && Y == Nothing 
        getweights(ConstWeight(), x, y)
    elseif Y == Nothing
        @. 1 / Ref(wtsop(X))(x) ^ N
    elseif X == Nothing 
        @. 1 / Ref(wtsop(Y))(y) ^ N
    else
        @. 1 / Ref(wtsop(S))(Ref(wtsop(X))(x), Ref(wtsop(Y))(y)) ^ N
    end
end
getweights(::ConstWeight, x, y) = @. const_weight(x, y)
getweights(::RootXWeight, x, y) = @. 1 / sqrt(x)
getweights(::RootYWeight, x, y) = @. 1 / sqrt(y)
getweights(::RootXYWeight, x, y) = @. 1 / sqrt(x + y)
function getweights(::RootLogXWeight, x, y) 
    @. 1 / sqrt(max(log(x), eps(x)))
end
function getweights(::RootLogYWeight, x, y) 
    @. 1 / sqrt(max(log(y), eps(y)))
end
getweights(::RootExpXWeight, x, y) = @. 1 / sqrt(exp(x))
getweights(::RootExpYWeight, x, y) = @. 1 / sqrt(exp(y))
getweights(::XWeight, x, y) = @. 1 / x
getweights(::YWeight, x, y) = @. 1 / y
getweights(::XYWeight, x, y) = @. 1 / (x + y)
function getweights(::LogXWeight, x, y) 
    @. 1 / max(log(x), eps(x))
end
function getweights(::LogYWeight, x, y) 
    @. 1 / max(log(y), eps(y))
end
getweights(::ExpXWeight, x, y) = @. 1 / exp(x)
getweights(::ExpYWeight, x, y) = @. 1 / exp(y)
getweights(::SqXWeight, x, y) = @. 1 / x ^ 2
getweights(::SqYWeight, x, y) = @. 1 / y ^ 2
getweights(::SqXYWeight, x, y) = @. 1 / (x + y) ^ 2
function getweights(::SqLogXWeight, x, y) 
    @. 1 / max(log(x), eps(x)) ^ 2
end
function getweights(::SqLogYWeight, x, y) 
    @. 1 / max(log(y), eps(y)) ^ 2
end
getweights(::SqExpXWeight, x, y) = @. 1 / exp(2 * x)
getweights(::SqExpYWeight, x, y) = @. 1 / exp(2 * y)

"""
    wtsop(::WSum)
    wtsop(::Proportional)
    wtsop(::Logarithmic)
    wtsop(::Exponential)

Real operation function. 
"""
wtsop(::WSum) = +
wtsop(::Proportional) = identifty
wtsop(::Logarithmic) = log
wtsop(::Exponential) = exp


# """
#     WFN::Dict{String, Function}

# The default dictionary that maps names of weight function to the actual functions
# """
# const WEIGHT_FUNCTION = Dict{String, Function}(
#     "1"     => const_weight,
#     "1/√x"  => (x, y) -> 1/sqrt(x),
#     "1/x^0.5"  => (x, y) -> 1/sqrt(x),
#     "1/x^(1/2)"  => (x, y) -> 1/sqrt(x),
#     "1/x"   => (x, y) -> 1/x,
#     "1/x²"  => (x, y) -> 1/x^2,
#     "1/x^2"  => (x, y) -> 1/x^2,
#     "1/√y"  => (x, y) -> 1/sqrt(y),
#     "1/y^0.5"  => (x, y) -> 1/sqrt(y),
#     "1/y^(1/2)"  => (x, y) -> 1/sqrt(y),
#     "1/y"   => (x, y) -> 1/y,
#     "1/y²"  => (x, y) -> 1/y^2,
#     "1/y^2"  => (x, y) -> 1/y^2,
#     "1/√(x+y)"  => (x, y) -> 1/sqrt(x+y),
#     "1/(x+y)^0.5"  => (x, y) -> 1/sqrt(x+y),
#     "1/(x+y)^(1/2)"  => (x, y) -> 1/sqrt(x+y),
#     "1/(x+y)"   => (x, y) -> 1/(x+y),
#     "1/(x+y)²"  => (x, y) -> 1/(x+y)^2,
#     "1/(x+y)^2"  => (x, y) -> 1/(x+y)^2,
#     "1/√|log(x)|"  => (x, y) -> 1/sqrt(abs(log(x))),
#     "1/|log(x)|^0.5"  => (x, y) -> 1/sqrt(abs(log(x))),
#     "1/|log(x)|^(1/2)"  => (x, y) -> 1/sqrt(abs(log(x))),
#     "1/|log(x)|"   => (x, y) -> 1/abs(log(x)),
#     "1/log(x)²"  => (x, y) -> 1/abs(log(x))^2,
#     "1/log(x)^2"  => (x, y) -> 1/log(x)^2,
#     "1/√|log(y)|"  => (x, y) -> 1/sqrt(abs(log(y))),
#     "1/|log(y)|^0.5"  => (x, y) -> 1/sqrt(abs(log(y))),
#     "1/|log(y)|^(1/2)"  => (x, y) -> 1/sqrt(abs(log(y))),
#     "1/|log(y)|"   => (x, y) -> 1/abs(log(y)),
#     "1/log(y)²"  => (x, y) -> 1/log(y)^2,
#     "1/log(y)^2"  => (x, y) -> 1/log(y)^2,
#     "1/√exp(x)"  => (x, y) -> 1/exp(x/2),
#     "1/exp(x/2)"  => (x, y) -> 1/exp(x/2),
#     "1/√e^x"  => (x, y) -> 1/exp(x/2),
#     "1/√eˣ"  => (x, y) -> 1/exp(x/2),
#     "1/e^(x/2)"  => (x, y) -> 1/exp(x/2),
#     "1/exp(x)"   => (x, y) -> 1/exp(x),
#     "1/e^x"   => (x, y) -> 1/exp(x),
#     "1/eˣ"   => (x, y) -> 1/exp(x),
#     "1/exp(2x)"  => (x, y) -> 1/exp(2x),
#     "1/e^(2x)"  => (x, y) -> 1/exp(2x),
#     "1/e^2x"  => (x, y) -> 1/exp(2x),
#     "1/e²ˣ"  => (x, y) -> 1/exp(2x),
#     "1/√exp(y)"  => (x, y) -> 1/exp(y/2),
#     "1/exp(y/2)"  => (x, y) -> 1/exp(y/2),
#     "1/√e^y"  => (x, y) -> 1/exp(y/2),
#     "1/√eʸ"  => (x, y) -> 1/exp(y/2),
#     "1/e^(y/2)"  => (x, y) -> 1/exp(y/2),
#     "1/exp(y)"   => (x, y) -> 1/exp(y),
#     "1/e^y"   => (x, y) -> 1/exp(y),
#     "1/eʸ"   => (x, y) -> 1/exp(y),
#     "1/exp(2y)"  => (x, y) -> 1/exp(2y),
#     "1/e^(2y)"  => (x, y) -> 1/exp(2y),
#     "1/e^2y"  => (x, y) -> 1/exp(2y),
#     "1/e²ʸ"  => (x, y) -> 1/exp(2y)

# )

# const WEIGHT_NAME = Dict{String, Tuple{String, String}}(
#     "1"     => ("1", "1"),
#     "1/√x"  => ("1/√x", "1/x^(1/2)"),
#     "1/x^0.5" => ("1/√x", "1/x^(1/2)"),
#     "1/x^(1/2)" => ("1/√x", "1/x^(1/2)"),
#     "1/x"   => ("1/x", "1/x"),
#     "1/x²"  => ("1/x²", "1/x^2"),
#     "1/x^2" => ("1/x²", "1/x^2"),
#     "1/√y"  => ("1/√y", "1/y^(1/2)"),
#     "1/y^0.5" => ("1/√y", "1/y^(1/2)"),
#     "1/y^(1/2)" => ("1/√y", "1/y^(1/2)"),
#     "1/y"   => ("1/y", "1/y"),
#     "1/y²"  => ("1/y²", "1/y^2"),
#     "1/y^2" => ("1/y²", "1/y^2"),
#     "1/√(x+y)"  => ("1/√(x+y)", "1/(x+y)^(1/2)"),
#     "1/(x+y)^0.5" => ("1/√(x+y)", "1/(x+y)^(1/2)"),
#     "1/(x+y)^(1/2)" => ("1/√(x+y)", "1/(x+y)^(1/2)"),
#     "1/(x+y)"   => ("1/(x+y)", "1/(x+y)"),
#     "1/(x+y)²"  => ("1/(x+y)²", "1/(x+y)^2"),
#     "1/(x+y)^2" => ("1/(x+y)²", "1/(x+y)^2"),
#     "1/√|log(x)|"  => ("1/√|log(x)|", "1/|log(x)|^(1/2)"),
#     "1/|log(x)|^0.5"  => ("1/√|log(x)|", "1/|log(x)|^(1/2)"),
#     "1/|log(x)|^(1/2)"  => ("1/√|log(x)|", "1/|log(x)|^(1/2)"),
#     "1/|log(x)|"   => ("1/|log(x)|", "1/|log(x)|"),
#     "1/log(x)²"  => ("1/log(x)²", "1/log(x)^2"),
#     "1/log(x)^2"  => ("1/log(x)²", "1/log(x)^2"),
#     "1/√|log(y)|"  => ("1/|log(y)|", "1/|log(y)|^(1/2)"),
#     "1/|log(y)|^0.5"  => ("1/√|log(y)|", "1/|log(y)|^(1/2)"),
#     "1/|log(y)|^(1/2)"  => ("1/√|log(y)|", "1/|log(y)|^(1/2)"),
#     "1/|log(y)|"   => ("1/|log(y)|", "1/|log(y)|"),
#     "1/log(y)²"  => ("1/log(y)²", "1/log(y)^2"),
#     "1/log(y)^2"  => ("1/log(y)²", "1/log(y)^2"),
#     "1/√exp(x)"  => ("1/√eˣ", "1/e^(x/2)"),
#     "1/exp(x/2)"  => ("1/√eˣ", "1/e^(x/2)"),
#     "1/√e^x"  => ("1/√eˣ", "1/e^(x/2)"),
#     "1/√eˣ"  => ("1/√eˣ", "1/e^(x/2)"),
#     "1/e^(x/2)"  => ("1/√eˣ", "1/e^(x/2)"),
#     "1/exp(x)"   => ("1/eˣ", "1/e^x"),
#     "1/e^x"   => ("1/eˣ", "1/e^x"),
#     "1/eˣ"   => ("1/eˣ", "1/e^x"),
#     "1/exp(2x)"  => ("1/e²ˣ", "1/e^2x"),
#     "1/e^(2x)"  => ("1/e²ˣ", "1/e^2x"),
#     "1/e^2x"  => ("1/e²ˣ", "1/e^2x"),
#     "1/e²ˣ"  => ("1/e²ˣ", "1/e^2x"),
#     "1/√exp(y)"  => ("1/√eʸ", "1/e^(y/2)"),
#     "1/exp(y/2)"  => ("1/√eʸ", "1/e^(y/2)"),
#     "1/√e^y"  => ("1/√eʸ", "1/e^(y/2)"),
#     "1/√eʸ"  => ("1/√eʸ", "1/e^(y/2)"),
#     "1/e^(y/2)"  => ("1/√eʸ", "1/e^(y/2)"),
#     "1/exp(y)"   => ("1/eʸ", "1/e^y"),
#     "1/e^y"   => ("1/eʸ", "1/e^y"),
#     "1/eʸ"   => ("1/eʸ", "1/e^y"),
#     "1/exp(2y)"  => ("1/e²ʸ", "1/e^2y"),
#     "1/e^(2y)"  => ("1/e²ʸ", "1/e^2y"),
#     "1/e^2y"  => ("1/e²ʸ", "1/e^2y"),
#     "1/e²ʸ"  => ("1/e²ʸ", "1/e^2y")

# )
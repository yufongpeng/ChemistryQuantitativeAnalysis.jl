abstract type CurveType end 
struct Proportional <: CurveType end
struct Linear <: CurveType end
struct QuadraticProportional <: CurveType end
struct Quadratic <: CurveType end
struct Logarithmic <: CurveType end
struct Exponential <: CurveType end
struct Power <: CurveType end

wdoc = """
    struct ComposedWeight{X, Y, S, N} end 

    const ConstWeight = ComposedWeight{Nothing, Nothing, WSum, 1}
    const RootXWeight = ComposedWeight{Proportional, Nothing, WSum, 1//2}
    const RootYWeight = ComposedWeight{Nothing, Proportional, WSum, 1//2}
    const RootXYWeight = ComposedWeight{Proportional, Proportional, WSum, 1//2}
    const XWeight = ComposedWeight{Proportional, Nothing, WSum, 1}
    const YWeight = ComposedWeight{Nothing, Proportional, WSum, 1}
    const XYWeight = ComposedWeight{Proportional, Proportional, WSum, 1}
    const SqXWeight = ComposedWeight{Proportional, Nothing, WSum, 2}
    const SqYWeight = ComposedWeight{Nothing, Proportional, WSum, 2}
    const SqXYWeight = ComposedWeight{Proportional, Proportional, WSum, 2}

    const RootLogXWeight = ComposedWeight{Logarithmic, Nothing, WSum, 1//2}
    const RootLogYWeight = ComposedWeight{Nothing, Logarithmic, WSum, 1//2}
    const LogXWeight = ComposedWeight{Logarithmic, Nothing, WSum, 1}
    const LogYWeight = ComposedWeight{Nothing, Logarithmic, WSum, 1}
    const SqLogXWeight = ComposedWeight{Logarithmic, Nothing, WSum, 2}
    const SqLogYWeight = ComposedWeight{Nothing, Logarithmic, WSum, 2}

    const RootExpXWeight = ComposedWeight{Exponential, Nothing, WSum, 1//2}
    const RootExpYWeight = ComposedWeight{Nothing, Exponential, WSum, 1//2}
    const ExpXWeight = ComposedWeight{Exponential, Nothing, WSum, 1}
    const ExpYWeight = ComposedWeight{Nothing, Exponential, WSum, 1}
    const SqExpXWeight = ComposedWeight{Exponential, Nothing, WSum, 2}
    const SqExpYWeight = ComposedWeight{Nothing, Exponential, WSum, 2}

Weight object to generate weight function. 

* `X::Type{<: CurveType}`: operation for x. `Nothing` indicates no x weight.
* `Y::Type{<: CurveType}`: operation for y. `Nothing` indicates no y weight.
* `S::Type`: operation for combining x and y.
* `N::Real`: degree of final polynomial or rational function after combing x and y. 

`CurveType` is restructed to `Proportional`, `Exponential`, `Logarithmic`.

Real function can be constructed using `getwts`.
"""
@doc wdoc
struct ComposedWeight{X, Y, S, N} end 

"""
    struct WSum end 

Operation for summing x and y.
"""
struct WSum end 

@doc wdoc
const ConstWeight = ComposedWeight{Nothing, Nothing, WSum, 1}
@doc wdoc
const RootXWeight = ComposedWeight{Proportional, Nothing, WSum, 1//2}
@doc wdoc
const RootYWeight = ComposedWeight{Nothing, Proportional, WSum, 1//2}
@doc wdoc
const RootXYWeight = ComposedWeight{Proportional, Proportional, WSum, 1//2}
@doc wdoc
const XWeight = ComposedWeight{Proportional, Nothing, WSum, 1}
@doc wdoc
const YWeight = ComposedWeight{Nothing, Proportional, WSum, 1}
@doc wdoc
const XYWeight = ComposedWeight{Proportional, Proportional, WSum, 1}
@doc wdoc
const SqXWeight = ComposedWeight{Proportional, Nothing, WSum, 2}
@doc wdoc
const SqYWeight = ComposedWeight{Nothing, Proportional, WSum, 2}
@doc wdoc
const SqXYWeight = ComposedWeight{Proportional, Proportional, WSum, 2}

@doc wdoc
const RootLogXWeight = ComposedWeight{Logarithmic, Nothing, WSum, 1//2}
@doc wdoc
const RootLogYWeight = ComposedWeight{Nothing, Logarithmic, WSum, 1//2}
@doc wdoc
const LogXWeight = ComposedWeight{Logarithmic, Nothing, WSum, 1}
@doc wdoc
const LogYWeight = ComposedWeight{Nothing, Logarithmic, WSum, 1}
@doc wdoc
const SqLogXWeight = ComposedWeight{Logarithmic, Nothing, WSum, 2}
@doc wdoc
const SqLogYWeight = ComposedWeight{Nothing, Logarithmic, WSum, 2}

@doc wdoc
const RootExpXWeight = ComposedWeight{Exponential, Nothing, WSum, 1//2}
@doc wdoc
const RootExpYWeight = ComposedWeight{Nothing, Exponential, WSum, 1//2}
@doc wdoc
const ExpXWeight = ComposedWeight{Exponential, Nothing, WSum, 1}
@doc wdoc
const ExpYWeight = ComposedWeight{Nothing, Exponential, WSum, 1}
@doc wdoc
const SqExpXWeight = ComposedWeight{Exponential, Nothing, WSum, 2}
@doc wdoc
const SqExpYWeight = ComposedWeight{Nothing, Exponential, WSum, 2}

"""
    const_weight(x, y) = 1

Constant weight function
"""
const_weight(x::S, y::T) where {S, T} = promote_type(S, T)(1.0)

"""
    getwts(::ComposedWeight, x, y)

Get weight values from weight object.
"""
function getwts(::ComposedWeight{X, Y, S, N}, x, y) where {X, Y, S, N}
    if X == Nothing && Y == Nothing 
        getwts(ConstWeight(), x, y)
    elseif Y == Nothing
        @. 1 / Ref(wtsop(X))(x) ^ N
    elseif X == Nothing 
        @. 1 / Ref(wtsop(Y))(y) ^ N
    else
        @. 1 / Ref(wtsop(S))(Ref(wtsop(X))(x), Ref(wtsop(Y))(y)) ^ N
    end
end
getwts(::ConstWeight, x, y) = @. const_weight(x, y)
getwts(::RootXWeight, x, y) = @. 1 / sqrt(x)
getwts(::RootYWeight, x, y) = @. 1 / sqrt(y)
getwts(::RootXYWeight, x, y) = @. 1 / sqrt(x + y)
function getwts(::RootLogXWeight, x, y) 
    @. 1 / sqrt(max(log(x), eps(x)))
end
function getwts(::RootLogYWeight, x, y) 
    @. 1 / sqrt(max(log(y), eps(y)))
end
getwts(::RootExpXWeight, x, y) = @. 1 / sqrt(exp(x))
getwts(::RootExpYWeight, x, y) = @. 1 / sqrt(exp(y))
getwts(::XWeight, x, y) = @. 1 / x
getwts(::YWeight, x, y) = @. 1 / y
getwts(::XYWeight, x, y) = @. 1 / (x + y)
function getwts(::LogXWeight, x, y) 
    @. 1 / max(log(x), eps(x))
end
function getwts(::LogYWeight, x, y) 
    @. 1 / max(log(y), eps(y))
end
getwts(::ExpXWeight, x, y) = @. 1 / exp(x)
getwts(::ExpYWeight, x, y) = @. 1 / exp(y)
getwts(::SqXWeight, x, y) = @. 1 / x ^ 2
getwts(::SqYWeight, x, y) = @. 1 / y ^ 2
getwts(::SqXYWeight, x, y) = @. 1 / (x + y) ^ 2
function getwts(::SqLogXWeight, x, y) 
    @. 1 / max(log(x), eps(x)) ^ 2
end
function getwts(::SqLogYWeight, x, y) 
    @. 1 / max(log(y), eps(y)) ^ 2
end
getwts(::SqExpXWeight, x, y) = @. 1 / exp(2 * x)
getwts(::SqExpYWeight, x, y) = @. 1 / exp(2 * y)

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
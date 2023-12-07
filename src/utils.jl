"""
    get_formula(cal::MultipleCalibration)
    get_formula(type::Bool, zero::Bool)

Get a `FormulaTerm` based on `type` and `zero` or parameters from `cal`.

See `MultipleCalibration` for detail description of `type` and `zero`.
"""
get_formula(cal::MultipleCalibration) = get_formula(cal.type, cal.zero)
get_formula(type::Bool, zero::Bool) = if type 
    zero ? @formula(y ~ 0 + x) : @formula(y ~ x)
else
    zero ? @formula(y ~ 0 + x + x ^ 2) : @formula(y ~ x + x ^ 2)
end

"""
    isd_of(method::MethodTable{A}, analyte::B) where {A, B <: A}

Return internal standard of `analyte` based on `method`.
"""
function isd_of(method::MethodTable{A}, analyte::B) where {A, B <: A}
    aid = findfirst(==(analyte), method.analyte_map.analytes)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    iid = method.analyte_map.isd[aid]
    iid > 0 ? method.analyte_map.analytes[iid] : nothing
end
"""
    isisd(method::MethodTable{A}, analyte::B) where {A, B <: A}

Return if `analyte` is a internal standard based on `method`.
"""
function isisd(method::MethodTable{A}, analyte::B) where {A, B <: A}
    aid = findfirst(==(analyte), method.analyte_map.analytes)
    isnothing(aid) && throw(ArgumentError("Analyte $analyte is not in the method"))
    method.analyte_map.isd[aid] < 0
end
"""
    find_analyte(dt::RowDataTable{A}, analyte::B) where {A, B <: A}
    find_analyte(dt::ColumnDataTable{A}, analyte::B) where {A, B <: A}

Return the index of the first element of `dt.analytes` for which the element equals to `analyte`.
"""
find_analyte(dt::RowDataTable{A}, analyte::B) where {A, B <: A} = findfirst(==(analyte), dt.analytes)
find_analyte(dt::ColumnDataTable{A}, analyte::B) where {A, B <: A} = findfirst(==(analyte), dt.analytes)

"""
    find_sample(dt::RowDataTable, sample::Symbol)
    find_sample(dt::ColumnDataTable, sample::Symbol)

Return the index of the first element of `dt.samples` for which the element equals to `sample`.
"""
find_sample(dt::RowDataTable, sample::Symbol) = findfirst(==(sample), dt.samples)
find_sample(dt::ColumnDataTable, sample::Symbol) = findfirst(==(sample), dt.samples)

"""
    get_analyte(dt::RowDataTable, id::Int)
    get_analyte(dt::ColumnDataTable, id::Int)
    get_analyte(dt::RowDataTable{A}, analyte::B) where {A, B <: A}
    get_analyte(dt::ColumnDataTable{A}, analyte::B) where {A, B <: A}

Get data belonging to `analyte` or `dt.analytes[id]` as a `Vector`.
"""
get_analyte(dt::RowDataTable, id::Int) = collect(getproperties(dt.table[id], Tuple(dt.sample_name)))
function get_analyte(dt::RowDataTable{A}, analyte::B) where {A, B <: A}
    id = findfirst(==(analyte), dt.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    [getproperty(dt.table, p)[id] for p in dt.sample_name]
end
get_analyte(dt::ColumnDataTable, id::Int) = getproperty(dt.table, dt.analyte_name[id])
function get_analyte(dt::ColumnDataTable{A}, analyte::B) where {A, B <: A}
    id = findfirst(==(analyte), dt.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    getproperty(dt.table, dt.analyte_name[id])
end

"""
    get_sample(dt::RowDataTable, id::Int)
    get_sample(dt::ColumnDataTable, id::Int)
    get_sample(dt::RowDataTable, sample::Symbol)
    get_sample(dt::ColumnDataTable, sample::Symbol)

Get data belonging to `sample` or `dt.samples[id]` as a `Vector`.
"""
get_sample(dt::RowDataTable, id::Int) = collect(getproperties(dt.table[id], Tuple(dt.sample_name)))
function get_sample(dt::RowDataTable, sample::Symbol)
    id = findfirst(==(sample), dt.samples)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    collect(getproperties(dt.table[id], Tuple(dt.sample_name)))
end
get_sample(dt::ColumnDataTable, id::Int) = getproperty(dt.table, dt.analyte_name[id])
function get_sample(dt::ColumnDataTable, sample::Symbol)
    id = findfirst(==(sample), dt.samples)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    getproperty(dt.table, dt.analyte_name[id])
end

function critical_point(cal::MultipleCalibration)
    β = cal.model.model.pp.beta0
    c, b, a = cal.zero ? (0, β...) : β
    -b / 2a
end

"""
    cal_range(cal::MultipleCalibration)

Return calibration range as a `Tuple` (lloq, uloq).
"""
cal_range(cal::MultipleCalibration) = (lloq(cal), uloq(cal))
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
formula_repr(cal::SingleCalibration; digits = nothing, sigdigits = 4) = "y = $(format_number(1/first(get_analyte(cal.caltable.conctable, cal.isd)); digits, sigdigits))x"
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

default_analyte_fn(x) = x
default_analyte_fn(x::AbstractString) = String(x)

vectorize(x) = [x]
vectorize(x::AbstractVector) = x
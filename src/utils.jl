"""
    get_formula(cal::MultipleCalibration)
    get_formula(type::Bool, zero::Bool)

Get a `FormulaTerm` based on `type` and `zero` or parameters from `cal`.

See `MultipleCalibration` for the meaning of `type` and `zero`.
"""
get_formula(cal::MultipleCalibration) = get_formula(cal.type, cal.zero)
get_formula(type::Bool, zero::Bool) = if type 
    zero ? @formula(y ~ 0 + x) : @formula(y ~ x)
else
    zero ? @formula(y ~ 0 + x + x ^ 2) : @formula(y ~ x + x ^ 2)
end

"""
    find_analyte(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A}
    find_analyte(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A}
    find_analyte(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A}

Return the index the first element of `tbl.analytes` for which the element equals to `analyte`.
"""
find_analyte(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A} = findfirst(==(analyte), tbl.analytes)
find_analyte(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A} = findfirst(==(analyte), tbl.analytes)
find_analyte(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A} = findfirst(==(analyte), tbl.analytes)

"""
    find_sample(tbl::RowAnalysisTable, sample::Symbol)
    find_sample(tbl::ColumnAnalysisTable, sample::Symbol)
    find_sample(tbl::SampleWrapperTable, sample::Symbol)

Return the index the first element of `tbl.samples` for which the element equals to `sample`.
"""
find_sample(tbl::RowAnalysisTable, sample::Symbol) = findfirst(==(sample), tbl.samples)
find_sample(tbl::ColumnAnalysisTable, sample::Symbol) = findfirst(==(sample), tbl.samples)
find_sample(tbl::SampleWrapperTable, sample::Symbol) = findfirst(==(sample), tbl.samples)

"""
    get_analyte(tbl::RowAnalysisTable, id::Int)
    get_analyte(tbl::ColumnAnalysisTable, id::Int)
    get_analyte(tbl::SampleWrapperTable, id::Int)
    get_analyte(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A}
    get_analyte(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A}
    get_analyte(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A}

Get data belonging to `analyte` or `tbl.analytes[id]` as a `Vector`.
"""
get_analyte(tbl::RowAnalysisTable, id::Int) = collect(getproperties(tbl.table[id], Tuple(tbl.sample_name)))
function get_analyte(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A}
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    collect(getproperties(tbl.table[id], Tuple(tbl.sample_name)))
end
get_analyte(tbl::ColumnAnalysisTable, id::Int) = getproperty(tbl.table, tbl.analyte_name[id])
function get_analyte(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A}
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    getproperty(tbl.table, tbl.analyte_name[id])
end
get_analyte(tbl::SampleWrapperTable, id::Int) = get_analyte(tbl.analysistable, id)
get_analyte(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A} = get_analyte(tbl.analysistable, analyte)

"""
    get_sample(tbl::RowAnalysisTable, id::Int)
    get_sample(tbl::ColumnAnalysisTable, id::Int)
    get_sample(tbl::SampleWrapperTable, id::Int)
    get_sample(tbl::RowAnalysisTable, sample::Symbol)
    get_sample(tbl::ColumnAnalysisTable, sample::Symbol)
    get_sample(tbl::SampleWrapperTable, sample::Symbol)

Get data belonging `to sample` or `tbl.samples[id]` as a `Vector`.
"""
get_sample(tbl::RowAnalysisTable, id::Int) = collect(getproperties(tbl.table[id], Tuple(tbl.sample_name)))
function get_sample(tbl::RowAnalysisTable, sample::Symbol)
    id = findfirst(==(sample), tbl.samples)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    collect(getproperties(tbl.table[id], Tuple(tbl.sample_name)))
end
get_sample(tbl::ColumnAnalysisTable, id::Int) = getproperty(tbl.table, tbl.analyte_name[id])
function get_sample(tbl::ColumnAnalysisTable, sample::Symbol)
    id = findfirst(==(sample), tbl.samples)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    getproperty(tbl.table, tbl.analyte_name[id])
end
get_sample(tbl::SampleWrapperTable, id::Int) = get_sample(tbl.analysistable, id)
get_sample(tbl::SampleWrapperTable, sample::Symbol) = get_sample(tbl.analysistable, sample)

"""
    find_isd(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A}
    find_isd(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A}
    find_isd(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A}

Return the index the first element of `tbl.analytes` for which the element is internal standard of `analyte`.
"""
function find_isd(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A}
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    getproperty(tbl.table, tbl.isd_map)[id]
end

function find_isd(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A}
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    tbl.isd_map[id]
end
find_isd(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A} = find_isd(tbl.analysistable, analyte)

"""
    get_isd(tbl::RowAnalysisTable, id::Int)
    get_isd(tbl::ColumnAnalysisTable, id::Int)
    get_isd(tbl::SampleWrapperTable, id::Int)
    get_isd(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A}
    get_isd(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A}
    get_isd(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A}

Get data belonging to internal standard of `analyte` or `tbl.analytes[id]` as a `Vector`. To be noticed, `id` is the index of internal standard itself.
"""
get_isd(tbl::RowAnalysisTable, id::Int) = id < 1 ? ones(Float64, length(tbl.samples)) : get_analyte(tbl, id)
function get_isd(tbl::RowAnalysisTable{A}, analyte::B) where {A, B <: A}
    id = find_isd(tbl, analyte)
    id < 1 && return ones(Float64, length(tbl.samples))
    get_analyte(tbl, id)
end
get_isd(tbl::ColumnAnalysisTable, id::Int) = id < 1 ? ones(Float64, length(tbl.samples)) : get_analyte(tbl, id)
function get_isd(tbl::ColumnAnalysisTable{A}, analyte::B) where {A, B <: A}
    id = find_isd(tbl, analyte)
    id < 1 && return ones(Float64, length(tbl.samples))
    get_analyte(tbl, id)
end
get_isd(tbl::SampleWrapperTable, id::Int) = get_isd(tbl.analysistable, id)
get_isd(tbl::SampleWrapperTable{A}, analyte::B) where {A, B <: A} = get_isd(tbl.analysistable, analyte)

"""
    set_isd!(tbl::RowAnalysisTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A}
    set_isd!(tbl::ColumnAnalysisTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A}
    set_isd!(tbl::SampleWrapperTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A}
    set_isd!(tbl::CalWrapperTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A}

Set internal standard of `analyte` to `isd`.
"""
function set_isd!(tbl::RowAnalysisTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A}
    isnothing(tbl.isd_map) && throw(ArgumentError("The table was constructed without internal standards. Please reconstruct with isd information."))
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    iid = findfirst(==(isd), tbl.analytes)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the table"))
    getproperty(tbl.table, tbl.isd_map)[id] = iid
    tbl
end
function set_isd!(tbl::ColumnAnalysisTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A}
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    iid = findfirst(==(isd), tbl.analytes)
    isnothing(iid) && throw(ArgumentError("Analyte $isd is not in the table"))
    tbl.isd_map[id] = iid
    tbl
end
set_isd!(tbl::SampleWrapperTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A} = (set_isd!(tbl.analysistable, analyte, isd); tbl)
function set_isd!(tbl::CalWrapperTable{A}, analyte::B, isd::C) where {A, B <: A, C <: A}
    set_isd!(tbl.conctable, analyte, isd)
    isnothing(tbl.signaltable) || isnothing(tbl.signaltable.isd_map) || set_isd!(tbl.signaltable, analyte, isd)
    tbl
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

Return string representation of `cal.weight` or `weight`. "none" for 0, "1/√x" for -0.5, "1/x" for -1, "1/x²" for -2, 
"x^\$`weight`" for positive `weight`, and "x^(\$`weight`)" for negative `weight`
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
    "x^($weight)"
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
    parse(Float64, first(match(r"x^\(*(-*\d*)\)*", weight)))
end

"""
    formula_repr(cal::SingleCalibration; digits = nothing, sigdigits = 4)
    formula_repr(cal::MultipleCalibration; digits = nothing, sigdigits = 4)

Return string representation of formula with specified `digits` and `sigdigits`. See `format_number`.
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

Return string representation of formula for text file output.
"""
formula_repr_utf8(cal::AbstractCalibration) = replace(formula_repr(cal), "x²" => "x^2")
"""
    weight_repr_utf8(cal::AbstractCalibration)

Return string representation of weight for text file output.
"""
weight_repr_utf8(cal::AbstractCalibration) = replace(weight_repr(cal), "x²" => "x^2", "√x" => "x^0.5")
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
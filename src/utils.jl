get_formula(cal::MultipleCalibration) = get_formula(cal.type, cal.zero)
get_formula(type::Bool, zero::Bool) = if type 
    zero ? @formula(y ~ 0 + x) : @formula(y ~ x)
else
    zero ? @formula(y ~ 0 + x + x ^ 2) : @formula(y ~ x + x ^ 2)
end

function get_analyte(tbl::RowAnalysisTable{T}, analyte::T) where T
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    collect(getproperties(tbl.table[id], tbl.sample_name))
end

function get_analyte(tbl::ColumnAnalysisTable{T}, analyte::T) where T
    id = findfirst(==(analyte), tbl.analytes)
    isnothing(id) && throw(ArgumentError("Analyte $analyte is not in the table"))
    getproperty(tbl.table, tbl.analyte_name[id])
end

function get_sample(tbl::RowAnalysisTable{T}, sample::Symbol) where T
    id = findfirst(==(sample), tbl.samples)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    collect(getproperties(tbl.table[id], tbl.sample_name))
end

function get_sample(tbl::ColumnAnalysisTable{T}, sample::Symbol) where T
    id = findfirst(==(sample), tbl.samples)
    isnothing(id) && throw(ArgumentError("Sample $sample is not in the table"))
    getproperty(tbl.table, tbl.analyte_name[id])
end

function critical_point(cal::MultipleCalibration)
    β = cal.model.model.pp.beta0
    c, b, a = cal.zero ? (0, β...) : β
    -b / 2a
end

cal_range(project::Project) = cal_range(project.calibration)
cal_range(cal::MultipleCalibration) = (lloq(cal), uloq(cal))
lloq(project::Project) = lloq(project.calibration)
uloq(project::Project) = uloq(project.calibration)
lloq(cal::MultipleCalibration) = (cal.type || last(cal.model.model.pp.beta0) < 0) ? cal.table.x[findfirst(cal.table.include)] : max(cal.table.x[findfirst(cal.table.include)], critical_point(cal))
uloq(cal::MultipleCalibration) = (cal.type || last(cal.model.model.pp.beta0) > 0) ? cal.table.x[findlast(cal.table.include)] : min(cal.table.x[findlast(cal.table.include)], critical_point(cal))


function weight_repr(cal::MultipleCalibration)
    cal.weight in [-0.5, -1, -2] || (cal.weight = 0)
    weight_repr(cal.weight)
end
weight_repr(weight::Number) = if weight == -0.5
    "1/√x"
elseif weight == -1
    "1/x"
elseif weight == -2
    "1/x²"
else
    "none"
end

weight_value(weight) = if weight == "1/√x"
    -0.5
elseif weight == "1/x"
    -1
elseif weight == "1/x²"
    -2
else
    0
end

formula_repr(cal::SingleCalibration) = "y = $(round(; sigdigits = 4))x"
function formula_repr(cal::AbstractCalibration)
    β = cal.model.model.pp.beta0
    cal.type && cal.zero && return "y = $(round(β[1]; sigdigits = 4))x"
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    if cal.type
        string("y = ", format_number(β[1]), op[1], abs(format_number(β[2])), "x")
    elseif cal.zero
        string("y = ", format_number(β[1]), "x", op[1], abs(format_number(β[2])), "x²")
    else
        string("y = ", format_number(β[1]), op[1], abs(format_number(β[2])), "x", op[2], abs(format_number(β[3])), "x²")
    end
end

formula_repr_utf8(cal::AbstractCalibration) = replace(formula_repr(cal), "x²" => "x^2")
weight_repr_utf8(cal::AbstractCalibration) = replace(weight_repr(cal), "x²" => "x^2", "√x" => "x^0.5")
format_number(x; digits) = format_number2int(round(x; digits))
format_number(x; sigdigits = 4) = format_number2int(round(x; sigdigits))
format_number2int(x) = 
    x == round(x) ? round(Int, x) : x

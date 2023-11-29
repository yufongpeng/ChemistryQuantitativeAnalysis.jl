
inv_predict(project::Project, tbl) = inv_predict(project.calibration, tbl)
inv_predict_cal!(project::Project) = (inv_predict_cal!.(project.calibration); project)
inv_predict_cal!(cal::SingleCalibration) = cal
function inv_predict_cal!(cal::MultipleCalibration)
    cal.table.x̂ .= inv_predict(cal, cal.table.y)
    cal
end
inv_predict(cal::MultipleCalibration, tbl::AbstractAnalysisTable) = inv_predict(cal, get_analyte(tbl, cal.source.analytes[cal.analyte]))
inv_predict(cal::SingleCalibration, tbl::AbstractAnalysisTable) = inv_predict(cal, get_analyte(tbl, cal.source.analytes[cal.analyte]))
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
inv_predict(cal::SingleCalibration, y::AbstractArray) = y ./ cal.β

accuracy(project::Project, tbl = project.calibration.table) = accuracy(project.calibration, tbl)
accuracy(cal::MultipleCalibration, tbl = cal.table.y) = accuracy(inv_predict(cal, tbl), tbl.x)
accuracy(cal::SingleCalibration, tbl) = accuracy(inv_predict(cal, tbl), tbl.x)
accuracy!(project::Project) = (accuracy!(project.calibration); project)
function accuracy!(cal::MultipleCalibration)
    cal.table.accuracy .= accuracy(cal.table.x̂, cal.table.x)
    cal
end
accuracy!(cal::SingleCalibration) = cal
accuracy(x̂::AbstractVector, x::AbstractVector) = @. x̂ / x

inv_predict_accuracy! = accuracy! ∘ inv_predict_cal!

function calibration(tbl::Table, source::AbstractAnalysisTable; 
                    analyte = 1,
                    isd = 0,
                    type = true, 
                    zero = false, 
                    weight = 0
                    )
    id = findall(x -> isa(x, Number), tbl.y)
    tbl = tbl[id]
    table = :id in propertynames(tbl) ? tbl : Table((; id = collect(1:length(tbl)), ), tbl)
    table = :include in propertynames(tbl) ? table : Table(table; include = trues(length(table)))
    f = get_formula(type, zero)
    model = calfit(table, f, type, zero, weight)
    xlevel = unique(table.x)
    table = Table(; id = table.id, level = [findfirst(x -> i == x, xlevel) for i in table.x], y = table.y, x = table.x, x̂ = zeros(Float64, length(table)), accuracy = zeros(Float64, length(table)), include = table.include)
    cal = Calibration(analyte, isd, type, zero, weight, f, source, table, model)
    inv_predict_accuracy!(cal)
    cal
end

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

function calfit!(cal::MultipleCalibration)
    cal.model = calfit(cal.table, cal.formula, cal.type, cal.zero, cal.weight)
    cal
end

function project(cal::AbstractAnalysisTable, sample = nothing;
                isd = 0,
                type = true, 
                zero = false, 
                weight = 0
                )
    Project(
        map(enumerate(cal.analytes)) do (id, analyte)
            tbl = Table(;id = cal.samples, x = get_analyte(cal, analyte))
            calibration(tbl, cal; analyte = id, isd, type, zero, weight)
        end,
        sample
    )
end
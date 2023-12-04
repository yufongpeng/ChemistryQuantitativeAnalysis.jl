"""
    inv_predict(project::Project, tbl::SampleWrapperTable)
    inv_predict(cal::SingleCalibration)
    inv_predict(cal::MultipleCalibration)
    inv_predict(cal::SingleCalibration, y::AbstractArray)
    inv_predict(cal::MultipleCalibration, y::AbstractArray)

Apply inverse prediction to all analytes in `tbl` or analyte specified in `cal` with `cal.table.y` or `y` as signals. The result will be a `SampleWrapperTable` or a `Vector`.
"""
function inv_predict(project::Project, tbl::SampleWrapperTable)
    cal_id = [i > 0 ? findfirst(c -> ==(project.caltable.conctable.analytes[c.analyte], tbl.analytes[i]), project.calibration) : 0 for i in tbl.cal_map]
    replace!(cal_id, nothing => 0)
    cs = Vector{Vector{Float64}}(undef, length(cal_id))
    Threads.@threads for i in eachindex(cal_id)
        cal_id[i] > 0 || (cs[i] = zeros(Float64, length(tbl.samples)); continue)
        cal = project.calibration[cal_id[i]]
        cs[i] = inv_predict(cal, tbl.analysistable; analyte = tbl.analytes[i])
    end
    SampleWrapperTable(tbl.cal_map, fill_sample_prediction!(tbl.analysistable, cs))
end

inv_predict(cal::SingleCalibration) = get_analyte(cal.caltable.conctable, cal.isd)
inv_predict(cal::MultipleCalibration) = inv_predict(cal, cal.table.y)

function inv_predict(cal::AbstractCalibration, tbl::AbstractAnalysisTable; analyte = cal.caltable.conctable.analytes[cal.analyte])
    cal.isd > 0 || return inv_predict(cal, get_analyte(tbl, analyte))
    inv_predict(cal, get_analyte(tbl, analyte) ./ get_analyte(tbl, cal.caltable.conctable.analytes[cal.isd]))
end
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
inv_predict(cal::SingleCalibration, y::AbstractArray) = y .* first(get_analyte(cal.caltable.conctable, cal.isd))

"""
    inv_predict!(project::Project, tbl::SampleWrapperTable)
    inv_predict!(cal::SingleCalibration)
    inv_predict!(cal::MultipleCalibration)

Inplace inverse prediction which replacing `project.results` or `cal.table.x̂` and returning `project.results` and `cal` respectively.
"""
function inv_predict!(project::Project, tbl::SampleWrapperTable)
    project.resulttable = inv_predict(project, tbl)
end
inv_predict!(cal::SingleCalibration) = cal
function inv_predict!(cal::MultipleCalibration)
    cal.table.x̂ .= inv_predict(cal, cal.table.y)
    cal
end

"""
    inv_predict_sample!(project::Project)
    inv_predict_sample!(project::Project, tbl::SampleWrapperTable)

Apply `inv_predict!` to `project.sampletable`. If a new `SampleWrapperTable` is provided, it updates `project.sampletable` first.
"""
inv_predict_sample!(project::Project) = inv_predict!(project, project.sampletable)
function inv_predict_sample!(project::Project, tbl::SampleWrapperTable)
    project.sampletable = tbl
    inv_predict!(project, project.sampletable)
end

function fill_sample_prediction!(tbl::ColumnAnalysisTable, prediction::Vector{Vector{Float64}})
    result = deepcopy(tbl)
    for (a, c) in zip(tbl.analyte_name, prediction)
        getproperty(result.table, a) .= c
    end
    result
end
function fill_sample_prediction!(tbl::RowAnalysisTable, prediction::Vector{Vector{Float64}})
    result = deepcopy(tbl)
    for (i, p) in enumerate(tbl.sample_name)
        id = findfirst(==(p), propertynames(tbl.table))
        getproperty(result.table, propertynames(tbl.table)[id]) .= getindex.(prediction, i)
    end
    result
end 

"""
    inv_predict_cal!(project::Project)

Apply `inv_predict!` to each calibration curve and return `project`.
"""
inv_predict_cal!(project::Project) = (inv_predict!.(project.calibration); project)

"""
    accuracy(project::Project, tbl = project.calibration.table) 
    accuracy(cal::MultipleCalibration, tbl = cal.table)
    accuracy(cal::SingleCalibration, tbl)
    accuracy(x̂::AbstractVector, x::AbstractVector)
    accuracy!(project::Project)
    accuracy!(cal::MultipleCalibration)
    accuracy!(cal::SingleCalibration)

Calculate accuracy. Return the calculated accuracy for non-mutating version, and input object for mutating version. Nothing will be changed for `SingleCalibration` because it does not store accuracy.
"""
accuracy(project::Project, tbl = project.calibration.table) = accuracy(project.calibration, tbl)
accuracy(cal::MultipleCalibration, tbl = cal.table) = accuracy(inv_predict(cal, tbl.y), tbl.x)
accuracy(cal::SingleCalibration, tbl) = accuracy(inv_predict(cal, tbl.y), tbl.x)
accuracy(x̂::AbstractVector, x::AbstractVector) = @. x̂ / x

"""
    accuracy(project::Project, tbl = project.calibration.table) 
    accuracy(cal::MultipleCalibration, tbl = cal.table)
    accuracy(cal::SingleCalibration, tbl)
    accuracy(x̂::AbstractVector, x::AbstractVector)
    accuracy!(project::Project)
    accuracy!(cal::MultipleCalibration)
    accuracy!(cal::SingleCalibration)

Calculate accuracy. Return the calculated accuracy for non-mutating version, and input object for mutating version. Nothing will be changed for `SingleCalibration` because it does not store accuracy.
"""
accuracy!(project::Project) = (accuracy!(project.calibration); project)
function accuracy!(cal::MultipleCalibration)
    cal.table.accuracy .= accuracy(cal.table.x̂, cal.table.x)
    cal
end
accuracy!(cal::SingleCalibration) = cal

"""
    const inv_predict_accuracy! = accuracy! ∘ inv_predict!

Apply `inv_predict!` and `accuracy!` subsequantly to an `AbstractCalibration`.
"""
const inv_predict_accuracy! = accuracy! ∘ inv_predict!
"""
    const inv_predict_cal_accuracy! = accuracy! ∘ inv_predict!

Apply `inv_predict_cal!` and `accuracy!` subsequantly to a `Project`.
"""
const inv_predict_cal_accuracy! = accuracy! ∘ inv_predict_cal!

"""
    calibration(tbl::Table, caltable::CalWrapperTable; analyte = 1, isd = 0, type = true, zero = false, weight = 0)
    calibration(project::Project{A}, analyte::B; id = project.caltable.signaltable.samples,
                        type = true, zero = false, weight = 0) where {A, B <: A}
    calibration(caltable::CalWrapperTable{A}, analyte::B; id = caltable.signaltable.samples, type = true, zero = false, weight = 0) where {A, B <: A}

Create `MultipleCalibration`.

* `tbl`: `TypedTable.Table` which will be stored as `table`. It is clean-up calibration data for points selection.
* `caltable`: `CalWrapperTable`, calibration data which will be stored as `caltable`. 
* `project`: `Project`.

# Keyword arguments
* `analyte`: `Int` which will be stored as `analyte`.
* `isd`: `Int` which will be stored as `isd`.
* `type`: `Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`), which will be stored as `type`.
* `zero`: `Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`), which will be stored as `zero`.
* `weight`: `Float64` represents the exponential applying to each element of `x` as a weighting vector, which will be stored as `weight`.
"""
function calibration(tbl::Table, caltable::CalWrapperTable; 
                    analyte = 1,
                    isd = 0,
                    type = true, 
                    zero = false, 
                    weight = 0
                    )
    id = findall(x -> isa(x, Number), tbl.y)
    tbl = tbl[id]
    xlevel = unique(tbl.x)
    table = :id in propertynames(tbl) ? tbl : Table((; id = collect(1:length(tbl)), ), tbl)
    table = :level in propertynames(tbl) ? table : Table(table; level = [findfirst(x -> i == x, xlevel) for i in table.x])
    table = :include in propertynames(tbl) ? table : Table(table; include = trues(length(table)))
    f = get_formula(type, zero)
    model = calfit(table, f, type, zero, weight)
    table = Table(; id = table.id, level = table.level, y = table.y, x = table.x, x̂ = zeros(Float64, length(table)), accuracy = zeros(Float64, length(table)), include = table.include)
    cal = MultipleCalibration(analyte, isd, type, zero, weight, f, caltable, table, model)
    inv_predict_accuracy!(cal)
    cal
end
calibration(project::Project{A}, analyte::B; 
                    id = project.caltable.signaltable.samples,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B <: A} = calibration(project.caltable, analyte; id, type, zero, weight)
function calibration(caltable::CalWrapperTable{A}, analyte::B; 
                    id = caltable.signaltable.samples,
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, B <: A}

    ord = sortperm(caltable.level_map)
    iid = find_isd(caltable.conctable, analyte)
    aid = find_analyte(caltable.conctable, analyte)
    level = caltable.level_map[ord]
    conc = get_analyte(caltable.conctable, analyte)
    table = Table(; 
                    id = id[ord], 
                    level = level, 
                    x = map(level) do l
                        conc[find_sample(caltable.conctable, Symbol(l))]
                    end, 
                    y = (get_analyte(caltable.signaltable, analyte) ./ get_isd(caltable.signaltable, iid))[ord], 
                    x̂ = zeros(Float64, length(id)), 
                    accuracy = zeros(Float64, length(id)),
                    include = trues(length(id)))
    f = get_formula(type, zero)
    model = calfit(table, f, type, zero, weight)
    cal = MultipleCalibration(aid, iid, type, zero, Float64(weight), f, caltable, table, model)
    inv_predict_accuracy!(cal)
    cal
end

"""
    calfit(tbl, formula, type, zero, weight)
    calfit(cal::MultipleCalibration)
    calfit!(project::Project)
    calfit!(cal::MultipleCalibration)

Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

Field `model` will be mutated for mutating version. Calling `calfit!` on a `Project` will apply `calfit!` to all calibration curves.

It returns `GLM` object for non-mutating version and input object for mutating version.
"""
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
calfit(cal::MultipleCalibration) = calfit(cal.table, cal.formula, cal.type, cal.zero, cal.weight)

"""
    calfit(tbl, formula, type, zero, weight)
    calfit(cal::MultipleCalibration)
    calfit!(project::Project)
    calfit!(cal::MultipleCalibration)

Fit a `GLM` model based on provided `formula`, `type`, `zero` and `weight` or parameters from calibration curves. 

Field `model` will be mutated for mutating version. Calling `calfit!` on a `Project` will apply `calfit!` to all calibration curves.

It returns `GLM` object for non-mutating version and input object for mutating version.
"""
calfit!(project::Project) = (calfit!.(project.calibration); project)
function calfit!(cal::MultipleCalibration)
    cal.model = calfit(cal.table, cal.formula, cal.type, cal.zero, cal.weight)
    cal
end

"""
    switch_isd!(project::Project{A}, analyte::B, isd::C) where {A, B <: A, C <: A}

Switch internal standard of `analyte` to `isd`.
"""
function switch_isd!(project::Project{A}, analyte::B, isd::C) where {A, B <: A, C <: A}
    aid = find_analyte(project.caltable.conctable, analyte)
    id = findfirst(x -> ==(x.analyte, aid), project.calibration)
    isnothing(id) && throw(ArgumentError("No fitted calibration data for $analyte"))
    set_isd!(project.caltable, analyte, isd)
    project.calibration[id] = calibration(project.caltable, analyte; 
            type = project.calibration[id].type,
            zero = project.calibration[id].zero,
            weight = project.calibration[id].weight)
    isnothing(project.sampletable) || set_isd!(project.sampletable, analyte, isd)
    @warn "Call `inv_predict` on samples to get the updated results."
    project
end

"""
    project(cal::CalWrapperTable{A, T}, sample = nothing, result = nothing;
                    type = true, 
                    zero = false, 
                    weight = 0
                    ) where {A, T}

Create `Project` from `cal`, and optionally `sample`, `result` with specified calibration parameters. See "MultipleCalibration" for detail description of keyword arguments.
"""
function project(cal::CalWrapperTable{A, T}, sample = nothing, result = nothing;
                type = true, 
                zero = false, 
                weight = 0
                ) where {A, T}
    Project(
        length(cal.conctable.samples) > 1 ? map(cal.conctable.analytes) do analyte
            calibration(cal, analyte; type, zero, weight)
        end : map(enumerate(cal.conctable.analytes)) do (id, analyte)
            SingleCalibration(id, cal.isd, cal)
        end
        ,
        cal,
        sample,
        result
    )
end
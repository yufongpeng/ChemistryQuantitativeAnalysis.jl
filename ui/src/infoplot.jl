function viewinfo(cal::MultipleCalibration, at, method::AnalysisMethod; 
                rel_sig = :relative_signal,
                est_conc = :estimated_concentration, 
                lloq_multiplier = 4//3, dev_acc = 0.15)
    tbl = cal.table
    color1 = [i ? "honeydew" : "darkseagreen" for i in tbl.include]
    loq_level = tbl.level[findfirst(tbl.include)]
    ft = @view tbl[tbl.level .> loq_level]
    lt = @view tbl[tbl.level .<= loq_level]
    color3 = vcat(
            [i ? (abs(j - 1) <= lloq_multiplier * dev_acc ? "honeydew" : "lightpink") : "darkseagreen" for (i, j) in zip(lt.include, lt.accuracy)],
            [i ? (abs(j - 1) <= dev_acc ? "honeydew" : "lightpink") : "darkseagreen" for (i, j) in zip(ft.include, ft.accuracy)]
            )
    (isnothing(at) || isempty(at)) && return begin
        PlotlyJS.plot(
            table(
                header = attr(values = ["Sample", "Level", "Y", "X", "Predicted X", "Accuracy"],
                            line_color = "darkgreen",
                            fill_color = "limegreen",
                            align = "center",
                            font = attr(color = "white")),
                cells = attr(
                    values = [tbl.id, tbl.level, round.(tbl.y; sigdigits = 4), round.(tbl.x; sigdigits = 4), round.(tbl.x̂; sigdigits = 4), round.(tbl.accuracy; sigdigits = 4)],
                    line_color = "darkgreen",
                    fill_color = vcat(repeat([c], 5), [i]),
                    align = "right"
                )
            )
        )
    end
    j = findfirst(==(first(cal.analyte)), method.analyte)
    js = findall(==(j), method.analytetable.calibration)
    id = convert(Vector{Any}, tbl.id)
    level = convert(Vector{Any}, tbl.level)
    y = convert(Vector{Any}, round.(tbl.y; sigdigits = 4))
    x = convert(Vector{Any}, round.(tbl.x; sigdigits = 4))
    x̂ = convert(Vector{Any}, round.(tbl.x̂; sigdigits = 4))
    acc = convert(Vector{Any}, round.(tbl.accuracy; sigdigits = 4))
    cs = [id, level, y, x, x̂, acc]
    lloq, uloq = dynamic_range(cal)
    color2 = deepcopy(color1)
    for j in js
        analyte = method.analyte[j]
        ay = getanalyte(getproperty(at, rel_sig), analyte)
        ax = getanalyte(getproperty(at, est_conc), analyte)
        push!(cs[1], analyte)
        push!(color1, "rgb(235, 193, 238)")
        push!(color2, "rgb(235, 193, 238)")
        push!(color3, "rgb(235, 193, 238)")
        for v in @view cs[2:end]
            push!(v, "")
        end
        append!(cs[1], sampleobj(at))
        append!(cs[2], repeat([""], length(sampleobj(at))))
        append!(cs[3], round.(ay; sigdigits = 4))
        append!(cs[4], repeat([""], length(sampleobj(at))))
        append!(cs[5], round.(ax; sigdigits = 4))
        append!(cs[6], repeat([""], length(sampleobj(at))))
        c2 = [(i >= lloq && i <= uloq) ? "honeydew" : "lightpink" for i in ax]
        append!(color1, repeat(["honeydew"], length(sampleobj(at))))
        append!(color2, c2)
        append!(color3, repeat(["honeydew"], length(sampleobj(at))))
    end
    PlotlyJS.plot(
            table(
                header = attr(values = ["Sample", "Level", "Y", "X", "Predicted X", "Accuracy"],
                            line_color = "darkgreen",
                            fill_color = "limegreen",
                            align = "center",
                            font = attr(color = "white")),
                cells = attr(
                    values = cs,
                    line_color = "darkgreen",
                    fill_color = vcat(repeat([color1], 4),[color2, color3]),
                    align = "right"
                )
            )
        )

end
function view_cal(tbl::Table; lloq_multiplier = 4//3, dev_acc = 0.15)
    c = [i ? "honeydew" : "darkseagreen" for i in tbl.include]
    loq_level = tbl.level[findfirst(tbl.include)]
    ft = @view tbl[tbl.level .> loq_level]
    lt = @view tbl[tbl.level .<= loq_level]
    i = vcat(
            [i ? (abs(j - 1) <= lloq_multiplier * dev_acc ? "honeydew" : "lightpink") : "darkseagreen" for (i, j) in zip(lt.include, lt.accuracy)],
            [i ? (abs(j - 1) <= dev_acc ? "honeydew" : "lightpink") : "darkseagreen" for (i, j) in zip(ft.include, ft.accuracy)]
            )
    Plotly.plot(
        table(
            header = attr(values = ["Calibration", "level", "y", "x", "predicted x", "accuracy"],
                        line_color = "darkgreen",
                        fill_color = "limegreen",
                        align = "center"),
            cells = attr(
                values = [tbl.id, tbl.level, round.(tbl.y; sigdigits = 4), round.(tbl.x; sigdigits = 4), round.(tbl.x̂; sigdigits = 4), round.(tbl.accuracy; sigdigits = 4)],
                line_color = "darkgreen",
                fill_color = vcat(repeat([c], 5), [i]),
                align = "right"
            )
        )
    )
end

function view_sample(at::AnalysisTable, analyte; lloq, uloq, lloq_multiplier = 4//3, dev_acc = 0.15)
    x = getanalyte(at.estimated_concentration, analyte)
    y = getanalyte(at.relative_signal, analyte)
    id = at.sample
    c = [(i >= lloq * (1 - dev_acc * lloq_multiplier) && i <= uloq * (1 + dev_acc)) ? "honeydew" : "lightpink" for i in x]
    Plotly.plot(
        table(
            header = attr(values = ["Sample", "y", "predicted x"],
                        line_color = "darkgreen",
                        fill_color = "limegreen",
                        align = "center"),
            cells = attr(
                values = [id, round.(y; sigdigits = 4), round.(x; sigdigits = 4)],
                line_color = "darkgreen",
                fill_color = [c],
                align = "right"
            )
        )
    )
end
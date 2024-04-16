"""
    interactive_calibrate!(batch::Batch; 
                signal = batch.method.signal,
                rel_sig = :relative_signal,
                est_conc = :estimated_concentration, 
                root = pwd(), 
                lloq_multiplier = 4//3, dev_acc = 0.15,
                fig_attr = Dict{Symbol, Any}(:resolution => (1350, 900)), 
                axis_attr = cal -> Dict(:title => string(first(cal.analyte)), :xlabel => "Concentration (nM)", :ylabel => "Abundance", :titlesize => 20), 
                plot_attr = cal -> Dict(
                                :scatter => Dict(
                                    :color => [:blue, :red], 
                                    :inspector_label => (self, i, p) -> string("id: ", cal.table.id[i], "\nlevel: ", cal.table.level[i], "\naccuracy: ", round(cal.table.accuracy[i]; sigdigits = 4))), 
                                :line => Dict(:color => :chartreuse))
                )

Interactively calibrate signal and concentration.

# Arguments
* `batch`: `Batch`.
* `signal`: `Symbol`, signal data property.
* `rel_sig`: `Symbol`, relative signal data property.
* `est_conc`: `Symbol`, estimated concentration data property.
* `root`: root directory for saving objects.
* `dev_acc`: allowed accuracy deviation.
* `lloq_multiplier`: multiplier for `dev_acc` at LLOQ level.
* `fig_attr`: `Dict` contains attributes of `Figure`.
* `axis_attr`: a function map each calibtation curve to a `Dict` containing its attributes of `Axis`.
* `plot_attr`: a function map each calibtation curve to a `Dict` containing its attributes of plots (`Line` or `Scatter`).
"""
function interactive_calibrate!(batch::Batch; 
                signal = batch.method.signal,
                rel_sig = :relative_signal,
                est_conc = :estimated_concentration, 
                root = pwd(), 
                lloq_multiplier = 4//3, dev_acc = 0.15,
                fig_attr = Dict{Symbol, Any}(:resolution => (1350, 900)), 
                axis_attr = cal -> Dict(:title => string(first(cal.analyte)), :xlabel => "Concentration (nM)", :ylabel => "Abundance", :titlesize => 20), 
                plot_attr = cal -> Dict(
                                :scatter => Dict(
                                    :color => [:blue, :red], 
                                    :inspector_label => (self, i, p) -> string("id: ", cal.table.id[i], 
                                                                                "\nlevel: ", cal.table.level[i], 
                                                                                "\naccuracy: ", round(cal.table.accuracy[i]; sigdigits = 4))
                                                ), 
                                :line => Dict(:color => :chartreuse))
                )
    isempty(batch.calibration) && init_calibration!(batch)
    update_quantification!(batch; signal, rel_sig, est_conc)
    fig = Figure(; fig_attr...)
    axis_attrs = map(axis_attr, batch.calibration)
    plot_attrs = map(plot_attr, batch.calibration)
    i = 1
    info = Window()
    function draw()
        label_analyte = Label(fig, string(first(batch.calibration[i].analyte)); halign = :left, width = 250)
        button_left = Button(fig, label = "<")
        button_right = Button(fig, label = ">")
        label_r2 = Label(fig, "R² = $(round(r2(batch.calibration[i].model); sigdigits = 4))"; halign = :left)
        label_formula = Label(fig, formula_repr(batch.calibration[i]); halign = :left)
        menu_type = Menu(fig, options = ["linear", "quadratic"], default = batch.calibration[i].type ? "linear" : "quadratic", tellwidth = false)
        menu_zero = Menu(fig, options = ["ignore (0, 0)", "include (0, 0)"], default = batch.calibration[i].zero ? "include (0, 0)" : "ignore (0, 0)", tellwidth = false)
        default_w = weight_repr(batch.calibration[i])
        menu_wt = Menu(fig, options = ["none", "1/√x", "1/x", "1/x²"], default = default_w)
        menu_zoom = Menu(fig, options = string.(0:length(unique(batch.calibration[i].table.x))), default = "0", tellwidth = false)
        # menu_show = Menu(fig, options = ["Cal", "Sample"], default = "Cal"; halign = :left, tellwidth = false)
        menu_export = Menu(fig, options = ["Fig", "Cal", "Fig&Cal"], default = "Cal"; halign = :left, tellwidth = false)
        ax = Axis(fig[1, 1]; axis_attrs[i]...)
        sc = scatter!(ax, batch.calibration[i].table.x, batch.calibration[i].table.y; get_point_attr(plot_attrs[i], batch.calibration[i])...)
        DataInspector(sc)
        xlevel = unique(batch.calibration[i].table.x)
        xscale = -reduce(-, extrema(xlevel))
        yscale = -reduce(-, extrema(batch.calibration[i].table.y))
        xrange = Table(; x = collect(LinRange(extrema(xlevel)..., convert(Int, reduce(-, extrema(xlevel)) ÷ maximum(xlevel[1:end - 1] .- xlevel[2:end]) * 100))))
        ln = lines!(ax, xrange.x, predict(batch.calibration[i].model, xrange); get!(plot_attrs[i], :line, Dict(:color => :chartreuse))...)
        tall = Menu(fig, options = ["All plots", "This plot"], default = "All plots")
        objs = Dict(:axis => ax, :scatter => sc, :line => ln)
        menu_obj = Menu(fig, options = collect(keys(objs)), default = "axis", halign = :left)
        button_confirm = Button(fig, label = "confirm", halign = :left)    
        textbox_attr = Textbox(fig, placeholder = "attribute", tellwidth = false, halign = :left)
        textbox_value = Textbox(fig, placeholder = "value (julia expression)", tellwidth = false, halign = :left)
        button_show = Button(fig, label = "Table")
        button_export = Button(fig, label = "export")
        button_save = Button(fig, label = "Save batch"; halign = :left)
        button_saveas = Button(fig, label = "Save batch as")
        lzoom = Label(fig, "Zoom")
        ltype = Label(fig, "Type")
        lzero = Label(fig, "Zero")
        lweight = Label(fig, "Weight")
        lps = Label(fig, "Plot setting", halign = :left, tellwidth = false)
        lc = [label_analyte, button_left, button_right, label_r2, label_formula, lzoom, menu_zoom, 
            ltype, menu_type, lzero, menu_zero, lweight, menu_wt, lps, tall, menu_obj, textbox_attr, button_confirm, textbox_value, 
            menu_export, button_show, button_export, button_save, button_saveas]
        fig[1, 2] = vgrid!(
            label_analyte, 
            hgrid!(button_left, button_right),
            label_r2,
            label_formula,
            hgrid!(lzoom, menu_zoom),
            hgrid!(ltype, menu_type), 
            hgrid!(lzero, menu_zero),
            hgrid!(lweight, menu_wt),
            hgrid!(lps, tall),
            menu_obj, 
            hgrid!(textbox_attr, button_confirm),
            textbox_value, 
            hgrid!(button_show, menu_export, button_export), 
            hgrid!(button_save, button_saveas);
            tellheight = false
        )
        xr = dynamic_range(batch.calibration[i])
        xr = xr .+ (xr[2] - xr[1]) .* (-0.05, 0.05)
        yr = signal_range(batch.calibration[i]) .* ((1 - dev_acc * lloq_multiplier), (1 + dev_acc))
        yr = yr .+ (yr[2] - yr[1]) .* (-0.05, 0.05)
        limits!(ax, xr, yr)
        body!(info, viewinfo(batch.calibration[i], batch.data, batch.method; rel_sig, est_conc, lloq_multiplier, dev_acc))
        #display(view_sample(sample; lloq = batch.calibration[i].table.x[findfirst(batch.calibration[i].table.include)], uloq = batch.calibration[i].table.x[findlast(batch.calibration[i].table.include)], lloq_multiplier, dev_acc))
        # Main.vscodedisplay(batch.calibration[i].table[batch.calibration[i].table.include])
        # fig[1, 3] = vgrid!(map(s -> Label(fig, s; halign = :left), split(sprint(showtable, batch.calibration[i].table), "\n"))...; tellheight = false, width = 250)
        function update!()
            update_calibration!(batch.calibration[i], batch.method)
            update_quantification!(batch)
            ln.input_args[2][] = predict(batch.calibration[i].model, xrange)
            label_r2.text = "R² = $(round(r2(batch.calibration[i].model); sigdigits = 4))"
            label_formula.text = formula_repr(batch.calibration[i])
            #sample.x̂ .= inv_predict(batch.calibration[i], sample)
            body!(info, viewinfo(batch.calibration[i], batch.data, batch.method; rel_sig, est_conc, lloq_multiplier, dev_acc))
        end
        function update_analyte!()
            for x in lc
                delete!(x)
            end
            delete!(ax)
            draw()
        end
        on(button_left.clicks) do s
            j = max(i - 1, firstindex(batch.calibration))
            if j != i
                i = j
                update_analyte!()
            end
        end
        on(button_right.clicks) do s
            j = min(i + 1, lastindex(batch.calibration))
            if j != i
                i = j
                update_analyte!()
            end
        end
        on(events(ax).mousebutton) do event
            if event.action == Mouse.press
                plot, id = pick(ax)
                if id != 0 && plot == sc
                    if event.button == Mouse.left
                        batch.calibration[i].table.include[id] = !batch.calibration[i].table.include[id]
                        delete!(ax, sc)
                        sc = scatter!(ax, batch.calibration[i].table.x, batch.calibration[i].table.y; get_point_attr(plot_attrs[i], batch.calibration[i])...)
                        DataInspector(sc)
                        update!()
                    end
                end
            end
            return Consume(false)
        end
        on(menu_type.selection) do s
            batch.calibration[i].type = s == "linear"
            update!()
        end
        on(menu_zero.selection) do s
            batch.calibration[i].zero = s == "include (0, 0)"
            update!()
        end
        on(menu_wt.selection) do s
            batch.calibration[i].weight = s == "none" ? 0 : weight_value(s) # CQA bug
            update!()
        end
        on(menu_zoom.selection) do s
            s = parse(Int, s)
            if s == 0
                #autolimits!(ax)
                xr = dynamic_range(batch.calibration[i])
                xr = xr .+ (xr[2] - xr[1]) .* (-0.05, 0.05)
                yr = signal_range(batch.calibration[i]) .* ((1 - dev_acc * lloq_multiplier), (1 + dev_acc))
                yr = yr .+ (yr[2] - yr[1]) .* (-0.05, 0.05)
                limits!(ax, xr, yr)
            else
                x_value = xlevel[s] 
                id = findall(==(x_value), batch.calibration[i].table.x)
                y_value = batch.calibration[i].table.y[id]
                Δy = length(y_value) == 1 ? 0.2 * y_value[1] : -reduce(-, extrema(y_value))
                yl = extrema(y_value) .+ (-Δy, Δy)
                Δx = Δy * xscale / yscale
                xl = x_value .+ (-Δx, Δx)
                limits!(ax, xl, yl)
            end
        end
        on(button_confirm.clicks) do s
            if menu_obj.selection[] == :scatter
                attr = Symbol(textbox_attr.stored_string[])
                isnothing(attr) && return
                if tall.selection[] == "All plots"
                    for j in eachindex(plot_attrs)
                        plot_attrs[j][:scatter][attr] = eval(Meta.parse(textbox_value.stored_string[]))
                    end
                else
                    plot_attrs[i][:scatter][attr] = eval(Meta.parse(textbox_value.stored_string[]))
                end
                delete!(ax, sc)
                sc = scatter!(ax, batch.calibration[i].table.x, batch.calibration[i].table.y; get_point_attr(plot_attrs[i], batch.calibration[i])...)
                DataInspector(sc)
                return
            elseif menu_obj.selection[] == :line
                attr = Symbol(textbox_attr.stored_string[])
                if tall.selection[] == "All plots"
                    for j in eachindex(plot_attrs)
                        plot_attrs[j][:line][attr] = eval(Meta.parse(textbox_value.stored_string[]))
                    end
                else
                    plot_attrs[i][:line][attr] = eval(Meta.parse(textbox_value.stored_string[]))
                end
            elseif menu_obj.selection[] == :axis
                attr = Symbol(textbox_attr.stored_string[])
                if tall.selection[] == "All plots"
                    for j in eachindex(axis_attrs)
                        axis_attrs[j][attr] = eval(Meta.parse(textbox_value.stored_string[]))
                    end
                else
                    axis_attrs[i][attr] = eval(Meta.parse(textbox_value.stored_string[]))
                end
            end
            x = getproperty(objs[menu_obj.selection[]], Symbol(textbox_attr.stored_string[]))[]
            if length(vectorize(x)) > 1 
                setproperty!(objs[menu_obj.selection[]], Symbol(textbox_attr.stored_string[]), repeat([eval(Meta.parse(textbox_value.stored_string[]))], length(x)))
            else
                setproperty!(objs[menu_obj.selection[]], Symbol(textbox_attr.stored_string[]), eval(Meta.parse(textbox_value.stored_string[]))) 
            end
        end
        on(button_show.clicks) do s
            if !active(info) 
                info = Window()
                body!(info, viewinfo(batch.calibration[i], batch.data, batch.method; rel_sig, est_conc, lloq_multiplier, dev_acc))
            end
                #Main.vscodedisplay(batch.calibration[i].table[batch.calibration[i].table.include])
        end
        on(button_export.clicks) do s
            if menu_export.selection[] == "Fig" || menu_export.selection[] == "Fig&Cal"
                save_dialog("Save figure as", nothing, ["*.png"]; start_folder = root) do f
                    f == "" || save(f, fig; update = false)
                end
            end
            if menu_export.selection[] == "Cal" || menu_export.selection[] == "Fig&Cal"
                save_dialog("Save calibration as", nothing, ["*.csv"]; start_folder = root) do f
                    f == "" || CSV.write(f, Table(; formula = [formula_repr_utf8(batch.calibration[i])], weight = [weight_repr_utf8(batch.calibration[i])], LLOQ = [format_number(lloq(batch.calibration[i]))], ULOQ = [format_number(uloq(batch.calibration[i]))],  r_squared = [format_number(r2(batch.calibration[i].model))]))
                end
            end
            #=
            if menu_show_export.selection[] == "All"
                save_dialog("Save as", nothing; start_folder = pwd()) do f
                    f == "" && return
                    mkpath(f)
                    j = i
                    for id in eachindex(plot_attrs)
                        i = id
                        update_analyte!()
                        sleep(1)
                        save(joinpath(f, "plot_$i.png"), fig; update = false)
                        CSV.write(joinpath(f, "cal_$i.csv"), Table(; formula = [formula_repr_utf8(batch.calibration[i])], weight = [weight_repr_utf8(batch.calibration[i])], LLOQ = [format_number(lloq(batch.calibration[i]))], ULOQ = [format_number(uloq(batch.calibration[i]))],  r_squared = [format_number(r2(batch.calibration[i].model))]))
                    end
                    i = j
                    update_analyte!()
                end
            end
            =#
        end
        on(button_saveas.clicks) do s
            save_dialog("Save as", nothing, ["*.batch"]; start_folder = root) do f
                f == "" || ChemistryQuantitativeAnalysis.write(f, batch)
            end
        end
        on(button_save.clicks) do s
            if endswith(root, ".batch")
                ask_dialog("Save batch?\n$root") && ChemistryQuantitativeAnalysis.write(root, batch)
            else
                save_dialog("Save as", nothing, ["*.batch"]; start_folder = root) do f
                    f == "" || ChemistryQuantitativeAnalysis.write(f, batch)
                end
            end
        end
    end
    draw()
    wait(display(fig))
    close(info)
end

#get_point_attr(plot_attr::Dict, incl::Bool) = NamedTuple(k => incl ? v[1] : v[2] for (k, v) in get!(plot_attr, :scatter, Dict(:color => [:blue, :red])))
#get_point_attr(plot_attr::Dict, incl::BitVector) = NamedTuple(k => isa(v, Vector) ? map(inc -> inc ? v[1] : v[2], incl) : v for (k, v) in get!(plot_attr, :scatter, Dict(:color => [:blue, :red])))
#get_point_attr(plot_attr::Dict, incl::Vector{Bool}) = NamedTuple(k => isa(v, Vector) ? map(inc -> inc ? v[1] : v[2], incl) : v for (k, v) in get!(plot_attr, :scatter, Dict(:color => [:blue, :red])))
get_point_attr(plot_attr::Dict, cal::MultipleCalibration) = NamedTuple(k => isa(v, Vector) ? map(inc -> inc ? v[1] : v[2], cal.table.include) : v for (k, v) in get!(plot_attr, :scatter, Dict(:color => [:blue, :red])))
#get_axis_attr(axis_attr::Dict, cal::MultipleCalibration) = NamedTuple(k => v isa Function ? v(cal) : v for (k, v) in axis_attr)

vectorize(x::AbstractVector) = x
vectorize(x) = [x]
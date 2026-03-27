"""
    ChemistryQuantitativeAnalysis.write(file::String, object; delim = "\\t")

Write `object` into 
* ".ical" for `InternalCalibrator`
* ".ecal" for `ExternalCalibrator`
* ".sdt" for `SampleDataTable`
* ".adt" for `AnalyteDataTable`
* ".at" for `AnalysisTable`
* ".am" for `AnalysisMethod`
* ".batch" for `Batch`

`delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

See README.md for the structure of all files.
"""
function write(file::String, batch::Batch; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[delim]\n", escape_string(string(delim)))
    end
    write(joinpath(file, "method.am"), batch.method; delim)
    mkpath(joinpath(file, "calibrator"))
    ft = isnothing(batch.method.signaltable) ? "ical" : "ecal"
    for (i, c) in enumerate(batch.calibrator)
        write(joinpath(file, "calibrator", "$i.$ft"), c; delim)
    end
    isnothing(batch.data) || write(joinpath(file, "data.at"), batch.data; delim)
end

write_datatable(file::String, tbl::AnalyteDataTable; delim = '\t') = write(file * ".adt", tbl; delim)
write_datatable(file::String, tbl::SampleDataTable; delim = '\t') = write(file * ".sdt", tbl; delim)

function write(file::String, tbl::AnalyteDataTable; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[delim]\n", escape_string(string(delim)), "\n\n[Analyte]\n", analytecol(tbl), "\n\n[Sample]\n", join(samplename(tbl), "\n"))
    end
    CSV.write(joinpath(file, "table.txt"), table(tbl); delim)
end

function write(file::String, tbl::SampleDataTable; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[delim]\n", escape_string(string(delim)), "\n\n[Analyte]\n", join(analytename(tbl), "\n"), "\n\n[Sample]\n", samplecol(tbl))
    end
    CSV.write(joinpath(file, "table.txt"), table(tbl); delim)
end

function write(file::String, tbl::AnalysisTable; delim = '\t')
    mkpath(file)
    rm.(readdir(file; join = true); recursive = true)
    for (i, (k, v)) in enumerate(pairs(tbl))
        write_datatable(joinpath(file, "$(i - 1)_$k"), v; delim)
    end
end

function write(file::String, tbl::AnalysisMethod; delim = '\t')
    mkpath(file)
    write_datatable(joinpath(file, "$(tbl.nom_conc)"), tbl.conctable; delim)
    isnothing(tbl.signaltable) || write_datatable(joinpath(file, "$(tbl.signal)"), tbl.signaltable; delim)
    id = nothing
    if isa(tbl.signaltable, SampleDataTable)
        id = findfirst(x -> getproperty(tbl.signaltable, x) == tbl.pointlevel, propertynames(tbl.signaltable))
    end
    open(joinpath(file, "config.txt"), "w+") do config
        isnothing(id) ? Base.write(config, "[signal]\n", tbl.signal, "\n\n[rel_sig]\n", tbl.rel_sig, "\n\n[est_conc]\n", tbl.est_conc, 
            "\n\n[nom_conc]\n", tbl.nom_conc, "\n\n[acc]\n", tbl.acc, "\n\n[delim]\n", escape_string(string(delim)), "\n\n[pointlevel]\n", join(tbl.pointlevel, "\n")) : 
            Base.write(config, "[signal]\n", tbl.signal, "\n\n[rel_sig]\n", tbl.rel_sig, "\n\n[est_conc]\n", tbl.est_conc, "\n\n[nom_conc]\n", 
                        tbl.nom_conc, "\n\n[acc]\n", tbl.acc, "\n\n[delim]\n", escape_string(string(delim)), "\n\n[levelname]\n", propertynames(tbl.signaltable)[id])
    end
    CSV.write(joinpath(file, "analytetable.txt"), tbl.analytetable; delim)
end

function write(file::String, cal::ExternalCalibrator; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[analyte]\n", string(cal.analyte), "\n\n[isd]\n", string(cal.isd), 
                encode_config(cal.model)..., 
                "\n\n[delim]\n", escape_string(string(delim)))
    end
    CSV.write(joinpath(file, "table.txt"), cal.table; delim)
end

function write(file::String, cal::InternalCalibrator; delim = '\t')
    mkpath(file)
    open(joinpath(file, "config.txt"), "w+") do config
        Base.write(config, "[analyte]\n", string(cal.analyte), "\n\n[conc]\n", string(cal.conc))
    end
end

"""
    mkbatch(file::String; 
            delim = '\t',
            data_table = SampleDataTable,
            signal_table = SampleDataTable,
            conc_table = SampleDataTable,
            data_config = Dict{Symbol, Any}(), 
            signal_config = Dict{Symbol, Any}(), 
            conc_config = Dict{Symbol, Any}(), 
            method_config = Dict{Symbol, Any}()
            )

Create a template directory of a batch. See "README.md" for available keys and values for config.

The default table wrapper is `SampleDataTable`, whose default sample column name is "Sample", and default sample column name for levels is "Level".
The default analyte column name for `AnalyteDataTable` is "Analyte".

The priorities of default analytes are `conc_config[:Analyte]`, `signal_config[:Analyte]`, `data_config[:Analyte]`, and lastly `["Analyte1", "Analyte2"]`.
The default samples for data are `["S1", "S2"]`. 
    
The default calibration points are `["C1", "C2", "C3", "C4", "C5"]`. The default calibration levels are `[1, 2, 3, 4, 5]`.

"""
function mkbatch(file::String; 
                delim = '\t',
                data_table = SampleDataTable,
                signal_table = SampleDataTable,
                conc_table = SampleDataTable,
                data_config = Dict{Symbol, Any}(), 
                signal_config = Dict{Symbol, Any}(), 
                conc_config = Dict{Symbol, Any}(), 
                method_config = Dict{Symbol, Any}()
                )
    mkpath(file)
    default_analyte = get(conc_config, :Analyte, get(signal_config, :Analyte, get(data_config, :Analyte, ["Analyte1", "Analyte2"])))
    default_sample = ["S1", "S2"]
    default_point = ["C1", "C2", "C3", "C4", "C5"]
    default_level = [1, 2, 3, 4, 5]
    default_config_c = [
        Dict(:delim => delim, :Analyte => default_analyte, :Sample => "Sample"), 
        Dict(:delim => delim, :Analyte => default_analyte, :Sample => "Sample"), 
        Dict(:delim => delim, :Analyte => default_analyte, :Sample => "Level")
    ]
    default_config_r = [
        Dict(:delim => delim, :Analyte => "Analyte", :Sample => default_sample), 
        Dict(:delim => delim, :Analyte => "Analyte", :Sample => default_point), 
        Dict(:delim => delim, :Analyte => "Analyte", :Sample => default_level)
    ]
    table = Vector{Table}(undef, 3)
    confs = Vector{String}(undef, 3)
    data_config = convert(Dict{Symbol, Any}, data_config)
    signal_config = convert(Dict{Symbol, Any}, signal_config)
    conc_config = convert(Dict{Symbol, Any}, conc_config)
    method_config = convert(Dict{Symbol, Any}, method_config)
    for (i, (tt, config, default_c, default_r)) in enumerate(zip([data_table, signal_table, conc_table], [data_config, signal_config, conc_config], default_config_c, default_config_r))
        default = tt == SampleDataTable ? default_c : default_r
        for (k, v) in default
            get!(config, k, v)
        end
        if tt == SampleDataTable
            #table[i] = string(config[:Sample], config[:delim], join(config[:Analyte], config[:delim]))
            config[:Analyte_merge] = join(config[:Analyte], "\n")
            config[:Sample_merge] = config[:Sample]
        else
            #table[i] = string(config[:Analyte], config[:delim], join(config[:Sample], config[:delim]))
            config[:Sample_merge] = join(config[:Sample], "\n")
            config[:Analyte_merge] = config[:Analyte]
        end
        confs[i] = string("[delim]\n", escape_string(string(config[:delim])), 
                            "\n\n[Analyte]\n", config[:Analyte_merge], 
                            "\n\n[Sample]\n", config[:Sample_merge])
    end
    if data_table == SampleDataTable
        #table[1] = string(table[1], "\n", join((string(x, data_config[:delim], join(repeat(["1.0"], length(data_config[:Analyte])), data_config[:delim])) for x in default_sample), "\n"))
        table[1] = Table(; Symbol(data_config[:Sample]) => default_sample, (Symbol.(data_config[:Analyte]) .=> Ref(repeat([1.0], length(default_sample))))...)
    else
        table[1] = Table(; Symbol(data_config[:Analyte]) => default_analyte, (Symbol.(data_config[:Sample]) .=> Ref(repeat([1.0], length(default_analyte))))...)
    end
    if signal_table == AnalyteDataTable && conc_table == AnalyteDataTable
        default_level = parse.(Int, string.(conc_config[:Sample]))
        default_point = signal_config[:Sample]
        pointlevel = [i > lastindex(default_level) ? default_level[end] : default_level[i] for i in eachindex(default_point)]
    elseif signal_table == AnalyteDataTable
        default_point = signal_config[:Sample]
        default_level = collect(eachindex(default_point))
        pointlevel = default_level
    elseif conc_table == AnalyteDataTable
        default_level = parse.(Int, string.(conc_config[:Sample]))
        default_point = map(default_level) do s
            string("C", s)            
        end
        pointlevel = default_level
    end
    default_method = if signal_table == SampleDataTable
        Dict(:signal => "area", :rel_sig => "relative_signal", :est_conc => "estimated_concentration", :nom_conc => "nominal_concentration", :acc => "accuracy", :delim => delim, :levelname => "Level")
    else
        Dict(:signal => "area", :rel_sig => "relative_signal", :est_conc => "estimated_concentration", :nom_conc => "nominal_concentration", :acc => "accuracy", :delim => delim, :pointlevel => pointlevel)
    end
    for (k, v) in default_method
        get!(method_config, k, v)
    end
    if haskey(method_config, :pointlevel)
        if length(method_config[:pointlevel]) != length(pointlevel) || any(x -> !in(x, default_level), method_config[:pointlevel])
            @warn "Replace poinlevel with default: $(pointlevel)"
            method_config[:pointlevel] = pointlevel
        end
    end
    if signal_table == SampleDataTable
        confms = string("[signal]\n", method_config[:signal], "\n\n[rel_sig]\n", method_config[:rel_sig], "\n\n[est_conc]\n", method_config[:est_conc], 
                        "\n\n[nom_conc]\n", method_config[:nom_conc], "\n\n[acc]\n", method_config[:acc], "\n\n[delim]\n", escape_string(string(method_config[:delim])), 
                        "\n\n[levelname]\n", method_config[:levelname])
        table[2] = Table(; Symbol(signal_config[:Sample]) => default_point, Symbol(method_config[:levelname]) => default_level, (Symbol.(signal_config[:Analyte]) .=> Ref(repeat([1.0], length(default_point))))...)
    else
        confms = string("[signal]\n", method_config[:signal], "\n\n[rel_sig]\n", method_config[:rel_sig], "\n\n[est_conc]\n", method_config[:est_conc], 
                        "\n\n[nom_conc]\n", method_config[:nom_conc], "\n\n[acc]\n", method_config[:acc], "\n\n[delim]\n", escape_string(string(method_config[:delim])), 
                        "\n\n[pointlevel]\n", join(method_config[:pointlevel], "\n"))
        table[2] = Table(; Symbol(signal_config[:Analyte]) => default_analyte, (Symbol.(signal_config[:Sample]) .=> Ref(repeat([1.0], length(default_analyte))))...)
    end
    if conc_table == SampleDataTable
        table[3] = Table(; Symbol(conc_config[:Sample]) => default_level, (Symbol.(conc_config[:Analyte]) .=> Ref(repeat([1.0], length(default_level))))...)
    else
        table[3] = Table(; Symbol(conc_config[:Analyte]) => default_analyte, (Symbol.(conc_config[:Sample]) .=> Ref(repeat([1.0], length(default_analyte))))...)
    end
    signal = method_config[:signal]
    nom_conc = method_config[:nom_conc]
    data_name = string(0, "_", signal, data_table == SampleDataTable ? ".sdt" : ".adt")
    mkpath(joinpath(file, "data.at", data_name))
    write(joinpath(file, "data.at", data_name, "config.txt"), confs[1])
    CSV.write(joinpath(file, "data.at", data_name, "table.txt"), table[1]; delim = data_config[:delim])
    cp = joinpath(file, "method.am", conc_table == SampleDataTable ? "$nom_conc.sdt" : "$nom_conc.adt")
    mkpath(cp)
    sp = joinpath(file, "method.am", signal_table == SampleDataTable ? "$signal.sdt" : "$signal.adt")
    mkpath(sp)
    write(joinpath(sp, "config.txt"), confs[2])
    CSV.write(joinpath(sp, "table.txt"), table[2]; delim = signal_config[:delim])
    write(joinpath(cp, "config.txt"), confs[3])
    CSV.write(joinpath(cp, "table.txt"), table[3]; delim = conc_config[:delim])
    CSV.write(joinpath(file, "method.am", "analytetable.txt"), Table(; analyte = default_analyte, isd = repeat([0], length(default_analyte)), std = collect(eachindex(default_analyte))); delim = method_config[:delim])
    write(joinpath(file, "method.am", "config.txt"), confms)
    write(joinpath(file, "config.txt"), string("[delim]\n", escape_string(string(delim))))
end

function print_summary(io::IO, cal::ExternalCalibrator{A}) where A
    print(io, "ExternalCalibrator{$A} of ", cal.analyte, " with ", length(unique(cal.table.level[cal.table.include])), " levels and ", length(findall(cal.table.include)), " points")
end

function print_summary(io::IO, cal::InternalCalibrator{A}) where A
    print(io, "InternalCalibrator{$A} of ", cal.analyte, " with single level")
end

function show(io::IO, ::MIME"text/plain", cal::ExternalCalibrator{A, N, T}) where {A, N, T}
    print(io, "ExternalCalibrator{$A, $N, $(shorten_type_repr(T))} with ", length(unique(cal.table.level[cal.table.include])), " levels and ", length(findall(cal.table.include)), " points:\n")
    print(io, "∘ STD: ", cal.analyte, "\n")
    print(io, "∘ ISD: ", cal.isd, "\n")
    print(io, "∘ Model: ")
    print_summary(io, cal.model)
    print(io, "\n")
    print(io, "∘ Equation: ", formula_repr(cal), "\n")
    print(io, "∘ R²: ", r2(cal.machine))
end

function show(io::IO, ::MIME"text/plain", cal::InternalCalibrator{A, N}) where {A, N}
    print(io, "InternalCalibrator{$A, $N} with single level:\n")
    print(io, "∘ STD (ISD): ", cal.analyte, "\n")
    print(io, "∘ Concentration: ", cal.conc)
end

function print_summary(io::IO, model::CalibrationModel{T}) where T
    print(io, "CalibrationModel | ")
    print(io, human_name(T()))
    print(io, " | Weight ", human_name(getweight(model)))
end

function show(io::IO, ::MIME"text/plain", model::CalibrationModel{T}) where T
    print(io, "CalibrationModel\n")
    print(io, "∘ Type: ")
    print(io, human_name(T()))
    print(io, "\n∘ Weight: ", human_name(getweight(model)))
end

function show(io::IO, model::LsqFitMachine) 
    print(io, "LsqFitMachine\n")
    print(io, model.fit)
end

function show(io::IO, ::MIME"text/plain", model::LsqFitMachine) 
    print(io, "LsqFitMachine\n")
    show(io, MIME"text/plain"(), model.fit)
end


function print_summary(io::IO, batch::Batch{A, T}) where {A, T}
    print(io, "Batch{$A} with ", count(<(0), batch.method.analytetable.isd), " internal standards out of ", length(batch.method.analytetable.analyte), " analytes")
end

function show(io::IO, ::MIME"text/plain", batch::Batch)
    print_summary(io, batch)
    print(io, ":\n∘ Analytes:\n")
    show(io, MIME"text/plain"(), batch.method.analytetable)
    print(io, "\n\n∘ Calibrators:")
    if length(batch.calibrator) > 10
        for c in @view batch.calibrator[1:5]
            print(io, "\n ")
            print_summary(io, c)
        end
        print(io, "\n ⋮")
        for c in @view batch.calibrator[end - 4:end]
            print(io, "\n ")
            print_summary(io, c)
        end
    else
        for c in batch.calibrator
            print(io, "\n ")
            print_summary(io, c)
        end
    end
    print(io, "\n\n∘ Data:\n")
    show(io, MIME"text/plain"(), batch.data)
end

function shorten_type_repr(T)
    t = repr(T)
    length(t) > 50 ? replace(t, r"\{.*\}" => "{…}") : t
end

function print_summary(io::IO, tbl::SampleDataTable{A, S, N, T}) where {A, S, N, T}
    print(io, "SampleDataTable{$A, $S, $N, $(shorten_type_repr(T))} with ", length(analyteobj(tbl)), " analytes and ", length(sampleobj(tbl)), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::SampleDataTable)
    print_summary(io, tbl)
    println(io, ":")
    show(io, MIME"text/plain"(), table(tbl))
end

function print_summary(io::IO, tbl::AnalyteDataTable{A, S, N, T}) where {A, S, N, T}
    print(io, "AnalyteDataTable{$A, $S, $N, $(shorten_type_repr(T))} with ", length(analyteobj(tbl)), " analytes and ", length(sampleobj(tbl)), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::AnalyteDataTable)
    print_summary(io, tbl)
    println(io, ":")
    show(io, MIME"text/plain"(), table(tbl))
end

function print_summary(io::IO, tbl::AnalysisTable{A, S, T}) where {A, S, T}
    print(io, "AnalysisTable{$A, $S, $(shorten_type_repr(T))} with ", length(analyteobj(tbl)), " analytes and ", length(sampleobj(tbl)), " samples")
end

function show(io::IO, ::MIME"text/plain", tbl::AnalysisTable)
    print_summary(io, tbl)
    print(io, ":")
    if length(tbl) > 2
        ks = propertynames(tbl)
        print(io, "\n∘ ", ks[1], " | ")
        show(io, MIME"text/plain"(), tbl[ks[1]])
        print(io, "\n   ⋮")
        print(io, "\n∘ ", ks[end], " | ")
        show(io, MIME"text/plain"(), tbl[ks[end]])
    else
        for (k, v) in pairs(tbl)
            print(io, "\n∘ ", k, " | ")
            show(io, MIME"text/plain"(), v)
        end
    end

end

function print_summary(io::IO, tbl::AnalysisMethod{A}) where A
    if isnothing(tbl.signaltable)
        print(io, "AnalysisMethod{$A} with ", count(<(0), tbl.analytetable.isd), " internal standards out of ", length(tbl.analytetable.analyte), " analytes")
    else
        print(io, "AnalysisMethod{$A} with ", count(<(0), tbl.analytetable.isd), " internal standards out of ", length(tbl.analytetable.analyte), " analytes, ", length(sampleobj(tbl.conctable)), " levels and ", length(sampleobj(tbl.signaltable)), " points")
    end
end

function show(io::IO, ::MIME"text/plain", tbl::AnalysisMethod)
    print_summary(io, tbl)
    print(io, ":\n∘ Analytes:\n")
    show(io, MIME"text/plain"(), tbl.analytetable)
    print(io, "\n\n∘ Level: ", join(tbl.pointlevel, ", "), "\n")
    print(io, "∘ Concentration: \n")
    show(io, MIME"text/plain"(), table(tbl.conctable))
    print(io, "\n\n∘ Signal: \n")
    show(io, MIME"text/plain"(), isnothing(tbl.signaltable) ? nothing : table(tbl.signaltable))
    print(io, "\n\n∘ Data keys: \n")
    print(io,  join([tbl.signal, tbl.rel_sig, tbl.est_conc, tbl.nom_conc, tbl.acc], ", "), "\n")
end

"""
    human_name(::CurveType)
    human_name(::ComposedWeight)

Names for human.
"""
human_name(::T) where {T <: CurveType} = repr(T)
human_name(::Proportional) = "Proportional"
human_name(::Linear) = "Linear"
human_name(::QuadraticProportional) = "QuadraticProportional"
human_name(::Quadratic) = "Quadratic"
human_name(::Logarithmic) = "Logarithmic"
human_name(::Exponential) = "Exponential"
human_name(::Power) = "Power"

human_name(x::T) where {T <: ComposedWeight} = repr(x)
human_name(::ConstWeight) = "1"
human_name(::RootXWeight) = "1/√x"
human_name(::RootYWeight) = "1/√y"
human_name(::RootXYWeight) = "1/√(x+y)"
human_name(::RootLogXWeight) = "1/√ln(x)"
human_name(::RootLogYWeight) = "1/√ln(y)"
human_name(::RootExpXWeight) = "1/√eˣ"
human_name(::RootExpYWeight) = "1/√eʸ"
human_name(::XWeight) = "1/x"
human_name(::YWeight) = "1/y"
human_name(::XYWeight) = "1/(x+y)"
human_name(::LogXWeight) = "1/ln(x)"
human_name(::LogYWeight) = "1/ln(y)"
human_name(::ExpXWeight) = "1/eˣ"
human_name(::ExpYWeight) = "1/eʸ"
human_name(::SqXWeight) = "1/x²"
human_name(::SqYWeight) = "1/y²"
human_name(::SqLogXWeight) = "1/ln(x)²"
human_name(::SqLogYWeight) = "1/ln(y)²"
human_name(::SqExpXWeight) = "1/e²ˣ"
human_name(::SqExpYWeight) = "1/e²ʸ"

"""
    human_name_ascii(::CurveType)
    human_name_ascii(::ComposedWeight)

Names for text output.
"""
human_name_ascii(x::T) where {T <: CurveType} = human_name(x)
human_name_ascii(x::T) where {T <: ComposedWeight} = human_name(x)
human_name_ascii(::RootXWeight) = "1/x^0.5"
human_name_ascii(::RootYWeight) = "1/y^0.5"
human_name_ascii(::RootXYWeight) = "1/(x+y)^0.5"
human_name_ascii(::RootLogXWeight) = "1/ln(x)^0.5"
human_name_ascii(::RootLogYWeight) = "1/ln(y)^0.5"
human_name_ascii(::RootExpXWeight) = "1/e^(x/2)"
human_name_ascii(::RootExpYWeight) = "1/e^(y/2)"
human_name_ascii(::ExpXWeight) = "1/e^x"
human_name_ascii(::ExpYWeight) = "1/e^y"
human_name_ascii(::SqXWeight) = "1/x^2"
human_name_ascii(::SqYWeight) = "1/y^2"
human_name_ascii(::SqLogXWeight) = "1/ln(x)^2"
human_name_ascii(::SqLogYWeight) = "1/ln(y)^2"
human_name_ascii(::SqExpXWeight) = "1/e^(2x)"
human_name_ascii(::SqExpYWeight) = "1/e^(2y)"

"""
    encode_config(model::CalibrationModel)

Encode keyword arguments of `mkcalmodel` for specific calibration model type into config file.
"""
function encode_config(model::CalibrationModel{T}) where T 
    ("\n\n[model]\n", string(model))
end

"""
    formula_repr(cal::InternalCalibrator; digits = nothing, sigdigits = 4)
    formula_repr(cal::ExternalCalibrator; digits = nothing, sigdigits = 4)

Return string representation of formula of `cal` with specified `digits` and `sigdigits`. See `format_number`.
"""
formula_repr(cal::InternalCalibrator; digits = nothing, sigdigits = 4) = "y = $(format_number(1/cal.conc; digits, sigdigits))x"
formula_repr(cal::ExternalCalibrator; digits = nothing, sigdigits = 4) = formula_repr(cal.model, cal.machine; digits, sigdigits)
function formula_repr(model::CalibrationModel, machine; digits = nothing, sigdigits = 4)
    β = coef(machine)
    formula_repr(model, β; digits, sigdigits)
end

formula_repr(::CalibrationModel{Proportional}, β::AbstractVector; digits = nothing, sigdigits = 4) = "y = $(format_number(β[1]; digits, sigdigits))x"
function formula_repr(::CalibrationModel{Linear}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "x")
end
function formula_repr(::CalibrationModel{QuadraticProportional}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), "x", op[1], abs(format_number(β[2]; digits, sigdigits)), "x²")
end
function formula_repr(::CalibrationModel{Quadratic}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "x", op[2], abs(format_number(β[3]; digits, sigdigits)), "x²")
end
function formula_repr(::CalibrationModel{Logarithmic}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " - " : " + "
    end
    string("y = ", format_number(β[1]; digits, sigdigits), op[1], abs(format_number(β[2]; digits, sigdigits)), "log(x)")
end
function formula_repr(::CalibrationModel{Exponential}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? " -" : ""
    end
    string("y = ", format_number(β[1]; digits, sigdigits), "e ^ (", op[1], abs(format_number(β[2]; digits, sigdigits)), "x)")
end
function formula_repr(::CalibrationModel{Power}, β::AbstractVector; digits = nothing, sigdigits = 4) 
    op = map(β[2:end]) do b
        b < 0 ? "-" : "" 
    end
    string("y = ", format_number(β[1]; digits, sigdigits), "x ^ ", string(isempty(op[1]) ? "" : "(", op[1], abs(format_number(β[2]; digits, sigdigits)), isempty(op[1]) ? "" : ")"))
end
"""
    formula_repr_ascii(cal::AbstractCalibrator; digits = nothing, sigdigits = 4)

Return string representation of formula of `cal` for text file output.
"""
formula_repr_ascii(cal::AbstractCalibrator; digits = nothing, sigdigits = 4) = replace(formula_repr(cal; digits, sigdigits), "x²" => "x^2")

"""
    format_number(x; digits = nothing, sigdigits = 4)

Return string representation of number with specified `digits`.

If `digits` is `nothing`, the function uses `sigdigits` instead.
"""
format_number(x; digits = nothing, sigdigits = 4) = isnothing(digits) ? format_number2int(round(x; sigdigits)) : format_number2int(round(x; digits))
format_number2int(x) = x == round(x) ? round(Int, x) : x
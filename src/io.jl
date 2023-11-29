
function read_calibration(file::String, T, source::AbstractAnalysisTable, )
    config_kw, config_vl = eachline(joinpath(file, "config.csv"))
    config = Dict{Symbol, Any}()
    for (k, v) in zip(split(config_kw, ","), split(config_vl, ","))
        if k == "type" || k == "zero"
            v = v == "TRUE" || v == "True" || v == "true" || v == ""
        elseif k == "weight"
            v = parse(Float64, v)
        elseif k == "analyte" || k == "isd"
            v = parse(Int, v)
        elseif k == "beta"
            v = parse(Float64, v)
        else
            continue
        end
        get!(config, Symbol(k), v)
    end
    endswith(file, ".scal") && return SingleCalibration(config[:analyte], config[:isd], source, config[:beta])
    tbl = CSV.read(joinpath(file, "calibration.csv"), T)
    calibration(tbl, source; config...)
end

function read_analysistable(file::String, T; anayte_fn = identity)
    tbl = CSV.read(joinpath(file, "calibration.csv"), T)
    config = filter!(!isempty, readlines(joinpath(file, "config.txt")))
    config = replace.(config, Ref(" " => ""))
    if pop!(config) == "R"
        sample_name = Symbol(pop!(config))
        analyte_name = Symbol.(config)
        RowAnalysisTable(sample_name, analyte_name, analyte_fn.(analyte_name), tbl)
    else
        sample_name, analyte_name... = config
        analyte_name = Symbol.(analyte_name)
        ColumnAnalysisTable(Symbol(sample_name), analyte_name, analyte_fn.(analyte_name), tbl)
    end
end

function read_project(file::String, T; analyte_fn = identity)
    source = read_analysistable(file, T; analyte_fn)
    fs = findfirst(f -> endswith(f, ".sample"), readdir(file))
    Project([read_calibration(joinpath(file, f), T, source) for f in readdir(file) if endswith(f, ".cal") || endswith(f, ".scal")],
        isnothing(fs) ? nothing : read_analysistable(joinpath(file, fs), T; analyte_fn)
    )
end

function show(io::IO, cal::MultipleCalibration)
    print(io, "Calibration of $(cal.source.analytes[cal.analyte]) with multiple levels")
end

function show(io::IO, cal::SingleCalibration)
    print(io, "Calibration of $(cal.source.analytes[cal.analyte]) with single level")
end

function show(io::IO, ::MIME"text/plain", cal::MultipleCalibration),
    print(io, cal, ":\n")
    print(io, "∘ ISD: ", cal.isd > 0 ? cal.source.analytes[cal.isd] : nothing, "\n"),
    print(io, "∘ Type: ", cal.type ? "quadratic" : "linear")
    print(io, "∘ (0, 0): ", cal.zero ? "included\n" : "ommited\n")
    print(io, "∘ Weight: ", weight_repr(cal.weight), "\n")
    print(io, "∘ Formula", formula_repr(cal))
end

function show(io::IO, ::MIME"text/plain", cal::SingleCalibration)
    print(io, cal, ":\n")
    print(io, "∘ ISD: ", cal.isd > 0 ? cal.source.analytes[cal.isd] : nothing, "\n")
    print(io, "∘ Formula", formula_repr(cal))
end

function show(io::io::IO, project::Project)
    string("Project with $(length(project.calibration)) analytes")
end
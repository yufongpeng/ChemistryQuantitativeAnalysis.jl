struct AnalyteN
    N::Int
end
show(io::IO, analyte::AnalyteN) = print(io, "Analyte", analyte.N)
function AnalyteN(name::AbstractString) 
    nm = match(r"Analyte(\d*)$", name)
    isnothing(nm) && throw(ArgumentError("Invalid name"))
    (n, ) = nm 
    AnalyteN(parse(Int, n))
end

struct AnalyteG1
    name::String
end
struct AnalyteG2
    name::String
end
struct AnalyteOther
    name::String
end
const AnalyteTest = Union{AnalyteG1, AnalyteG2, AnalyteOther}
show(io::IO, analyte::AnalyteTest) = print(io, analyte.name)
function (::Type{AnalyteTest})(name::String)
    g = match(r"^G(\d)\(.*\)$", name)
    isnothing(g) && return AnalyteOther(name)
    g = parse(Int, first(g))
    g == 1 ? AnalyteG1(name) : g == 2 ? AnalyteG2(name) : AnalyteOther(name)
end

analytes = typedmap(AnalyteTest, ["G1(drug_a)", "G2(drug_a)", "G1(drug_b)", "G2(drug_b)"])
conc = Float64[0, 1, 2, 5, 10, 20, 50, 100]
signal1 = vcat(0.001, Float64[1, 2, 5, 10, 20, 50, 100], 0.005, [1, 2, 5, 10, 20, 50, 100] .+ 0.1, 0.002, [1, 2, 5, 10, 20, 50, 100] .- 0.1)
signal2 = vcat(0.002, Float64[1, 2, 5, 10, 20, 50, 100] .^ 2, 0.001, [1, 2, 5, 10, 20, 50, 100] .^ 2 .+ 0.1, 0.005, [1, 2, 5, 10, 20, 50, 100] .^ 2 .- 0.1)

conctable = SampleDataTable(
    Table(;
        var"level" = collect(0:7), 
        var"G1(drug_a)" = conc,
        var"G1(drug_b)" = conc .* 10, 
        var"G2(drug_a)" = repeat([50.0], 8),
        var"G2(drug_b)" = repeat([50.0], 8)),
    AnalyteTest,
    :level
)
signaltable = SampleDataTable(
    Table(;
        var"point" = reshape([string(a, "_", b) for (a, b) in Iterators.product(0:7, 1:3)], 24), 
        var"level" = repeat(0:7, 3),
        var"G1(drug_a)" = signal1,
        var"G2(drug_a)" = repeat([5.0], 24),
        var"G1(drug_b)" = signal2,
        var"G2(drug_b)" = repeat([2.0], 24)), 
    AnalyteTest,
    :point; 
    analytename = Symbol.(analytes)
)
method = AnalysisMethod(conctable, signaltable, :area, :level; analyte = analytes, isd = [2, -1, 4, -1], std = [1, -1, 3, -1])
cdata = AnalysisTable([:area], [
    SampleDataTable(
        Table(;
            var"Sample" = ["S1", "S2", "S3"], 
            var"G1(drug_a)" = Float64[6, 24, 54],
            var"G2(drug_a)" = Float64[5, 6, 6],
            var"G1(drug_b)" = Float64[200, 800, 9800],
            var"G2(drug_b)" = Float64[2, 2, 2]), 
        AnalyteTest,
        :Sample
        )
    ]
)
rdata = analysistable((:area => 
    AnalyteDataTable(
        Table(;
            var"Analyte" = analytes, 
            var"S1" = Float64[6, 6, 200, 2],
            var"S2" = Float64[24, 6, 800, 2],
            var"S3" = Float64[54, 6, 9800, 2]
            ), 
        :Analyte
        )
    , )
)
cbatch = Batch(method, cdata)
rbatch = Batch(method, rdata)
ebatch = Batch(method)
sdt = SampleDataTable(AnalyteN, :Sample, CSV.read(joinpath(datapath, "area.sdt", "table.txt"), Table; delim = '\t', stringtype = String))
adt = AnalyteDataTable(CSV.read(joinpath(datapath, "area.adt", "table.txt"), Table; delim = '\t', stringtype = String), :Analyte)
sbatch = Batch(sdt; conc_factor = repeat([1], 11))
abatch = Batch(adt; level = vcat(repeat(1:6; inner = 3), repeat([missing], 50)), ratio = [0.1, 0.2, 0.5, 1, 2, 5], conc_factor = repeat([1], 11))
sdt2 = CQA.read(joinpath(datapath, "area2.sdt"), Table)
adt2 = CQA.read(joinpath(datapath, "area2.adt"), Table)
sbatch2 = Batch(sdt2; level = :Level, ratio = :Ratio, conc_factor = 1:10)
abatch2 = Batch(adt2; level = vcat(repeat(1:6; inner = 3), repeat([missing], 50)), ratio = [0.1, 0.2, 0.5, 1, 2, 5], conc_factor = :Factor)
ical = CQA.analyze_calibrator!(InternalCalibrator(analyteobj(method.conctable)[3], 50.0))

initial_mc_c = CQA.read(joinpath(datapath, "initial_mc_c.batch"), DataFrame)
initial_mc_r = CQA.read(joinpath(datapath, "initial_mc_r.batch"), DataFrame)
initial_sc_c = CQA.read(joinpath(datapath, "initial_sc_c.batch"), DataFrame)
initial_sc_r = CQA.read(joinpath(datapath, "initial_sc_r.batch"), DataFrame)
save_mc_c = CQA.read(joinpath(datapath, "save_mc_c.batch"), DataFrame)
save_mc_r = CQA.read(joinpath(datapath, "save_mc_r.batch"), DataFrame)
save_sc_c = CQA.read(joinpath(datapath, "save_sc_c.batch"), DataFrame)
save_sc_r = CQA.read(joinpath(datapath, "save_sc_r.batch"), DataFrame)
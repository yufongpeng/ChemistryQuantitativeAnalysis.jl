using ChemistryQuantitativeAnalysis, TypedTables, DataFrames, Dictionaries
using Test, JET
import Base: show, convert
const CQA = ChemistryQuantitativeAnalysis

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
conc = Float64[1, 2, 5, 10, 20, 50, 100]
signal1 = vcat(Float64[1, 2, 5, 10, 20, 50, 100], [1, 2, 5, 10, 20, 50, 100] .+ 0.1, [1, 2, 5, 10, 20, 50, 100] .- 0.1)
signal2 = vcat(Float64[1, 2, 5, 10, 20, 50, 100] .^ 2, [1, 2, 5, 10, 20, 50, 100] .^ 2 .+ 0.1, [1, 2, 5, 10, 20, 50, 100] .^ 2 .- 0.1)

const datapath = joinpath(@__DIR__(), "data")
isapprox_nan(x::Float64, y::Float64) = isnan(x) && isnan(y) ? true : isapprox(x, y)
isapprox_nan(x::T, y::T) where T = true
test_show(x) = show(IOBuffer(), MIME"text/plain"(), x)
macro test_error(err, x)
    return quote
        try 
            $x
            false
        catch e
            isa(e, $err) ? true : false
        end
    end
end

macro test_error(x)
    return quote
        try 
            $x
            false
        catch e
            true
        end
    end
end

macro test_noerror(x)
    return quote
        try 
            $x
            true
        catch e
            false
        end
    end
end

@testset "ChemistryQuantitativeAnalysis.jl" begin
    @testset "Constructors" begin
        global conctable = ColumnDataTable(
            Table(;
                var"level" = collect(1:7), 
                var"G1(drug_a)" = conc,
                var"G1(drug_b)" = conc .* 10), 
            AnalyteTest,
            :level
        )
        global signaltable = ColumnDataTable(
            Table(;
                var"point" = reshape([string(a, "_", b) for (a, b) in Iterators.product(1:7, 1:3)], 21), 
                var"level" = repeat(1:7, 3),
                var"G1(drug_a)" = signal1,
                var"G2(drug_a)" = repeat([5.0], 21),
                var"G1(drug_b)" = signal2,
                var"G2(drug_b)" = repeat([2.0], 21)), 
            AnalyteTest,
            :point; 
            analytename = Symbol.(analytes)
        )
        global method = AnalysisMethod(conctable, signaltable, :area, :level; analyte = analytes, isd = [2, -1, 4, -1], calibration = [1, -1, 3, -1])
        global cdata = AnalysisTable([:area], [
            ColumnDataTable(
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
        global rdata = AnalysisTable([:area], [
            RowDataTable(
                Table(;
                    var"Analyte" = analytes, 
                    var"S1" = Float64[6, 6, 200, 2],
                    var"S2" = Float64[24, 6, 800, 2],
                    var"S3" = Float64[54, 6, 9800, 2]
                    ), 
                :Analyte
                )
            ]
        )
        global cbatch = Batch(method, cdata)
        global rbatch = Batch(method, rdata)
        global ebatch = Batch(method)
        @test getanalyte(cdata.area, 1) == getanalyte(rdata.area, 1) == getanalyte(ColumnDataTable(rdata.area, :Sample, Table), 1)
        @test getsample(cdata.area, 2) == getsample(rdata.area, 2) == getsample(RowDataTable(cdata.area, :Analyte, Table), 2)
        @test propertynames(cbatch) == (:method, :calibration, :data, :analyte, :isd, :nonisd, :point, :level)
        @test propertynames(cbatch.method) == (:analytetable, :signal, :pointlevel, :conctable, :signaltable, :analyte, :isd, :nonisd, :point, :level)
        @test propertynames(cdata) == (:area, )
        @test propertynames(cdata.area) == (:Sample, Symbol.(analytes)...)
        @test propertynames(rdata.area) == (:Analyte, Symbol.(["S1", "S2", "S3"])...)
    end
    @testset "Interface.jl" begin
        cdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 6
        @test cdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 6
        rdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 4
        @test rdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 4
        rdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 6
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) collect(cdata.area)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) collect(rdata.area)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) columns(cdata.area)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) columns(rdata.area)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) rows(cdata.area)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) rows(rdata.area)
        @test collect(eachanalyte(cdata.area)) == collect(eachanalyte(rdata.area))
        @test collect(eachsample(cdata.area)) == collect(eachsample(rdata.area))
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) insert!(cbatch.data, :true_concentration, deepcopy(cbatch.data.area))
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) set!(cbatch.data, :true_concentration, deepcopy(cbatch.data.area))
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) delete!(cbatch.data, :true_concentration)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) unset!(cbatch.data, :true_concentration)
    end
    @testset "Calibration" begin
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) accuracy(cbatch.calibration[1])
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) quantification(rbatch.calibration[2])
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) CQA.calibrate!(rbatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) CQA.validate!(cbatch)
        cbatch.calibration[2].type = false
        rbatch.calibration[AnalyteG1("G1(drug_b)")].type = false
        update_calibration!(cbatch, AnalyteG1("G1(drug_b)"))
        update_calibration!(rbatch, 2)
        @test all(isapprox.(cbatch.calibration[1].table.accuracy[1:3], [1, 1.1, 0.9]))
        # @test all(isapprox.(cbatch.calibration[2].table.accuracy[1:3], [1, sqrt(1.1), sqrt(0.9)]))
        # @test all(isapprox.(rbatch.calibration[1].table.accuracy[1:3], [1, 1.1, 0.9]))
        @test all(isapprox.(rbatch.calibration[2].table.accuracy[1:3], [1, sqrt(1.1), sqrt(0.9)]))
        @test isapprox(update_calibration!(SingleCalibration((analyteobj(method.conctable)[2], ), 100.0), method).conc, first(getanalyte(method.conctable, 2)))
    end
    @testset "Quantification" begin
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) set_relative_signal(rbatch.data, rbatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) set_inv_predict(cbatch.data, cbatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) set_quantification(rbatch.data, rbatch)
        # @test_call target_modules = (ChemistryQuantitativeAnalysis, ) update_quantification!(rbatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) set_quantification!(cbatch.data, cbatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) quantification(cbatch, cbatch.data)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) quantification(cbatch.calibration[1], cbatch.data.area, cbatch.calibration[1].analyte...)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) update_quantification!(ebatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) update_relative_signal!(ebatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) update_inv_predict!(ebatch)
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) update_accuracy!(ebatch)
        @test all(isapprox.(update_quantification!(cbatch).data.estimated_concentration.var"G1(drug_b)", getanalyte(update_quantification!(cbatch).data.estimated_concentration, AnalyteTest("G1(drug_b)"))))
        # @test all(isapprox_nan.(set_quantification!(rbatch.data, rbatch).estimated_concentration.S1, set_inv_predict!(set_relative_signal(rbatch.data, rbatch), rbatch).estimated_concentration.S1))
        # @test all(isapprox.(quantification(cbatch, cbatch.data).var"G1(drug_a)", inv_predict(cbatch.calibration[1], relative_signal(cbatch, cbatch.data))))
        # @test all(isapprox_nan.(quantification(rbatch, rbatch.data).S2, update_quantification!(rbatch).data.estimated_concentration.S2))
        # @test all(isapprox.(update_inv_predict!(update_relative_signal!(cbatch)).data.estimated_concentration.var"G1(drug_b)", set_quantification(cbatch.data, cbatch).estimated_concentration.var"G1(drug_b)"))
        # @test all(isapprox.(set_inv_predict(cbatch.data, cbatch).estimated_concentration.var"G1(drug_b)", cbatch.data.estimated_concentration.var"G1(drug_b)"))
        set!(cbatch.data, :true_concentration, deepcopy(cbatch.data.estimated_concentration))
        @test all(isapprox.(update_accuracy!(cbatch).data.accuracy.var"G1(drug_b)", [1.0, 1.0, 1.0]))
        @test all(isapprox.(set_accuracy!(cbatch.data).accuracy.var"G1(drug_b)", [1.0, 1.0, 1.0]))
    end
    @testset "Utils" begin
        @test getanalyte(cdata.area, AnalyteG1("G1(drug_b)")) == getanalyte(cdata.area, Symbol("G1(drug_b)"))
        @test getanalyte(rdata.area, AnalyteG1("G1(drug_b)")) == getanalyte(rdata.area, Symbol("G1(drug_b)"))
        @test getsample(cdata.area, "S2") == getsample(cdata.area, Symbol("S2"))
        @test getsample(rdata.area, Symbol("S2")) == getproperty(CQA.table(rdata.area), Symbol("S2"))
        @test all(isapprox.(dynamic_range(cbatch.calibration[1]), (1, 100)))
        @test_call target_modules = (ChemistryQuantitativeAnalysis, ) signal_range(rbatch.calibration[2])
        @test endswith(formula_repr_utf8(cbatch.calibration[2]), "x^2")
        @test weight_repr_utf8(cbatch.calibration[1]) == "none"
        @test weight_value("1/x^3") == -3
    end
    @testset "IO" begin
        global initial_mc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_mc_c.batch"), DataFrame)
        global initial_mc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_mc_r.batch"), DataFrame)
        global initial_sc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_sc_c.batch"), DataFrame)
        global initial_sc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_sc_r.batch"), DataFrame)
        ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_mc_c.batch", "data.at"), DataFrame)
        ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_mc_c.batch", "data.at", "0_area.dt"), DataFrame)
        ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_mc_c.batch", "method.mt"), DataFrame)
        ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_mc_c.batch", "calibration", "1.mcal"), DataFrame)
        ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_sc_c.batch", "calibration", "1.scal"), DataFrame)
        update_quantification!(initial_mc_c)
        update_quantification!(initial_mc_r)
        update_quantification!(initial_sc_c)
        update_quantification!(initial_sc_r)
        global save_mc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_mc_c.batch"), DataFrame)
        global save_mc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_mc_r.batch"), DataFrame)
        global save_sc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_sc_c.batch"), DataFrame)
        global save_sc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_sc_r.batch"), DataFrame)
        @test all(isapprox_nan.(collect(CQA.table(initial_mc_c.data.estimated_concentration)[1, :]), collect(CQA.table(save_mc_c.data.estimated_concentration)[1, :])))
        @test all(isapprox_nan.(collect(CQA.table(initial_mc_r.data.estimated_concentration)[1, :]), collect(CQA.table(save_mc_r.data.estimated_concentration)[1, :])))
        @test all(isapprox_nan.(collect(CQA.table(initial_sc_c.data.estimated_concentration)[1, :]), collect(CQA.table(save_sc_c.data.estimated_concentration)[1, :])))
        @test all(isapprox_nan.(collect(CQA.table(initial_sc_r.data.estimated_concentration)[1, :]), collect(CQA.table(save_sc_r.data.estimated_concentration)[1, :])))
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_mc_c.batch"), initial_mc_c)
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_mc_r.batch"), initial_mc_r)
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_sc_c.batch"), initial_sc_c)
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_sc_r.batch"), initial_sc_r)
        @test @test_noerror test_show(save_mc_c)
        # @test @test_noerror test_show(save_mc_r)
        # @test @test_noerror test_show(save_sc_c)
        @test @test_noerror test_show(save_sc_r)
        @test @test_noerror test_show(save_mc_c.calibration[1])
        # @test @test_noerror test_show(save_mc_r.calibration[1])
        # @test @test_noerror test_show(save_sc_c.calibration[1])
        @test @test_noerror test_show(save_sc_r.calibration[1])
        @test @test_noerror test_show(save_mc_c.method)
        # @test @test_noerror test_show(save_mc_r.method)
        # @test @test_noerror test_show(save_sc_c.method)
        # @test @test_noerror test_show(save_sc_r.method)
        @test @test_noerror test_show(save_mc_c.data)
        @test @test_noerror test_show(save_mc_r.data)
        # @test @test_noerror test_show(save_sc_c.data)
        # @test @test_noerror test_show(save_sc_r.data)
        mkbatch(joinpath(datapath, "new1.batch"); data_config = Dict(:Type => "R", :Sample => ["S0", "S1", "S2"]), signal_config = Dict(:Type => "R", :Sample => ["C0"]), conc_config = Dict(:Type => "R", :Sample => [1]))
        mkbatch(joinpath(datapath, "new2.batch"); data_config = Dict(:Analyte => ["A1", "A2"]), signal_config = Dict(:Analyte => ["A1", "A2", "A3"]), conc_config = Dict(:Type => "R", :Sample => [1, 2, 3, 4, 5, 6]))
        mkbatch(joinpath(datapath, "new3.batch"); data_config = Dict(:Type => "R", :Sample => ["S0", "S1", "S2"]), conc_config = Dict(:Analyte => ["A1", "A2", "A3"]), signal_config = Dict(:Type => "R", :Sample => ["C1", "C2", "C3", "C4", "C5", "C6"]))
        mkbatch(joinpath(datapath, "new4.batch"); method_config = Dict(:signal => "height", :levelname => "L"))
        n1 = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "new1.batch"), DataFrame)
        n2 = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "new2.batch"), DataFrame)
        n3 = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "new3.batch"), DataFrame)
        n4 = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "new4.batch"), DataFrame)
        @test sampleobj(n1.data) == ["S0", "S1", "S2"]
        @test length(unique(n2.calibration[1].table.level[n2.calibration[1].table.include])) == 6
        @test analyteobj(n3.data) == ["A1", "A2", "A3"]
        @test :L in propertynames(n4.method.signaltable)
    end
end
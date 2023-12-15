using ChemistryQuantitativeAnalysis, TypedTables, DataFrames
using Test
import Base: show

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
function Union{AnalyteG1, AnalyteG2, AnalyteOther}(name::String)
    g = match(r"^G(\d)\(.*\)$", name)
    isnothing(g) && return AnalyteOther(name)
    g = parse(Int, first(g))
    g == 1 ? AnalyteG1(name) : g == 2 ? AnalyteG2(name) : AnalyteOther(name)
end
analyte_names = ["G1(drug_a)", "G2(drug_a)", "G1(drug_b)", "G2(drug_b)"]
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
            DataFrame(
                "level" => collect(1:7), 
                "G1(drug_a)" => conc,
                "G1(drug_b)" => conc .* 10), 
            :level; 
            analytetype = AnalyteTest
        )
        global signaltable = ColumnDataTable(
            DataFrame(
                "point" => repeat(1:7, 3), 
                "G1(drug_a)" => signal1,
                "G2(drug_a)" => repeat([5.0], 21),
                "G1(drug_b)" => signal2,
                "G2(drug_b)" => repeat([2.0], 21)), 
            :point; 
            analytetype = AnalyteTest
        )
        global method = MethodTable(conctable, signaltable, :area, repeat(1:7, 3); analyte = AnalyteTest.(analyte_names), isd = [2, -1, 4, -1], calibration = [1, -1, 3, -1])
        global rdata = AnalysisTable([:area], [
            RowDataTable(
                DataFrame(
                    "Analyte" => ["G1(drug_a)", "G2(drug_a)", "G1(drug_b)", "G2(drug_b)"], 
                    "S1" => Float64[6, 6, 200, 2],
                    "S2" => Float64[24, 6, 800, 2],
                    "S3" => Float64[54, 6, 9800, 2]
                    ), 
                :Analyte; 
                analytetype = AnalyteTest
                )
            ]
        )
        global cdata = AnalysisTable([:area], [
            ColumnDataTable(
                DataFrame(
                    "Sample" => ["S1", "S2", "S3"], 
                    "G1(drug_a)" => Float64[6, 24, 54],
                    "G2(drug_a)" => Float64[5, 6, 6],
                    "G1(drug_b)" => Float64[200, 800, 9800],
                    "G2(drug_b)" => Float64[2, 2, 2]), 
                :Sample; 
                analytetype = AnalyteTest
                )
            ]
        )
        global rdata = AnalysisTable([:area], [
            RowDataTable(
                DataFrame(
                    "Analyte" => ["G1(drug_a)", "G2(drug_a)", "G1(drug_b)", "G2(drug_b)"], 
                    "S1" => Float64[6, 6, 200, 2],
                    "S2" => Float64[24, 6, 800, 2],
                    "S3" => Float64[54, 6, 9800, 2]
                    ), 
                :Analyte; 
                analytetype = AnalyteTest
                )
            ]
        )
        global cbatch = Batch(method, cdata)
        global rbatch = Batch(method, rdata)
        @test getanalyte(cdata.area, 1) == getanalyte(rdata.area, 1)
        @test getanalyte(cdata.area, 1) == getanalyte(ColumnDataTable(rdata.area, :Sample), 1)
        @test getsample(cdata.area, 2) == getsample(rdata.area, 2)
        @test getsample(RowDataTable(cdata.area, :Analyte), 2) == getsample(rdata.area, 2)
    end
    @testset "Interface.jl" begin
        cdata2 = AnalysisTable([:area], [
            ColumnDataTable(
                Table(cdata.area.table), 
                :Sample; 
                analytetype = AnalyteTest
                )
            ]
        )
        rdata2 = AnalysisTable([:area], [
            RowDataTable(
                Table(rdata.area.table), 
                :Analyte; 
                analytetype = AnalyteTest
                )
            ]
        )
        cdata.area[1, "G2(drug_a)"] = 6
        @test cdata.area[1, "G1(drug_a)"] == 6
        @test collect(cdata2.area) == collect(cdata2.area.table)
        @test collect(rdata2.area) == collect(rdata2.area.table)
        @test columns(cdata2.area) == columns(cdata2.area.table)
        @test columns(rdata2.area) == columns(rdata2.area.table)
        @test rows(cdata2.area) == rows(cdata2.area.table)
        @test rows(rdata2.area) == rows(rdata2.area.table)
        @test collect(eachanalyte(cdata.area)) == collect(eachanalyte(rdata.area))
        @test collect(eachsample(cdata.area)) == collect(eachsample(rdata.area))
    end
    @testset "Calibration" begin
        cbatch.calibration[2].type = false
        rbatch.calibration[2].type = false
        update_calibration!(cbatch, AnalyteG1("G1(drug_b)"))
        update_calibration!(rbatch, AnalyteG1("G1(drug_b)"))
        @test all(isapprox.(cbatch.calibration[1].table.accuracy[1:3], [1, 1.1, 0.9]))
        @test all(isapprox.(cbatch.calibration[2].table.accuracy[1:3], [1, sqrt(1.1), sqrt(0.9)]))
        @test all(isapprox.(rbatch.calibration[1].table.accuracy[1:3], [1, 1.1, 0.9]))
        @test all(isapprox.(rbatch.calibration[2].table.accuracy[1:3], [1, sqrt(1.1), sqrt(0.9)]))
    end
    @testset "Quantification" begin
        @test all(isapprox_nan.(set_quantification!(rbatch.data, rbatch).estimated_concentration.S1, update_inv_predict!(rbatch, set_relative_signal(rbatch.data, rbatch)).data.estimated_concentration.S1))
        @test all(isapprox.(quantification(cbatch, cbatch.data).var"G1(drug_a)", inv_predict(cbatch.calibration[1], relative_signal(cbatch, cbatch.data))))
        @test all(isapprox_nan.(quantification(rbatch, rbatch.data).S2, update_quantification!(rbatch, rbatch.data).data.estimated_concentration.S2))
        @test all(isapprox.(update_inv_predict!(update_relative_signal!(cbatch)).data.estimated_concentration.var"G1(drug_b)", set_quantification(cbatch.data, cbatch).estimated_concentration.var"G1(drug_b)"))
    end
    @testset "Utils" begin
        @test getanalyte(cdata.area, AnalyteG1("G1(drug_b)")) == getanalyte(cdata.area, Symbol("G1(drug_b)"))
        @test getanalyte(rdata.area, AnalyteG1("G1(drug_b)")) == getanalyte(rdata.area, Symbol("G1(drug_b)"))
        @test getsample(cdata.area, "S2") == getsample(cdata.area, Symbol("S2"))
        @test getsample(rdata.area, Symbol("S2")) == getproperty(rdata.area.table, Symbol("S2"))
        @test all(isapprox.(dynamic_range(cbatch.calibration[1]), (1, 100)))
        @test endswith(formula_repr_utf8(cbatch.calibration[2]), "x^2")
        @test weight_repr_utf8(cbatch.calibration[1]) == "none"
        @test weight_value("1/x^3") == -3
    end
    @testset "IO" begin
        global initial_mc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_mc_c.batch"), Table)
        global initial_mc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_mc_r.batch"), Table)
        global initial_sc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_sc_c.batch"), Table)
        global initial_sc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "initial_sc_r.batch"), Table)
        update_quantification!(initial_mc_c)
        update_quantification!(initial_mc_r)
        update_quantification!(initial_sc_c)
        update_quantification!(initial_sc_r)
        global save_mc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_mc_c.batch"), Table)
        global save_mc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_mc_r.batch"), Table)
        global save_sc_c = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_sc_c.batch"), Table)
        global save_sc_r = ChemistryQuantitativeAnalysis.read(joinpath(datapath, "save_sc_r.batch"), Table)
        @test all(isapprox_nan.(collect(initial_mc_c.data.estimated_concentration.table[1]), collect(save_mc_c.data.estimated_concentration.table[1])))
        @test all(isapprox_nan.(collect(initial_mc_r.data.estimated_concentration.table[1]), collect(save_mc_r.data.estimated_concentration.table[1])))
        @test all(isapprox_nan.(collect(initial_sc_c.data.estimated_concentration.table[1]), collect(save_sc_c.data.estimated_concentration.table[1])))
        @test all(isapprox_nan.(collect(initial_sc_r.data.estimated_concentration.table[1]), collect(save_sc_r.data.estimated_concentration.table[1])))
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_mc_c.batch"), initial_mc_c)
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_mc_r.batch"), initial_mc_r)
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_sc_c.batch"), initial_sc_c)
        @test @test_noerror ChemistryQuantitativeAnalysis.write(joinpath(datapath, "save_sc_r.batch"), initial_sc_r)
        @test @test_noerror test_show(save_mc_c)
        @test @test_noerror test_show(save_mc_r)
        @test @test_noerror test_show(save_sc_c)
        @test @test_noerror test_show(save_sc_r)
        @test @test_noerror test_show(save_mc_c.calibration[1])
        @test @test_noerror test_show(save_mc_r.calibration[1])
        @test @test_noerror test_show(save_sc_c.calibration[1])
        @test @test_noerror test_show(save_sc_r.calibration[1])
        @test @test_noerror test_show(save_mc_c.method)
        @test @test_noerror test_show(save_mc_r.method)
        @test @test_noerror test_show(save_sc_c.method)
        @test @test_noerror test_show(save_sc_r.method)
        @test @test_noerror test_show(save_mc_c.data)
        @test @test_noerror test_show(save_mc_r.data)
        @test @test_noerror test_show(save_sc_c.data)
        @test @test_noerror test_show(save_sc_r.data)
    end
end
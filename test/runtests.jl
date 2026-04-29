using ChemistryQuantitativeAnalysis, TypedTables, DataFrames, Dictionaries, CSV, StatsAPI
using Test
import Base: show, convert
const CQA = ChemistryQuantitativeAnalysis

# sdt, adt
# Analyte 1: Linear, 1000 * x, AX = A1 / 1000 
# Analyte 2: ISD
# Analyte 3: Linear, 2000 * x
# Analyte 4: ISD
# Analyte 5: Linear, 1500 * x
# Analyte 6: Quadratic, 10 * AX ^ 2 
# Analyte 7: Quadratic, A6 - 0.1
# Analyte 8: Logarithmic, 200log(AX) + 460
# Analyte 9: Exponential, EXP(2 * AX + 1)
# Analyte 9: Power, 200 AX ^ 1.5

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

macro test_noerror(err, x)
    return quote
        try 
            $x
            true
        catch e
            isa(e, $err) ? true : false
        end
    end
end

include("objects.jl")

@testset "ChemistryQuantitativeAnalysis.jl" begin
    @testset "Constructors" begin
        @test all(isapprox.(sbatch.method.conctable.Analyte1, collect(abatch.method.conctable[1])[2:end]))
        @test getanalyte(cdata.area, 1) == getanalyte(rdata.area, 1) == getanalyte(SampleDataTable(rdata.area, :Sample, Table), 1)
        @test getsample(cdata.area, 2) == getsample(rdata.area, 2) == getsample(AnalyteDataTable(cdata.area, :Analyte, Table), 2)
        @test propertynames(cbatch) == (:method, :calibrator, :data, :analyte, :isd, :nonisd, :std, :point, :level)
        @test propertynames(cbatch.method) == (:analytetable, :signal,  :rel_sig, :est_conc, :nom_conc, :acc, :pointlevel, :conctable, :signaltable, :analyte, :isd, :nonisd, :std, :point, :level)
        @test propertynames(cdata) == (:area, )
        @test propertynames(cdata.area) == (:Sample, Symbol.(analytes)...)
        @test propertynames(rdata.area) == (:Analyte, Symbol.(["S1", "S2", "S3"])...)
        @test all(isapprox.(sbatch2.method.conctable.Analyte7 .* 2, sbatch2.method.conctable.Analyte2 .* 7))
    end
    @testset "Interface.jl" begin
        cdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 4
        @test cdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 4
        cdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 6.0
        @test cdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 6.0
        rdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 4
        @test rdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 4
        rdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 6.0
        @test rdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 6.0
        @test rdata.area[1] == (Analyte = analytes[1], S1 = 6.0, S2 = 24.0, S3 = 54.0)
        get!(cdata, :area, cdata.area)
        get(cdata, :area, nothing)
        @test haskey(cdata, :area)
        cdata[:area] = copy(cdata[:area])
        rdata[:area] = copy(rdata[:area])
        @test all([v for (k, v) in zip(keys(rdata), values(rdata))] .== [v for v in rdata])
        @test cbatch.isd == cbatch.method.isd
        @test cbatch.nonisd == cbatch.method.nonisd
        @test cbatch.point == cbatch.method.point
        @test cbatch.level == cbatch.method.level
        @test collect(cdata.area)[1].var"G1(drug_a)" == collect(rdata.area)[1].var"S1"
        @test collect(columns(cdata.area))[2][1] == collect(rows(rdata.area))[1][2]
        @test collect(columns(rdata.area))[2][1] == collect(rows(cdata.area))[1][2]
        # getcolumn, columnnames
        @test collect(eachanalyte(cdata.area)) == collect(eachanalyte(rdata.area))
        @test collect(eachsample(cdata.area)) == collect(eachsample(rdata.area))
        @test all(isapprox.(insert!(cbatch.data, :nominal_concentration, deepcopy(cbatch.data.area)).nominal_concentration.var"G1(drug_a)", set!(cbatch.data, :nominal_concentration, deepcopy(cbatch.data.area)).nominal_concentration.var"G1(drug_a)"))
        @test !in(:nominal_concentration, propertynames(delete!(cbatch.data, :nominal_concentration)))
        @test !in(:nominal_concentration, propertynames(unset!(cbatch.data, :nominal_concentration)))
    end
    @testset "Calibration" begin
        calibrate!(cbatch)
        calibrate!(rbatch)
        calibrate(cbatch, 1)
        @test all(CQA.validate_calibrator(cbatch.calibrator[1]) .== CQA.validate_calibrator(rbatch.calibrator[1]))
        @test all(CQA.quantify_calibrator(rbatch.calibrator[2]) .== CQA.quantify_calibrator(cbatch.calibrator[2]))
        model_calibrator!(cbatch; weight = XWeight())
        @test cbatch.calibrator[1].model.weight == XWeight()
        model_calibrator!(cbatch, Table(; analyte = [AnalyteG1("G1(drug_a)"), AnalyteG1("G1(drug_b)")], weight = [ConstWeight(), XWeight()]), model = LinearCalibrator)
        model_calibrator!(rbatch, AnalyteG1("G1(drug_b)"); model = QuadraticCalibrator)
        @test rbatch.calibrator[2].model isa QuadraticCalibrator
        edit_method_calibrate!(sbatch, Table(; 
            std = [1, -1, 3, -1, 5, 6, 7, 8, 9, 10, 11], 
            isd = [2, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0],
            model = [LinearCalibrator, Nothing, ProportionalCalibrator, Nothing, LinearCalibrator, QuadraticCalibrator, QuadraticOriginCalibrator, LogarithmicCalibrator, ExponentialCalibrator, PowerCalibrator, Nothing]))
        edit_method_calibrate!(abatch, 
            [], 
            [2, 4], 
            [1 => 2, 3 => 4], 
            Table(; analyte = abatch.std, model = [LinearCalibrator, ProportionalCalibrator, LinearCalibrator, QuadraticCalibrator, QuadraticOriginCalibrator, LogarithmicCalibrator, ExponentialCalibrator, PowerCalibrator, Nothing]))
        sbatch.data.area.Analyte2 .*= 10
        sbatch.method.signaltable.Analyte2 .*= 10
        for v in columns(abatch.data.area)
            v[2] /= 10
        end
        for v in columns(abatch.method.signaltable)
            v[2] /= 10
        end
        edit_method!(sbatch; model = LinearCalibrator, signal_threshold = 1e-3)
        edit_method_calibrate!(sbatch, Table(; analyte = sbatch.std, model = [LinearCalibrator, ProportionalCalibrator, LinearCalibrator, QuadraticCalibrator, QuadraticOriginCalibrator, LogarithmicCalibrator, ExponentialCalibrator, PowerCalibrator, Nothing]))
        calibrate!(abatch, "Analyte1")
        @test rbatch.calibrator[1].analyte == cbatch.calibrator[1].analyte
        @test all(isapprox.(cbatch.calibrator[1].table.accuracy[4:6], [1, 1.1, 0.9]))
        @test all(isapprox.(rbatch.calibrator[2].table.accuracy[4:6], [1, sqrt(1.1), sqrt(0.9)]))
        assign_isd_calibrate!(sbatch, AnalyteN(1), AnalyteN(4))
        @test sbatch.calibrator[1].isd == AnalyteN(4)
        @test sbatch.method.analytetable.isd[1] == 4
        assign_isd_calibrate!(sbatch, AnalyteN(1), AnalyteN(2))
        @test sbatch.calibrator[1].isd == AnalyteN(2)
        @test sbatch.method.analytetable.isd[1] == 2
        assign_std_calibrate!(sbatch, AnalyteN(1), AnalyteN(3))
        @test sbatch.method.analytetable.std[1] == 3
        @test sbatch.method.analytetable.isd[1] == 4
        assign_std_calibrate!(sbatch, AnalyteN(1), AnalyteN(1), AnalyteN(2))
        @test last(sbatch.calibrator).analyte == AnalyteN(1)
        @test sbatch.method.analytetable.std[1] == 1
        @test sbatch.method.analytetable.isd[1] == 2
        assign_std_calibrate!(sbatch, AnalyteN(3), AnalyteN(3), AnalyteN(4))
        @test sbatch.method.analytetable.std[3] == 3
        @test sbatch.method.analytetable.isd[3] == 4
        pushfirst!(sbatch.calibrator, pop!(sbatch.calibrator))
        @test last(sbatch.calibrator).machine isa CQA.EmptyMachine
        # InternalCalibrator
        @test isapprox(CQA.quantify_calibrator(ical)[1], getanalyte(method.conctable, 3)[1])
        @test CQA.analyze_calibrator!(Batch(method, [ical])).calibrator[1] == ical
        @test model_calibrator!(ical, Nothing) == ical
    end
    @testset "Quantification" begin
        # empty batch
        @test @test_noerror @test_logs (:info, "There is no data!") quantify!(ebatch)
        @test @test_noerror @test_logs (:info, "There is no data!") quantify_relative_signal!(ebatch)
        @test @test_noerror @test_logs (:info, "There is no data!") quantify_inv_predict!(ebatch)
        @test @test_noerror @test_logs (:info, "There is no data!") validate!(ebatch)
        @test @test_noerror @test_logs (:info, "There is no data!") analyze!(ebatch)
        # Batch: quantify, quantify!, set_quantify!, set_quantify, set_inv_predict, relative_signal
        @test all(isapprox.(getanalyte(set_relative_signal(rbatch.data, rbatch).relative_signal, AnalyteTest("G1(drug_a)")), cbatch.data.area.var"G1(drug_a)" ./ cbatch.data.area.var"G2(drug_a)"))
        set_relative_signal!(cbatch.data, cbatch)
        @test all(isapprox.(getanalyte(relative_signal(rbatch, rbatch.data), AnalyteTest("G1(drug_a)")), getanalyte(cbatch.data.relative_signal, AnalyteTest("G1(drug_a)"))))
        set_inv_predict!(cbatch.data, cbatch)
        set_quantify!(rbatch.data, rbatch)
        @test all(isapprox.(set_inv_predict(cbatch.data, cbatch).estimated_concentration.var"G1(drug_a)", getanalyte(rbatch.data.estimated_concentration, AnalyteTest("G1(drug_a)"))))
        @test all(isapprox.(getanalyte(set_quantify(rbatch.data, rbatch).estimated_concentration, AnalyteTest("G1(drug_a)")), quantify(cbatch, cbatch.data).var"G1(drug_a)"))
        @test all(isapprox.(analyze!(cbatch).data.estimated_concentration.var"G1(drug_a)", quantify(cbatch.calibrator[1], cbatch.data.area, cbatch.calibrator[1].analyte, cbatch.calibrator[1].isd)))
        
        # Cal: inv_predict, quantify
        # @test quantify(ical) == inv_predict(ical)
        @test quantify(ical, cbatch.data.area, AnalyteTest("G1(drug_a)")) == inv_predict(ical, cbatch.data.relative_signal, AnalyteTest("G1(drug_a)"))
        # @test quantify(cbatch.calibrator[2]) == inv_predict(cbatch.calibrator[2])
        @test quantify(cbatch.calibrator[2], cbatch.data.area) == inv_predict(cbatch.calibrator[2], cbatch.data.relative_signal)

        cbatch.data[:nominal_concentration] = deepcopy(cbatch.data.estimated_concentration)
        @test all(isapprox.(validate!(cbatch).data.accuracy.var"G1(drug_b)", set_accuracy(cbatch.data, cbatch).accuracy.var"G1(drug_b)"))
        @test all(isapprox.(accuracy(cbatch, cbatch.data).var"G1(drug_b)", accuracy(cbatch.method, cbatch.data).var"G1(drug_b)"))
    end
    @testset "Utils" begin
        @test isstd(AnalyteN(1), sbatch.method)
        @test isisd(AnalyteN(2), sbatch.method)
        @test isdof(AnalyteN(1), sbatch.method) == AnalyteN(2)
        @test stdof(AnalyteN(5), sbatch.method) == AnalyteN(5)
        @test getanalyte(cdata.area, AnalyteG1("G1(drug_b)")) == getanalyte(cdata.area, Symbol("G1(drug_b)"))
        @test getanalyte(rdata.area, AnalyteG1("G1(drug_b)")) == getanalyte(rdata.area, Symbol("G1(drug_b)"))
        @test getsample(cdata.area, "S2") == getsample(cdata.area, Symbol("S2"))
        @test getsample(rdata.area, Symbol("S2")) == getproperty(CQA.table(rdata.area), Symbol("S2"))
        @test getcalibrator(cbatch, AnalyteTest("G1(drug_a)")).analyte == cbatch.std[1]
        @test getcalibrator(cbatch, Symbol("G1(drug_a)")).isd == getcalibrator(rbatch, findcalibrator(rbatch, AnalyteTest("G1(drug_a)"))).isd
        @test samplename(rdata) == samplename(rdata.area)
        @test analytename(cdata) == analytename(cdata.area)
        @test all(isapprox.(dynamic_range(cbatch.calibrator[1]), (1, 100)))
        @test all(isapprox.(signal_range(cbatch.calibrator[2]), (signal_lloq(cbatch.calibrator[2]), signal_uloq(cbatch.calibrator[2]))))
        @test isapprox(first(signal_range(ical)), signal_lloq(ical))
        @test isinf(signal_uloq(ical))
        @test signal_loq(ical) == 0
        @test isapprox(signal_lob(sbatch.calibrator[2]) + signal_lod(sbatch.calibrator[2]),  0.4945 * signal_loq(sbatch.calibrator[2]))
        for cal in sbatch.calibrator
            global c = cal
            @test @test_noerror signal_lob(c)
            @test @test_noerror signal_lod(c)
            @test @test_noerror signal_loq(c)
        end 
        @test CQA.cqaconvert(String, "1") == "1"
        @test CQA.cqaconvert(Real, 1) == Real(1)
        @test CQA.cqamap(Int, [1, 2, 3]) == [1, 2, 3]
        @test CQA.cqamap(Real, [1, 2, 3]) == Real[1, 2, 3]
        @test CQA.cqamap(Int, Real[1, 2, 3]) == [1, 2, 3]
        @test CQA.cqamap(Int, ["1", "2", "3"]) == [1, 2, 3]
        @test CQA.cqamap(string, [1, 2, 3]) == string.([1, 2, 3])
        @test CQA.table(CQA.table_convert(DataFrame, cdata.area)) isa DataFrame
        @test CQA.table(CQA.table_convert(DataFrame, sdt)) isa DataFrame
        @test CQA.table(CQA.table_convert(DataFrame, adt)) isa DataFrame
        @test last(first(findoutofrange(cbatch, AnalyteG1("G1(drug_a)"); limit = x -> (8, 50)))) == [1]
        @test last(first(findunderrange(cbatch; limit = x -> lloq(x) * 8))) == [1]
        @test last(first(findoverrange(cbatch; limit = x -> uloq(x) / 3))) == [3]
        markoutofrange!(cbatch, AnalyteG1("G1(drug_a)"); limit = x -> (8, 50))
        markunderrange!(rbatch; limit = x -> lloq(x) * 8)
        markoverrange!(rbatch; limit = x -> uloq(x) / 3, value = (x, y) -> Inf)
        @test isnan(getanalyte(cbatch.data.estimated_concentration, AnalyteG1("G1(drug_a)"))[1])
        @test isnan(getanalyte(rbatch.data.estimated_concentration, AnalyteG1("G1(drug_a)"))[1])
        @test isinf(getanalyte(rbatch.data.estimated_concentration, AnalyteG1("G1(drug_b)"))[3])
    end
    @testset "IO" begin
        for w in [ConstWeight,
                    RootXWeight,
                    RootYWeight,
                    RootXYWeight,
                    XWeight,
                    YWeight,
                    XYWeight,
                    SqXWeight,
                    SqYWeight,
                    SqXYWeight,
                    RootLogXWeight,
                    RootLogYWeight,
                    LogXWeight,
                    LogYWeight,
                    SqLogXWeight,
                    SqLogYWeight,
                    RootExpXWeight,
                    RootExpYWeight,
                    ExpXWeight,
                    ExpYWeight,
                    SqExpXWeight,
                    SqExpYWeight]
            global wf = w()
            @test @test_noerror CQA.getweights(wf, 1.0, 1.0)
            @test @test_noerror CQA.human_name(wf)
            @test @test_noerror CQA.human_name_ascii(wf)
        end
        for cal in sbatch.calibrator
            global c = cal
            @test @test_noerror CQA.formula_repr_ascii(c)
        end 
        @test @test_noerror CQA.formula_repr_ascii(ical)
        @test endswith(CQA.formula_repr_ascii(rbatch.calibrator[2]), "x^2")
        @test @test_noerror CQA.read(joinpath(datapath, "initial_mc_c.batch", "data.at"), DataFrame)
        @test @test_noerror CQA.read(joinpath(datapath, "initial_mc_c.batch", "data.at", "0_area.sdt"), DataFrame)
        @test @test_noerror CQA.read(joinpath(datapath, "initial_mc_r.batch", "data.at", "0_area.adt"), DataFrame)
        @test @test_noerror CQA.read(joinpath(datapath, "initial_mc_c.batch", "method.am"), DataFrame)
        @test @test_noerror CQA.read(joinpath(datapath, "save_mc_c.batch", "calibrator", "1.ecal"), DataFrame)
        @test @test_noerror CQA.read(joinpath(datapath, "save_sc_c.batch", "calibrator", "1.ical"), DataFrame)
        quantify!(sbatch)
        quantify!(abatch)
        calibrate!(initial_mc_c)
        calibrate!(initial_mc_r)
        calibrate!(initial_sc_c)
        calibrate!(initial_sc_r)
        quantify!(initial_mc_c)
        quantify!(initial_mc_r)
        quantify!(initial_sc_c)
        quantify!(initial_sc_r)
        # value validation
        @test isapprox(CQA.table(initial_mc_c.data.estimated_concentration)[1, 2], collect(CQA.table(sbatch.data.estimated_concentration)[1])[2] * 10)
        @test all(isapprox_nan.(collect(CQA.table(initial_mc_c.data.estimated_concentration)[1, :]), collect(CQA.table(save_mc_c.data.estimated_concentration)[1, :])))
        @test isapprox(CQA.table(initial_mc_r.data.estimated_concentration)[1, 3], collect(CQA.table(abatch.data.estimated_concentration)[1])[1] * 10)
        @test all(isapprox_nan.(collect(CQA.table(initial_mc_r.data.estimated_concentration)[1, :]), collect(CQA.table(save_mc_r.data.estimated_concentration)[1, :])))
        @test all(isapprox_nan.(collect(CQA.table(initial_sc_c.data.estimated_concentration)[1, :]), collect(CQA.table(save_sc_c.data.estimated_concentration)[1, :])))
        @test all(isapprox_nan.(collect(CQA.table(initial_sc_r.data.estimated_concentration)[1, :]), collect(CQA.table(save_sc_r.data.estimated_concentration)[1, :])))
        model_calibrator!(save_mc_c; model = PowerCalibrator)
        model_calibrator!(save_mc_r; weight = SqXWeight())
        @test all(isapprox.(StatsAPI.predict(save_mc_c.calibrator[1], Table(; x = [1, 2, 3])), StatsAPI.predict(save_mc_c.calibrator[1], [1, 2, 3])))
        @test all(isapprox.(StatsAPI.predict(save_mc_c.calibrator[1].machine, Table(; x = [1, 2, 3])), StatsAPI.predict(save_mc_c.calibrator[1].machine, [1, 2, 3])))
        calibrate!(save_mc_c)
        calibrate!(save_mc_r)
        calibrate!(save_sc_c)
        calibrate!(save_sc_r)
        quantify!(save_mc_c)
        quantify!(save_mc_r)
        quantify!(save_sc_c)
        quantify!(save_sc_r)
        # write
        @test @test_noerror SystemError CQA.write(joinpath(datapath, "save_mc_c.batch"), initial_mc_c)
        @test @test_noerror SystemError CQA.write(joinpath(datapath, "save_mc_r.batch"), initial_mc_r)
        @test @test_noerror SystemError CQA.write(joinpath(datapath, "save_sc_c.batch"), initial_sc_c)
        @test @test_noerror SystemError CQA.write(joinpath(datapath, "save_sc_r.batch"), initial_sc_r)
        # show
        @test @test_noerror test_show(save_mc_c)
        for _ in 1:10
            push!(save_sc_r.calibrator, save_sc_r.calibrator[1]) # Many calibrator
        end
        delete!(save_sc_r.data, :estimated_concentration) # <3 table
        @test @test_noerror test_show(save_sc_r)
        @test @test_noerror test_show(save_mc_c.calibrator[1])
        @test @test_noerror test_show(save_sc_r.calibrator[1])
        @test @test_noerror test_show(save_mc_c.method)
        @test @test_noerror test_show(save_sc_c.method)
        @test @test_noerror test_show(save_mc_c.data)
        @test @test_noerror test_show(save_mc_r.data)
        @test @test_noerror test_show(sbatch.calibrator[1].model)
        @test @test_noerror test_show(sbatch.calibrator[2].model)
        @test @test_noerror test_show(sbatch.calibrator[3].model)
        @test @test_noerror test_show(sbatch.calibrator[4].model)
        @test @test_noerror test_show(sbatch.calibrator[5].model)
        @test @test_noerror test_show(sbatch.calibrator[6].model)
        @test @test_noerror test_show(sbatch.calibrator[7].model)
        @test @test_noerror test_show(sbatch.calibrator[8].model)
        # mkbatch validation
        @test @test_noerror SystemError mkbatch(joinpath(datapath, "new1.batch"); data_table = AnalyteDataTable, signal_table = AnalyteDataTable, conc_table = AnalyteDataTable, data_config = Dict(:Sample => ["S0", "S1", "S2"]), signal_config = Dict(:Sample => ["C0"]), conc_config = Dict(:Sample => [1]))
        @test @test_noerror SystemError mkbatch(joinpath(datapath, "new2.batch"); conc_table = AnalyteDataTable, data_config = Dict(:Analyte => ["A1", "A2"]), signal_config = Dict(:Analyte => ["A1", "A2", "A3"]), conc_config = Dict(:Sample => [1, 2, 3, 4, 5, 6]))
        @test @test_noerror SystemError mkbatch(joinpath(datapath, "new3.batch"); data_table = AnalyteDataTable, signal_table = AnalyteDataTable, data_config = Dict(:Sample => ["S0", "S1", "S2"]), conc_config = Dict(:Analyte => ["A1", "A2", "A3"]), signal_config = Dict(:Sample => ["C1", "C2", "C3", "C4", "C5", "C6"]))
        @test @test_noerror SystemError mkbatch(joinpath(datapath, "new4.batch"); method_config = Dict(:signal => "height", :levelname => "L"))
        n1 = CQA.read(joinpath(datapath, "new1.batch"), DataFrame)
        n2 = CQA.read(joinpath(datapath, "new2.batch"), DataFrame)
        n3 = CQA.read(joinpath(datapath, "new3.batch"), DataFrame)
        n4 = CQA.read(joinpath(datapath, "new4.batch"), DataFrame)
        calibrate!(n2)
        @test sampleobj(n1.data) == ["S0", "S1", "S2"]
        @test length(unique(n2.calibrator[1].table.level[n2.calibrator[1].table.include])) == 6
        @test analyteobj(n3.data) == ["A1", "A2", "A3"]
        @test :L in propertynames(n4.method.signaltable)
    end
end
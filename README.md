# ChemistryQuantitativeAnalysis

[![Build Status](https://github.com/yufongpeng/ChemistryQuantitativeAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yufongpeng/ChemistryQuantitativeAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yufongpeng/ChemistryQuantitativeAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yufongpeng/ChemistryQuantitativeAnalysis.jl)

`ChemistryQuantitativeAnalysis.jl` is a package for quantitative analysis of chemicals based on tabular data. 

Graphical user interface is available for adjusting calibration curve. See [`ChemistryQuantitativeAnalysisUI.jl`](https://github.com/yufongpeng/ChemistryQuantitativeAnalysisUI.jl).

Command line interfaces are implemented in [`juliaquant`](https://github.com/yufongpeng/juliaquant) for non-programmers.

## Tabular data wrapper
This package provides two wrappers for data, `SampleDataTable{A, S, N, T}` and `AnalyteDataTable{A, S, N, T}` which are subtypes of `AbstractDataTable{A, S, N, T}`. `SampleDataTable` indicates that part of columns represent analytes, and all rows reprsent samples. `AnalyteDataTable` indicates that part of columns represent samples, and all rows represent analytes. Both types behave like the underlying `table` object. 
|Fields|`SampleDataTable{A, S, N, T}`|`AnalyteDataTable{A, S, N, T}`|
|----------|---------------------|------------------|
|`analytecol`|-|`Symbol`, column name whose elements are analytes.|
|`samplecol`|`Symbol`, column name whose elements are samples.|-|
|`analyte`|`Vector{A}`, analytes in user-defined types.|same|
|`sample`|`Vector{S}`, samples in user-defined types.|same|
|`table`|Tabular data of type `T`|same|

`SampleDataTable{A, S, N, T}` can be constructed in the following ways:
1. `SampleDataTable(analytetype::TypeOrFn, samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))`
2. `SampleDataTable(table, analytetype::TypeOrFn, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))`
3. `SampleDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))`
4. `SampleDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))`

`AnalyteDataTable{A, S, N, T}` can be constructed in the following ways:
1. `AnalyteDataTable(analytecol::Symbol, sampletype::TypeOrFn, table; samplename = setdiff(propertynames(table), [analytecol]))`
2. `AnalyteDataTable(table, analytecol::Symbol, sampletype::TypeOrFn; samplename = setdiff(propertynames(table), [analytecol]))`
3. `AnalyteDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol]))`
4. `AnalyteDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol]))`

`analaytename` and `samplename` will be converted to a vector of string, and then converted to the desired type using `cqaconvert(analytetype, analaytename)` and `cqaconvert(sampletype, samplename)`.

## AnalysisMethod
`AnalysisMethod{A, M, C, D}` is used for storing method, containing all analytes, their internal standards and calibration curve setting, and data for fitting calibration curve.
|Property|Description|
|----------|---------|
|`analytetable::M <: Table`|table with at least 4 columns, `analytes` identical to property `analytes`, `isd`, matching each analyte to index of its internal standard, and `std` matching each analyte to the index of calibration standard. `-1` indicates the analyte itself is internal standard, and `0` indicates no internal standard. `model` stores the calibration model. For example, a row `(analytes = AnalyteX, isd = 2, std = 3, model = LinearCalibrator(ConstWeight()))` means that internal standard of `AnalyteX` is the second analyte, and it will be external calibrated by the third analyte with ordinary linear regression. If isd and std are the same, internal calibrator is used.|
|`signal::Symbol`|type and key name of experimental acquisition data.|
|`rel_sig::Symbol`|key name of  relative signal.|
|`est_conc::Symbol`|key name of  estimated concentration.|
|`true_conc::Symbol`| key name of  true concentration.|
|`acc::Symbol`|key name of accuracy.|
|`pointlevel::Vector{Int}`| matching each point to level. It can be empty if there is only one level in `conctable`.|
|`conctable::C <: AbstractDataTable{A, Int}`| containing concentration data for each level. Sample names must be symbol or string of integers for multiple levels.|
|`signaltable::D <: AbstractDataTable{A, S}`| containig signal for each point. It can be `nothing` if signal data is unecessary.|
|`analyte::AbstractVector{A}`|analytes in user-defined types.|
|`isd::AbstractVector{A}`|internal standards.|
|`nonisd::AbstractVector{A}`|analytes which are not internal standards.|
|`std::AbstractVector{A}`|analytes which are calibration standards.|
|`point::AbstractVector{S}`|calibration points, identical to `signaltable.samples`. If `signaltable` is `nothing`, this value is `nothing` too.|
|`level::AbstractVector{Int}`|calibration levels, identical to `conctable.samples`.|

Constructors of `AnalysisMethod`:
1. `AnalysisMethod(analytetable::Table, signal::Symbol, rel_sig::Symbol, est_conc::Symbol, true_conc::Symbol, acc::Symbol, pointlevel::Vector{Int}, conctable::AbstractDataTable{A, Int}, signaltable::Union{AbstractDataTable, Nothing})`
2. `AnalysisMethod(conctable::AbstractDataTable, signaltable::SampleDataTable, signal, levelname; kwargs...)`
3. `AnalysisMethod(conctable::AbstractDataTable, signaltable::Nothing; kwargs...)`
4. `AnalysisMethod(conctable::AbstractDataTable, signaltable::Union{AbstractDataTable, Nothing}, signal, pointlevel::AbstractVector{Int}; kwargs...)`

Keyword arguments
* `rel_sig`: key name of relative signa. It defaults to `:relative_signal`.
* `est_conc`: key name of estimated concentration. It defaults to `:estimated_concentration`.
* `true_conc`: key name of true concentration. It defaults to `:true_concentration`. 
* `acc`: key name of accuracy. It defaults to `:accuracy`.
* Other keyword arguments will be columns in `analytetable`; when `analyte`, `isd` and `std` are not provided, it will use analyte in `conctable`. 

`levelname` is the column name for `pointlevel` if `signaltable` is a `SampleDataTable`.

## AnalysisTable
`AnalysisTable{A, S, T}` is basically a `Dictionary{Symbol, <: AbstractDataTable{T}}` which data can be extracted using proeperty syntax. For example, `at[:area] === at.area`.
|Field|Description|
|----------|---------|
|`analyte::Vector{A}`|analytes in user-defined types.|
|`sample::Vector{S}`|samples in user-defined types.|
|`tables::Dictionary{Symbol, <: AbstractDataTable{T}}`|a dictionary mapping data type to datatable.|

The key names are determined by an `AnalysisMethod`. 

`AnalysisTable{A, S, T}` can be constructed in the following ways:
1. `AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S}})`
2. `analysistable(iter)`

`iter` is an iterable iter of key-value `Pair`s (or other iterables of two elements, such as a two-tuples). Keys should be `Symbol`s, and values should be `AbstractDataTable`s.

## Calibrator
This package provides two calibrator types, `ExternalCalibrator{A, N, T}` and `InternalCalibrator{A, N}` which are subtypes of `AbstractCalibrator{A, N}`.

### ExternalCalibrator
This type is for external calibration with or without internal calibration. It can be constructed from a `AnalysisMethod{A, S}` containing calibration data, and analyte `A` using function `calibrate`.
|Field|Description|
|----------|-----------|
|`analyte`|analyte being quantified|
|`isd`|internal standard for which `nothing` indicates no internal standard|
|`table::TypedTable.Table`|the cleaned up calibration data, containing 7 columns|
|`model::AbstractCalibrationModel`|calibration model|
|`machine`|calibration machine. Either a `LsqFitMachine` or fitted `GLM` object|


The columns in `table`:
|Column|Description|
|------|-----------|
|`id`|Point name|
|`level`|The index of concentration level. The larger, the higher concentraion it represents.|
|`y`|Signal or relative signal|
|`x`|True concentraion|
|`x̂`|Predicted concentration|
|`accuracy`|Accuracy, i.e. `x̂/x`.|
|`include`|Whether this point is included or not|

To predict concentration, call `quantify_calibrator` and `quantify_calibrator!`. To calculate accuracy, call `validate_calibrator` and `validate_calibrator!`. `analyze_calibrator!` applys `quantify_calibrator!` and `validate_calibrator!` subsequently. 

### InternalCalibrator
This type contains data for single point internal calibration. 
|Field|Description|
|----------|-----------|
|`analyte::Tuple{A}`|the analyte with known concentration (internal standard).|
|`conc::Float64`|concentration of analyte.|

## Batch
`Batch{A, M, C, D}` represents a batch for quantitative analysis.
|Property|Description|
|----------|-----------|
|`method::M <: AnalysisMethod{A}`|method|
|`calibrator::C <: Union{AbstractVector{ExternalCalibrator{A}}, AbstractVector{InternallCalibrator{A}}}`|calibrators|
|`data::D <: Union{AnalysisTable{A}, Nothing`|Data for analysis|
|`analyte::AbstractVector{A}|analytes in user-defined types, identical to `method.analyte`.|
|`isd::AbstractVector{A}`|analytes which are internal standards, identical to `method.isd`.|
|`nonisd::AbstractVector{A}`|analytes which are not internal standards, identical to `method.nonstd`.|
|`std::AbstractVector{A}`|analytes which are calibration standards, identical to `method.std`.|
|`point::Union{Nothing, AbstractVector{S}}`|calibration points, identical to `method.point`.|
|`level::AbstractVector{Int}`|calibration levels, identical to `method.level`.|

Constructors for `Batch{A, M, C, D}`:
1. `Batch(method::M, calibrator::C, data::D = nothing)`
2. `Batch(method::AnalysisMethod, data = nothing; type = true, zero = false, weight = 0)`
3. `Batch(batch::Batch, at::AnalysisTable)`
4. `Batch(dt::AbstractDataTable; signal = :area, rel_sig = :relative_signal, est_conc = :estimated_concentration, true_conc = :true_concentration, acc = :accuracy, calid = r"Cal_(\d)_(\d*-*\d*)", order = "LR", level = nothing, ratio = nothing, dilution_factor = nothing, conc_factor = 1, parse_decimal = x -> replace(x, "-" => "."))`

The last method allows user to use encoded sample names or additional tabular data from `dt` to generate `AnalysisMethod`. 

Note that the constructor does not automatically fit any calibration curves if no calibrators are given. User can edit `analytetable` with `edit_method!`, `assign_isd!`, `assign_std!`, and then apply `calibrate!` to start calibration or in combination, `..._calibrate!`.

To calculate relative signal, concentration or accuracy and save the result, call `quantify_relative_signal!`, `quantify_inv_predict!` (in combination, `quantify!!`) and `validate!` (in combination with `quantify!, `analyze!`), respectively.

## Reading and writting data to disk
To use data on disk, user should create a directory in the following structure:
```
batch_name.batch
├──config.txt
├──method.am
│  ├──true_concentration.sdt
│  │  ├──config.txt
│  │  └──table.txt
│  ├──area.sdt
│  │  ├──config.txt
│  │  └──table.txt
│  ├──analytetable.txt
│  └──config.txt
├──calibrator
│  ├──1.ecal
│  └──2.ical
└──data.at
   ├──0_quantity1.sdt
   ├──1_quantity2.sdt 
   └──2_quantity3.sdt
```
There is a function `mkbatch` which creates a valid batch directory programatically.

Config files have the following general forms
```
[header1]
value

[header2]
value1
value2
value3
⋮
```
The first config file contains a header `delim` which determines the default delimiter for `table.txt` in this directory and subdirectories.

`data.at` and `calibrator` is not necessary for initializing a batch. The former can be added to the batch directly in julia, and the latter will be generated after calibration.

All `*.sdt` files can be replaced with `*.adt` files.

### *.sdt
All `*.sdt` files will be read as `SampleDataTable`.

Config file needs the following headers.
```
[delim]
\t

[Sample]
sample_col_name

[Analyte]
analyte_col_name_1
analyte_col_name_2
⋮
``` 

### *.adt
All `*.adt` files will be read as `AnalyteDataTable`.

Config file needs the following headers.
```
[delim]
\t

[Analyte]
analyte_col_name

[Sample]
sample_col_name_1
sample_col_name_2
⋮
``` 

### *.am
It must contain two `.sdt` or `.adt` files. A file must contain nominal concentration for each analyte and level. The sample names must be integers.
Another file is signal data for each analyte and calibration point. The file name is determined by `config.txt`.

Config file for `method.am` needs the following headers.
```
[signal]
area

[rel_sig]
relative_signal

[est_conc]
estimated_concentration

[nom_conc]
nominal_concentration

[acc]
accuracy

[delim]
\t

[levelname]
level

[pointlevel]
level_for_1st_point
level_for_2nd_point
⋮
```
`signal`, `rel_sig`, `est_conc`, `nom_conc`, and `acc` specifys which `.sdt` or `.adt` files containing corresponding data in the method directory or ascociated `**.at` directory. 
For the above file, `method.am/area.sdt` or `method.am/area.adt` will become `method.signaltable`; `method.am/nominal_concentration.sdt` or `method.am/nominal_concentration.adt` will become `method.conctable`.

`pointlevel` maps each point to level which should be integers.

`levelname` specifys the column representing property `pointlevel` of `AnalysisMethod`. It only works for which `signaltable` is `SampleDataTable`; otherwise, it falls back to use `pointlevel`.

`analytetable.txt` needs to contain analyte names, index of their internal standards, index of calibration standards, and calibration model. The column names are fixed for these four columns.
```
analytes isd   std model
analyte1 isd1  std1 LinearCalibrator(ConstWeight())
analyte2 isd2  std2 LinearCalibrator(ConstWeight())
⋮
```
The delimiter should be "\t", and the order of columns does not matter.

### *.at
It can contain multiple `*.sdt` or `*.adt`. The file names must start from an integer with `_` following the name, e.g. `0_area.sdt`. The integer is for the order of reading into `AnalysisTable`, and `name` will be the key. The name of data is determined in `method.am/config.txt`.

### Reading and writing Batch
To read a batch into julia, call `ChemistryQuantitativeAnalysis.read`.
```julia-repl
julia> batch = ChemistryQuantitativeAnalysis.read("batch_name.batch", T; table_type, analytetype, sampletype, numbertype, modeltype, delim)
```
`T` is the sink function for tabular data; it should create an object following `Tables.jl` interface. `table_type` is `T` parameter in the type signature of `Batch` which determines the underlying table type, `analytetype` is a concrete type for `analyte`, `sampletype` is a concrete type for `sample`, `numbertype` is a concrete type for numeric data, `modeltype` specifys calibration model type, and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist. 

For any `x` of type `analytetype` or `sampletype`, `x` equals `cqaconvert(type, string(x)))`. Additionally, `tryparse` have to be extended for `CSV` parsing:
* `tryparse(::Type{analytetype}, x::AbstractString)` is neccessary for `AnalyteDataTable`.
* `tryparse(::Type{sampletype}, x::AbstractString)` is neccessary for `SampleDataTable`.

To write batch to disk, call `ChemistryQuantitativeAnalysis.write`. There is a keyword argument `delim` controling delimiter of tables.
```julia-repl
julia> ChemistryQuantitativeAnalysis.write("batch_name.batch", batch; delim = '\t')
```
There will be a folder `calibrator` containing multiple `*.ecal` or `*.ical` folders. The former is for `ExternalCalibrator` and the latter is for `InternallCalibrator`.

## Examples
```julia
using ChemistryQuantitativeAnalysis, TypedTables, DataFrames
const CQA = ChemistryQuantitativeAnalysis
import Base: show, convert, tryparse

# Custom Analyte type
struct AnalyteG1
    name::String
end
struct AnalyteG2
    name::String
end
struct AnalyteOther
    name::String
end
const AnalyteTest = Union{AnalyteG1, AnalyteG2, AnalyteOther} # Use Union rather than abstarct type because csv parser only support concrete type
show(io::IO, analyte::AnalyteTest) = print(io, analyte.name)

# Analyte parser
function (::Type{AnalyteTest})(name::String)
    g = match(r"^G(\d)\(.*\)$", name)
    isnothing(g) && return AnalyteOther(name)
    g = parse(Int, first(g))
    g == 1 ? AnalyteG1(name) : g == 2 ? AnalyteG2(name) : AnalyteOther(name)
end
tryparse(::Type{AnalyteTest}, s::String) = AnalyteTest(s) # For reading data from disk

# Generate data
analytes = typedmap(AnalyteTest, ["G1(drug_a)", "G2(drug_a)", "G1(drug_b)", "G2(drug_b)"])
conc = Float64[0, 1, 2, 5, 10, 20, 50, 100]
signal1 = vcat(0.001, Float64[1, 2, 5, 10, 20, 50, 100], 0.005, [1, 2, 5, 10, 20, 50, 100] .+ 0.1, 0.002, [1, 2, 5, 10, 20, 50, 100] .- 0.1)
signal2 = vcat(0.002, Float64[1, 2, 5, 10, 20, 50, 100] .^ 2, 0.001, [1, 2, 5, 10, 20, 50, 100] .^ 2 .+ 0.1, 0.005, [1, 2, 5, 10, 20, 50, 100] .^ 2 .- 0.1)

# Create method
conctable = SampleDataTable(
      DataFrame(
         "level" => collect(0:7), 
         "G1(drug_a)" => conc,
         "G1(drug_b)" => conc .* 10, 
         "G2(drug_a)" => repeat([50.0], 8),
         "G2(drug_b)" => repeat([50.0], 8)),
      AnalyteTest,
      :level
)
signaltable = SampleDataTable(
      DataFrame(
         "point" => reshape([string(a, "_", b) for (a, b) in Iterators.product(0:7, 1:3)], 24), 
         "level" => repeat(0:7, 3),
         "G1(drug_a)" => signal1,
         "G2(drug_a)" => repeat([5.0], 24),
         "G1(drug_b)" => signal2,
         "G2(drug_b)" => repeat([2.0], 24)), 
      AnalyteTest,
      :point; 
      analytename = Symbol.(analytes)
)
method = AnalysisMethod(conctable, signaltable, :area, :level; analyte = analytes, isd = [2, -1, 4, -1], std = [1, -1, 3, -1])
cdata = AnalysisTable([:area], [
      SampleDataTable(
         DataFrame(
            "Sample" => ["S1", "S2", "S3"], 
            "G1(drug_a)" => Float64[6, 24, 54],
            "G2(drug_a)" => Float64[5, 6, 6],
            "G1(drug_b)" => Float64[200, 800, 9800],
            "G2(drug_b)" => Float64[2, 2, 2]), 
         AnalyteTest,
         :Sample
         )
      ]
)
rdata = analysistable((:area => 
      AnalyteDataTable(
         DataFrame(
            "Analyte" => analytes, 
            "S1" => Float64[6, 6, 200, 2],
            "S2" => Float64[24, 6, 800, 2],
            "S3" => Float64[54, 6, 9800, 2]
            ), 
         :Analyte
         )
      , )
)

nominal_conc = SampleDataTable(
         DataFrame(
            "Sample" => ["S1", "S2", "S3"], 
            "G1(drug_a)" => Float64[5, 25, 50],
            "G2(drug_a)" => Float64[5, 5, 5],
            "G1(drug_b)" => Float64[150, 300, 980],
            "G2(drug_b)" => Float64[2, 2, 2]), 
         AnalyteTest,
         :Sample
         )

# Create batch
cbatch = Batch(method, cdata)
rbatch = Batch(deepcopy(method), rdata)

# Calibration
calibrate!(cbatch)
calibrate!(rbatch)
cbatch.calibrator # a vector of `ExternalCalibrator`
model_calibrator!(cbatch, 3; model = QuadraticCalibrator)
model_calibrator!(rbatch, AnalyteG1("G1(drug_b)"); weight = XWeight())

# Quantification
quantify_relative_signal!(cbatch) # A new data `cbatch.data.relative_signal` is created.
quantify_inv_predict!(cbatch) # Fit `cbatch.data.relative_signal` into calibration curve to create `cbatch.data.estimated_concentration`.
quantify!(cbatch) # equivalent to `update_inv_predict!(update_relative_signal!(cbatch))`
cbatch.data[:nominal_concentration] = nominal_conc
validate!(cbatch)
analyze!(cbatch)

# Utils
analyteobj(cdata.area) # analytes of type `AnalyteTest`
sampleobj(cdata) # samples of type `String`
analytename(cdata) # analytes of type `Symbol`
samplename(cdata.area) # samples of type `Symbol`
cdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 6 # Set value using `dt[analyte, sample]`
cdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 6
collect(eachanalyte(cdata.area))
collect(eachsample(cdata.area))
getanalyte(cdata.area, AnalyteG1("G1(drug_b)")) # get data of `AnalyteG1("G1(drug_b)")`
getanalyte(cdata.area, 1) # get data of first analyte
getsample(cdata.area, "S2")
getcalibrator(cbatch, AnalyteTest("G1(drug_a)"))
dynamic_range(cbatch.calibrator[1])
signal_range(rbatch.calibrator[2])
signal_lob(rbatch.calibrator[2])
signal_lod(rbatch.calibrator[2])
signal_loq(rbatch.calibrator[2])
signal_lloq(rbatch.calibrator[2])
signal_uloq(rbatch.calibrator[2])
lloq(rbatch.calibrator[2])
uloq(rbatch.calibrator[2])
```
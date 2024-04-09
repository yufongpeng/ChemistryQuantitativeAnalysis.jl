# ChemistryQuantitativeAnalysis

[![Build Status](https://github.com/yufongpeng/ChemistryQuantitativeAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yufongpeng/ChemistryQuantitativeAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yufongpeng/ChemistryQuantitativeAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yufongpeng/ChemistryQuantitativeAnalysis.jl)

`ChemistryQuantitativeAnalysis.jl` is a package for quantitative analysis of chemicals based on tabular data. 

For command line interfaces, see [`juliaquant`](https://github.com/yufongpeng/juliaquant).

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
|`analytetable`|`M <: Table` with at least 3 columns, `analytes` identical to property `analytes`, `isd`, matching each analyte to index of its internal standard, and `calibration` matching each analyte to index of other analyte for fitting its calibration curve. `-1` indicates the analyte itself is internal standard, and `0` indicates no internal standard. For example, a row `(analytes = AnalyteX, isd = 2, calibration = 3)` means that internal standard of `AnalyteX` is the second analyte, and it will be quantified using calibration curve of the third analyte.|
|`signal`|`Symbol`, propertyname for extracting signal data from an `AnalysisTable`|
|`pointlevel`|`Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.|
|`conctable`|`C <: AbstractDataTable{A, Int}` containing concentration data for each level. Sample names must be symbol or string of integers for multiple levels. One level indicates using `SingleCalibration`.|
|`signaltable`|`D <: AbstractDataTable{A, S}` containig signal for each point. It can be `nothing` if signal data is unecessary.|
|`analyte`|`AbstractVector{A}`, analytes in user-defined types.|
|`isd`|`AbstractVector{<: A}` that each analytes are internal standards.|
|`nonisd`|`AbstractVector{<: A}` that each analytes are not internal standards.|same|
|`point`|`AbstractVector{S}`, calibration points, identical to `signaltable.samples`. If `signaltable` is `nothing`, this value is `nothing` too.|
|`level`|`AbstractVector{Int}`, calibration levels, identical to `conctable.samples`.|

Constructors of `AnalysisMethod`:
1. `AnalysisMethod(analytetable::Table, signal::Symbol, pointlevel::Vector{Int}, conctable::AbstractDataTable{A, Int}, signaltable::Union{AbstractDataTable, Nothing})`
2. `AnalysisMethod(conctable::AbstractDataTable{A, Int}, signaltable::SampleDataTable, signal::Symbol, levelname::Symbol; kwargs...)`
3. `AnalysisMethod(conctable::AbstractDataTable{A, Int}, signaltable::Nothing, signal::Symbol; kwargs...)`
4. `AnalysisMethod(conctable::AbstractDataTable, signaltable::Union{AbstractDataTable, Nothing}, signal::Symbol, pointlevel::AbstractVector{Int}; kwargs...)`

`kwargs` will be columns in `analytetable`. `levelname` is the column name for `pointlevel` if `signaltable` is a `SampleDataTable`.

## AnalysisTable
`AnalysisTable{A, S, T}` is basically a `Dictionary{Symbol, <: AbstractDataTable{T}}` which data can be extracted using proeperty syntax. For example, `at[:area] === at.area`.
|Field|Description|
|----------|---------|
|`analyte`|`Vector{A}`, analytes in user-defined types.|
|`sample`|`Vector{S}`, samples in user-defined types.|
|`tables`|`Dictionary{Symbol, <: AbstractDataTable{T}}`, a dictionary mapping data type to datatable.|

The key for signal data is determined by `method.signal`. Default names for relative signal, true concentration, estimated concentration and accuracy are `relative_signal`, `true_concentration`, `estimated_concentration` and `accuracy`.

`AnalysisTable{A, S, T}` can be constructed in the following ways:
1. `AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S}})`
2. `analysistable(iter)`

`iter` is an iterable iter of key-value `Pair`s (or other iterables of two elements, such as a two-tuples). Keys should be `Symbol`s, and values should be `AbstractDataTable`s.

## Calibration
This package provides two calibration types, `MultipleCalibration{A, N, T}` and `SingleCalibration{A, N}` which are subtypes of `AbstractCalibration{A, N}`.

### MultipleCalibration
This type fits and stores calibration curve. It can be created from a `AnalysisMethod{A, S}` containing calibration data, an analyte `A` using function `calibration`.
|Field|Description|
|----------|-----------|
|`analyte`|`Tuple{A, Any}`. First element is the analyte being quantified, and the second element is its internal standard for which `nothing` indicates no internal standard.|
|`type`|`Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`).|
|`zero`|`Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`).|
|`weight`|`Float64` represents the exponential applying to each element of `x` as a weighting vector.|
|`formula`|`FormulaTerm`, the formula for fitting calibration curve.|
|`table`|`TypedTable.Table`, the clean up calibration data, containing 7 columns.|
|`model`|`GLM` object|

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

To predict concentration, call `inv_predict`. To calculate accuracy, call `accuracy`. `type`, `zero`, and `weigtht` can be modified directly. To change internal standard, modify `analyte`. After any modification, call `update_calibration!` with method to update the `model`.

### SingleCalibration
This type contains data for single pont calibration. 
|Field|Description|
|----------|-----------|
|`analyte`|`Tuple{A}` is the analyte with known concentration (internal standard).|
|`conc`|`Float64`, concentration of analyte.|

## Batch
`Batch{A, M, C, D}` represents a batch for quantitative analysis.
|Property|Description|
|----------|-----------|
|`method`|`M <: AnalysisMethod{A}`, method.|
|`calibration`|`C <: Union{AbstractVector{MultipleCalibration{<: A}}, AbstractVector{SingleCalibration{<: A}}}`|
|`data`|Data for analysis, `D <: Union{AnalysisTable{<: A}, Nothing}`.|
|`analyte`|`AbstractVector{A}`, analytes in user-defined types, identical to `method.analytetable.analyte`.|
|`isd`|`AbstractVector{<: A}`, analytes which are internal standards, identical to `method.analytetable.analyte`.|
|`nonisd`|`AbstractVector{<: A}` that each analytes are not internal standards.|same|
|`point`|`AbstractVector{S}` or `Nothing`, calibration points, identical to `method.point`.|
|`level`|`AbstractVector{Int}`, calibration levels, identical to `method.level`.|

Constructors for `Batch{A, M, C, D}`:
1. `Batch(method::M, calibration::C, data::D = nothing)`
2. `Batch(method::AnalysisMethod, data = nothing; type = true, zero = false, weight = 0)`
3. `Batch(batch::Batch, at::AnalysisTable)`
4. `Batch(dt::AbstractDataTable; signal::Symbol = :area, calid = r"Cal_(\d)_(\d*-*\d*)", order = "LR", f2c = 1, parse_decimal = x -> replace(x, "-" => "."))`

The last method allows user to use encoded sample names to generate `AnalysisMethod`. Note that the returned batch does not have any calibration curves, which allows user to modify `analytetable`, and then apply `init_calibration!` to start calibrate.

To calculate relative signal, concentration or accuracy and save the result, call `update_relative_signal!`, `update_inv_predict!` (in combination, `update_quantification!`) and `update_accuracy!`, respectively.

## User Interface
Use the function `ui_init` which activates an environment for ui and run `interactive_calibrate!` on the batch.

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
├──calibration
│  ├──1.mcal
│  └──2.mcal
└──data.at
   ├──0_quantity1.sdt
   ├──1_quantity2.sdt 
   └──2_quantity3.sdt
```
There is a function `mkbatch` which create a valid batch directory programatically.

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
The header `delim` determines the default delimiter for `table.txt` in this directory and subdirectories.

`data.at` and `calibration` is not necessary for initializing a batch. The former can be added to the batch directly in julia, and the latter will be generated after calibration.

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
It must contain two `.sdt` or `.adt` files. `true_concentration.sdt` or `true_concentration.adt` contains true concentration for each analyte and level. The sample names must be integers.
Another file is signal data for each analyte and calibration point. The file name is determined by `config.txt`.

Config file for `method.am` needs the following  headers.
```
[signal]
area

[delim]
\t

[levelname]
level

[pointlevel]
level_for_1st_point
level_for_2nd_point
⋮
```
`signal` specifys which `.sdt` or `.adt` file serving as signal data. For the above file, `method.am/area.sdt` or `method.am/area.adt` will become `method.signaltable`.

`pointlevel` maps each point to level which should be integers.

`level` specifys the column representing property `pointlevel` of `AnalysisMethod`. It only works for which `signaltable` is `SampleDataTable`; otherwise, it falls back to use `pointlevel`.

`analytetable.txt` needs to contain analyte names, index of their internal standards, and index of of other analytes whose calibration curve is used. The column names are fixed for these three columns.
```
analytes isd   calibration other_information
analyte1 isd1  calibration_analyte_id1 other_information1
analyte2 isd2  calibration_analyte_id2 other_information2
⋮
```
The delimiter should be "\t", and the order of columns does not matter.

### *.at
It can contain multiple `*.sdt` or `*.adt`. The file names must start from an integer with `_` following the name, e.g. `0_area.sdt`. The integer is for the order of reading into `AnalysisTable`, and `name` will be the key. The name of signal data is determined in `method.am/config.txt`.

### Reading and writing Batch
To read a batch into julia, call `ChemistryQuantitativeAnalysis.read`.
```julia-repl
julia> batch = ChemistryQuantitativeAnalysis.read("batch_name.batch", T; table_type, analytetype, delim)
```
`T` is the sink function for tabular data; it should create an object following `Tables.jl` interface. `table_type` is `T` parameter in the type signature of `Batch` which determines the underlying table type, `analytetype` is a concrete type for `analyte`, `sampletype` is a concrete type for `sample`, and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist. 

For `analytetype` and `sampletype`, `string(cqaconvert(analytetype, x))` and `string(cqaconvert(sampletype, x))` should equal `x` if `x` is a valid string. Additionally, `tryparse` have to be extended for `CSV` parsing:
* `tryparse(::Type{analytetype}, x::AbstractString)` is neccessary for `AnalyteDataTable`.
* `tryparse(::Type{sampletype}, x::AbstractString)` is neccessary for `SampleDataTable`.

To write batch to disk, call `ChemistryQuantitativeAnalysis.write`. There is a keyword argument `delim` controling delimiter of tables.
```julia-repl
julia> ChemistryQuantitativeAnalysis.write("batch_name.batch", batch; delim = '\t')
```
There will be a folder `calibration` containing multiple `*.mcal` or `*.scal` folders. The former is for `MultipleCalibration` and the latter is for `SingleCalibration`.

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
conc = Float64[1, 2, 5, 10, 20, 50, 100]
signal1 = vcat(Float64[1, 2, 5, 10, 20, 50, 100], [1, 2, 5, 10, 20, 50, 100] .+ 0.1, [1, 2, 5, 10, 20, 50, 100] .- 0.1)
signal2 = vcat(Float64[1, 2, 5, 10, 20, 50, 100] .^ 2, [1, 2, 5, 10, 20, 50, 100] .^ 2 .+ 0.1, [1, 2, 5, 10, 20, 50, 100] .^ 2 .- 0.1)

# Create method
conctable = SampleDataTable(
   DataFrame(
         "level" => collect(1:7), 
         "G1(drug_a)" => conc,
         "G1(drug_b)" => conc .* 10), 
   :level; 
   analytetype = AnalyteTest
)
signaltable = SampleDataTable(
   DataFrame(
         "point" => reshape([string(a, "_", b) for (a, b) in Iterators.product(1:7, 1:3)], 21), 
         "level" => repeat(1:7, 3),
         "G1(drug_a)" => signal1,
         "G2(drug_a)" => repeat([5.0], 21),
         "G1(drug_b)" => signal2,
         "G2(drug_b)" => repeat([2.0], 21)), 
   :point; 
   analytetype = AnalyteTest
)
method = AnalysisMethod(conctable, signaltable, :area, :point; analyte = analytes, isd = [2, -1, 4, -1], calibration = [1, -1, 3, -1])

# Create sample data
cdata = AnalysisTable([:area], [
   SampleDataTable(
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
rdata = AnalysisTable([:area], [
   AnalyteDataTable(
         DataFrame(
            "Analyte" => analytes, 
            "S1" => Float64[6, 6, 200, 2],
            "S2" => Float64[24, 6, 800, 2],
            "S3" => Float64[54, 6, 9800, 2]
            ), 
         :Analyte
         )
   ]
) # Less efficient for quantification

# Create batch
cbatch = Batch(method, cdata)
rbatch = Batch(method, rdata)

# Calibration
cbatch.calibration # a vector of `MultipleCalibration`
cbatch.calibration[2].type = false # Use quadratic regression for the second analyte
rbatch.calibration[AnalyteG1("G1(drug_b)")].type = false # Use quadratic regression for `AnalyteG1("G1(drug_b)")`
update_calibration!(cbatch, 2)
update_calibration!(rbatch, AnalyteG1("G1(drug_b)"))

# Quantification
update_relative_signal!(cbatch) # A new data `cbatch.data.relative_signal` is created.
update_inv_predict!(cbatch) # Fit `cbatch.data.relative_signal` into calibration curve to create `cbatch.data.estimated_concentration`.
update_quantification!(cbatch) # equivalent to `update_inv_predict!(update_relative_signal!(cbatch))`

# Utils
analyteobj(cdata.area) # analytes of type `AnalyteTest`
sampleeobj(cdata) # samples of type `String`
analytename(cdata) # analytes of type `Symbol`
samplename(cdata.area) # samples of type `Symbol`
cdata.area[AnalyteTest("G2(drug_a)"), "S1"] = 6 # Set value using `dt[analyte, sample]`
cdata.area[AnalyteTest("G2(drug_a)"), "S1"] == 6
collect(eachanalyte(cdata.area))
collect(eachsample(cdata.area))
getanalyte(cdata.area, AnalyteG1("G1(drug_b)")) # get data of `AnalyteG1("G1(drug_b)")`
getanalyte(cdata.area, 1) # get data of first analyte
getsample(cdata.area, "S2")
dynamic_range(cbatch.calibration[1])
signal_range(rbatch.calibration[2])
lod(rbatch.calibration[2])
loq(rbatch.calibration[2])
formula_repr(cbatch.calibration[2])
```
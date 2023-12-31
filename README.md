# ChemistryQuantitativeAnalysis

[![Build Status](https://github.com/yufongpeng/ChemistryQuantitativeAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yufongpeng/ChemistryQuantitativeAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yufongpeng/ChemistryQuantitativeAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yufongpeng/ChemistryQuantitativeAnalysis.jl)

`ChemistryQuantitativeAnalysis.jl` is a package for quantitative analysis of chemicals based on tabular data.

## Tabular data wrapper
This package provides two basic wrappers, `ColumnDataTable{A, T}` and `RowDataTable{A, T}` which are subtypes of `AbstractDataTable{A, T}`. `ColumnDataTable` indicates that part of columns represent analytes, and all rows reprsent samples. `RowDataTable` indicates that part of columns represent samples, and all rows represent analytes. Both types have the same properties, but the actual meanings may be different. 
|Property|`ColumnDataTable{A, T}`|`RowDataTable{A, T}`|
|----------|---------------------|------------------|
|`analytename`|`Vector{Symbol}`, the column names that are analytes names|`Vector{Symbol}`, symbols transformed from column `analytecol`.|
|`samplename`|`Vector{Symbol}`, symbols transformed from column `samplecol`.|`Vector{Symbol}`, the column names that are sample names.|
|`analyte`|`Vector{A}`, analytes in user-defined types.|same|
|`sample`|`Vector`, the column `samplecol`.|`Vector{Symbol}`, the column names that are sample names..|
|`table`|Tabular data of type `T`|same|

`ColumnDataTable` can be created by `ColumnDataTable(table, samplecol; analytetype, analytename)`. By default, `analytename` includes all properties of `table` without `samplecol`. `RowDataTable` can be created by `RowDataTable(table, analytecol; analytetype, samplename)`. By default, `samplename` includes all properties of `table` without `analytecol`.
To add new samples to `ColumnDataTable{A, T}`, user can directly modify `table`; for `RowDataTable{A, T}`, user have to modify `samplename` as well. To add new analytes, user can directly modify `table` for `RowDataTable{A, T}`, and additionally modify `analyte` for `ColumnDataTable{A, T}`.

The package provides another two wrappers, `MethodTable{A, T}`, and `AnalysisTable{A, T} <: AbstractAnalysisTable{A, T}`.
### MethodTable
This type is used for storing method, containing all analytes, their internal standards and calibration curve setting, and data for fitting calibration curve.
|Property|Description|
|----------|---------|
|`analytetable`|`Table` with at least 3 columns, `analytes` identical to property `analytes`, `isd`, matching each analyte to index of its internal standard, and `calibration` matching each analyte to index of other analyte for fitting its calibration curve. `-1` indicates the analyte itself is internal standard, and `0` indicates no internal standard. For example, a row `(analytes = AnalyteX, isd = 2, calibration = 3)` means that internal standard of `AnalyteX` is the second analyte, and it will be quantified using calibration curve of the third analyte.|
|`signal`|`Symbol`, propertyname for extracting signal data from an `AnalysisTable`|
|`pointlevel`|`Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.|
|`conctable`|`AbstractDataTable{A, <: T}` containing concentration data for each level. Sample names must be symbol or string of integers for multiple levels. One level indicates using `SingleCalibration`.|
|`signaltable`|`AbstractDataTable{A, <: T}` containig signal for each point. It can be `nothing` if signal data is unecessary.|
|`analyte`|`Vector{A}`, analytes in user-defined types.|
|`isd`|`Vector{A}` that each analytes are internal standards.|
|`nonisd`|`Vector{A}` that each analytes are not internal standards.|same|
|`point`|`Vector`, calibration points, identical to `signaltable.samples`.|
|`level`|`Vector`, calibration levels, identical to `conctable.samples`.|

### AnalysisTable
`AnalysisTable{A, T}` is basically a `Dictionary{Symbol, <: AbstractDataTable{A, <: T}}` which data can be extracted using proeperty syntax. For example, `at.tables[:area] === at.area`.
|Property|Description|
|----------|---------|
|`analyte`|`Vector{A}`, analytes in user-defined types.|
|`sample`|`Vector`, sample names.|
|`tables`|`Dictionary{Symbol, <: AbstractDataTable{A, <: T}}`, a dictionary mapping data type to datatable.|
|Other|All keys of `tables`|

The key for signal data is determined by `method.signal`. Default names for relative signal, true concentration, estimated concentration and accuracy are `relative_signal`, `true_concentration`, `estimated_concentration` and `accuracy`.

## Calibration
This package provides two calibration types, `MultipleCalibration{A}` and `SingleCalibration{A}` which are subtypes of `AbstractCalibration{A}`.

### MultipleCalibration
This type fits and stores calibration curve. It can be created from a `MethodTable{A, T}` containing calibration data, an analyte `A` using function `calibration`.
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

To predict concentration, call `inv_predict` or `inv_predict!` for inplace replacement of `table.x̂`. To calculate accuracy, call `accuracy` or `acccuracy!` for inplace replacement of `table.accuracy`. `inv_predict_accuracy!` calls `inv_predict!` and `acccuracy!` subsequently. `type`, `zero`, and `weigtht` can be modified directly. To change internal standard, modify `analyte`. After any modification, call `update_calibration!` with method to update the `model`.

### SingleCalibration
This type contains data for single pont calibration. 
|Field|Description|
|----------|-----------|
|`analyte`|`Tuple{A}` is the analyte with known concentration (internal standard).|
|`conc`|`Float64`, concentration of analyte.|

## Batch
`Batch{A, T}` represents a batch for quantitative analysis where `A` is analyte type and `T` is table type.
|Property|Description|
|----------|-----------|
|`method`|`MethodTable{A, <: T}`, method.|
|`calibration`|`Vector{MultipleCalibration{<: A}}` or `Vector{SingleCalibration{<: A}}`|
|`data`|Data for analysis, `AnalysisTable{A, <: T}` or `Nothing`.|
|`analyte`|`Vector{A}`, analytes in user-defined types, identical to `method.analytetable.analyte`.|
|`isd`|`Vector{<: A}`, analytes which are internal standards.|
|`nonisd`|`Vector{A}` that each analytes are not internal standards.|same|
|`point`|`Vector{Symbol}`, calibration points, identical to `method.signaltable.samples`.|
|`level`|`Vector{Symbol}`, calibration levels, identical to `method.conctable.samples`.|

It can be created with only `method` and optionally `data`.

To predict concentration and calculate accuracy for calibration points, call `inv_predict!` and `acccuracy!`, respectively. To calculate relative signal, concentration or accuracy and save the result, call `update_relative_signal!`, `update_inv_predict!` (in combination, `update_quantification!`) and `update_accuracy!`, respectively. `inv_predict_accuracy!` calls `inv_predict_cal!` and `acccuracy!` subsequently. To change internal standard, call `set_isd!` with object `analyte` and `isd`.

## Reading and writting data to disk
To use data on disk, user should create a directory in the following structure:
```
batch_name.batch
├──config.txt
├──method.mt
│  ├──true_concentration.dt
│  │  ├──config.txt
│  │  └──table.txt
│  ├──area.dt
│  │  ├──config.txt
│  │  └──table.txt
│  ├──analytetable.txt
│  └──config.txt
├──calibration
│  ├──1.mcal
│  └──2.mcal
└──data.at
   ├──0_quantity1.dt
   ├──1_quantity2.dt 
   └──2_quantity3.dt
```
Config files have the following general forms
```
[property1]
value

[property2]
value1
value2
value3
.
.
.
```
The property `delim` determines the default delimiter for `table.txt` in this directory and subdirectories.

`data.at` and `calibration` is not necessary for initializing a batch. The former can be added to the batch directly in julia, and the latter will be generated after calibration.

### *.dt
All `*.dt` files will be read as `ColumnDataTable` or `RowDataTable`. They contain `config.txt` and `table.txt`.

Config file for `ColumnDataTable` needs the following properties.
```
[Type]
C

[delim]
\t

[Sample]
sample_col_name

[Analyte]
analyte_col_name_1
analyte_col_name_2
.
.
.
``` 
Config file for `RowDataTable` needs the following properties.
```
[Type]
R

[delim]
\t

[Analyte]
analyte_col_name

[Sample]
sample_col_name_1
sample_col_name_2
.
.
.
``` 

### *.mt
It must contain two `*dt` files. `true_concentration.dt` contains true concentration for each analyte and level. The sample names must be integers.
Another `*.dt` file is signal data for each analyte and calibration point. The file name is determined by `config.txt`.

Config file for `method.mt` needs the following  properties.
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
.
.
.
```
`signal` specifys which `.dt` file serving as signal data. For the above file, `method.mt/area.dt` will become `method.signaltable`.

`pointlevel` maps each point to level which should be integers.

`level` specifys the column representing property `pointlevel` of `MethodTable`. It only works for which `signaltable` is `ColumnDataTable`; otherwise, it falls back to use `pointlevel`.

`analytetable.txt` needs to contain analyte names, index of their internal standards, and index of of other analytes whose calibration curve is used. 
```
analytes isd   calibration
analyte1 isd1  calibration_analyte_id1
analyte2 isd2  calibration_analyte_id2
.
.
.
```
The delimiter should be "\t".

### *.at
It can contain multiple `*.dt`. The file names must start from an integer, `_` and `name.dt`, e.g. `0_area.dt`. The integer is for the order of reading into `AnalysisTable`, and `name` will be the key. The name of signal data is determined in `method.mt/config.txt`.

### Reading and writing Batch
To read a batch into julia, call `ChemistryQuantitativeAnalysis.read`.
```julia-repl
julia> batch = ChemistryQuantitativeAnalysis.read("batch_name.batch", T; table_type, analytetype, delim)
```
`T` is the sink function for tabular data; it should create an object following `Tables.jl` interface. `table_type` is `T` parameter in the type signature of `Batch` which determines the underlying table type, `analytetype` is a concrete type for `analyte` which msut have a method for string input, and `delim` specifies delimiter for tabular data if `config[:delim]` does not exist.

To write batch to disk, call `ChemistryQuantitativeAnalysis.write`. There is a keyword argument `delim` controling delimiter of tables.
```julia-repl
julia> ChemistryQuantitativeAnalysis.write("batch_name.batch", batch; delim = '\t')
```
There will be a folder `calibration` containing multiple `*.mcal` or `*.scal` folders. The former is for `MultipleCalibration` and the latter is for `SingleCalibration`.

## Examples
```julia
using ChemistryQuantitativeAnalysis, TypedTables, DataFrames
import Base: show

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
const AnalyteTest = Union{AnalyteG1, AnalyteG2, AnalyteOther}
show(io::IO, analyte::AnalyteTest) = print(io, analyte.name)

# Analyte parser
function Union{AnalyteG1, AnalyteG2, AnalyteOther}(name::String)
    g = match(r"^G(\d)\(.*\)$", name)
    isnothing(g) && return AnalyteOther(name)
    g = parse(Int, first(g))
    g == 1 ? AnalyteG1(name) : g == 2 ? AnalyteG2(name) : AnalyteOther(name)
end

# Generate data
analyte_names = ["G1(drug_a)", "G2(drug_a)", "G1(drug_b)", "G2(drug_b)"]
conc = Float64[1, 2, 5, 10, 20, 50, 100]
signal1 = vcat(Float64[1, 2, 5, 10, 20, 50, 100], [1, 2, 5, 10, 20, 50, 100] .+ 0.1, [1, 2, 5, 10, 20, 50, 100] .- 0.1)
signal2 = vcat(Float64[1, 2, 5, 10, 20, 50, 100] .^ 2, [1, 2, 5, 10, 20, 50, 100] .^ 2 .+ 0.1, [1, 2, 5, 10, 20, 50, 100] .^ 2 .- 0.1)

# Create method
conctable = ColumnDataTable(
   DataFrame(
         "level" => collect(1:7), 
         "G1(drug_a)" => conc,
         "G1(drug_b)" => conc .* 10), 
   :level; 
   analytetype = AnalyteTest
)
signaltable = ColumnDataTable(
   DataFrame(
         "point" => repeat(1:7, 3), 
         "G1(drug_a)" => signal1,
         "G2(drug_a)" => repeat([5.0], 21),
         "G1(drug_b)" => signal2,
         "G2(drug_b)" => repeat([2.0], 21)), 
   :point; 
   analytetype = AnalyteTest
)
method = MethodTable(conctable, signaltable, :area, :point; analyte = AnalyteTest.(analyte_names), isd = [2, -1, 4, -1], calibration = [1, -1, 3, -1])

# Create sample data
rdata = AnalysisTable([:area], [
   RowDataTable(
         DataFrame(
            "Analyte" => analyte_names, 
            "S1" => Float64[6, 6, 200, 2],
            "S2" => Float64[24, 6, 800, 2],
            "S3" => Float64[54, 6, 9800, 2]
            ), 
         :Analyte; 
         analytetype = AnalyteTest
         )
   ]
)
cdata = AnalysisTable([:area], [
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

# Create batch
cbatch = Batch(method, cdata)
rbatch = Batch(method, rdata)

# Calibration
cbatch.calibration # a vector of `MultipleCalibration`
cbatch.calibration[2].type = false # Use quadratic regression for the second analyte
rbatch.calibration[AnalyteG1("G1(drug_b)")].type = false # Use quadratic regression for AnalyteG1("G1(drug_b)")
update_calibration!(cbatch, 2)
update_calibration!(rbatch, AnalyteG1("G1(drug_b)"))

# Quantification
update_relative_signal!(cbatch) # A new data `cbatch.data.relative_signal` is created.
update_inv_predict!(cbatch) # Fit `cbatch.data.relative_signal` into calibration curve to create `cbatch.data.estimated_concentration`.
update_quantification!(cbatch) # equivalent to `update_inv_predict!(update_relative_signal!(cbatch))`

# Utils
cdata.area[1, "G2(drug_a)"] = 6
cdata.area[1, "G1(drug_a)"] == 6
collect(eachanalyte(cdata.area))
collect(eachsample(cdata.area))
getanalyte(cdata.area, AnalyteG1("G1(drug_b)"))
getanalyte(cdata.area, 1)
getsample(cdata.area, "S2")
dynamic_range(cbatch.calibration[1])
signal_range(rbatch.calibration[2])
formula_repr(cbatch.calibration[2])
```
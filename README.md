# QuantitativeAnalysis

[![Build Status](https://github.com/yufongpeng/QuantitativeAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yufongpeng/QuantitativeAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yufongpeng/QuantitativeAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yufongpeng/QuantitativeAnalysis.jl)

`QuantitativeAnalysis.jl` is a package for instrument calibration and analyte quantification based on tabular data.

## Tabular data wrapper
This package provides two basic wrappers, `ColumnDataTable{A, T}` and `RowDataTable{A, T}` which are subtypes of `AbstractDataTable{A, T}`. `ColumnDataTable` indicates that part of columns represent analytes, and all rows reprsent samples. `RowDataTable` indicates that part of columns represent samples, and all rows represent analytes. Both types have the same properties, but the actual meanings may be different. 
|Property|`ColumnDataTable{A, T}`|`RowDataTable{A, T}`|
|----------|---------------------|------------------|
|`sample_name`|`Symbol`, the column name that each element is sample name.|`Vector{Symbol}`, the column names that are sample names.|
|`analyte_name`|`Vector{Symbol}` stored in field `config`, the column names that are analytes names|`Symbol`, the column name that each element is analyte name.|
|`samples`|`Vector{Symbol}`, the column names that are sample names.|`Vector{Symbol}`, symbols transformed from column `sample_name`.|
|`analytes`|`Vector{A}` stored in field `config`, analytes in user-defined types.|`Vector{A}`, analytes in user-defined types.|
|`table`|Tabular data of type `T`|same|

`ColumnDataTable` can be created by `ColumnDataTable(table, sample_name; analyte_fn, analyte_name)`. By default, `analyte_name` includes all properties of `table` without `sample_name`. `RowDataTable` can be created by `RowDataTable(table, analyte_name; analyte_fn, sample_name)`. By default, `analyte_name` includes all properties of `table` without `sample_name`.
To add new samples to `ColumnDataTable{A, T}`, user can directly modify `table`; for `RowDataTable{A, T}`, user have to modify `sample_name` as well. To add new analytes, user can directly modify `table`, and modify `config` for `ColumnDataTable{A, T}` (`config` is a `TypedTables.Table`) and `analytes` for `RowDataTable{A, T}`.

The package provides another two wrappers, `MethodTable{A, T}`, and `AnalysisTable{A, T} <: AbstractAnalysisTable{A, T}`.
### MethodTable
This type is used for storing method, containing all analytes, their internal standards and calibration curve setting, and data for fitting calibration curve.
|Property|Description|
|----------|---------|
|`signal`|`Symbol`, propertyname for extracting signal data from an `AnalysisTable`|
|`analyte_map`|`Table` with 3 columns, `analytes` identical to property `analytes`, `isd`, matching each analyte to index of its internal standard, and `calibration` matching each analyte to index of other analyte for fitting its calibration curve. `-1` indicates the analyte itself is internal standard, and `0` indicates no internal standard. For example, a row `(analytes = AnalyteX, isd = 2, calibration = 3)` means that internal standard of `AnalyteX` is the second analyte, and it will be quantified using calibration curve of the third analyte.|
|`level_map`|`Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.|
|`conctable`|`AbstractDataTable{A, <: T}` containing concentration data for each level. Sample names must be symbol or string of integers for multiple levels. One level indicates using `SingleCalibration`.|
|`signaltable`|`AbstractDataTable{A, <: T}` containig signal for each point. It can be `nothing` if signal data is unecessary.|
|`analytes`|`Vector{A}`, analytes in user-defined types.|
|`isds`|`Vector{A}` that each analytes are internal standards.|
|`nonisds`|`Vector{A}` that each analytes are not internal standards.|same|
|`points`|`Vector{Symbol}`, calibration points, identical to `signaltable.samples`.|
|`levels`|`Vector{Symbol}`, calibration levels, identical to `conctable.samples`.|

### AnalysisTable
`AnalysisTable{A, T}` is basically a `Dictionary{Symbol, <: AbstractDataTable{A, <: T}}` which data can be extracted using proeperty syntax. For example, `at.tables[:area] === at.area`.
|Property|Description|
|----------|---------|
|`analytes`|`Vector{A}`, analytes in user-defined types.|
|`samples`|`Vector{Symbol}`, sample names.|
|`tables`|`Dictionary{Symbol, <: AbstractDataTable{A, <: T}}`, a dictionary mapping data type to datatable.|
|Other|All keys of `tables`|

The key for signal data is determined by `method.signal`. Default names for relative signal, true concentration, estimated concentration and accuracy are `relative_signal`, `true_concentration`, `estimated_concentration` and `accuracty`.

## Calibration
This package provides two calibration types, `MultipleCalibration{A}` and `SingleCalibration{A}` which are subtypes of `AbstractCalibration{A}`.

### MultipleCalibration
This type fits and stores calibration curve. It can be created from a `MethodTable{A, T}` containing calibration data, an analyte `A`.
|Attributes|Description|
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

To fit a new model, call `calfit` or `calfit!` for inplace substitution of `model`. To predict concentration, call `inv_predict` or `inv_predict!` for inplace replacement of `table.x̂`. To calculate accuracy, call `accuracy` or `acccuracy!` for inplace replacement of `table.accuracy`. `inv_predict_accuracy!` calls `inv_predict!` and `acccuracy!` subsequently. `type`, `zero`, and `weigtht` can be modified directly. To change internal standard, modify `analyte` and call `update_calibration!` with method (It is not recommended to change internal standard in this way, see `set_isd!` in [`Batch`](#batch) section). 

### SingleCalibration
This type contains data for single pont calibration. 
|Attributes|Description|
|----------|-----------|
|`analyte`|`Tuple{A}` is the analyte with known concentration (internal standard).|
|`conc`|`Float64`, concentration of analyte.|

## Batch
`Batch{A, T}` represents a batch for quantitative analysis where `A` is analyte type and `T` is table type.
|Attributes|Description|
|----------|-----------|
|`method`|`MethodTable{A, <: T}`, method.|
|`calibration`|`Vector{MultipleCalibration{<: A}}` or `Vector{SingleCalibration{<: A}}`|
|`data`|Data for analysis, `AnalysisTable{A, <: T}` or `Nothing`.|

It can be created with only `method` and optionally `data`.

To fit new models for all `calibration`, call `calfit!`. To predict concentration and calculate accuracy for calibration points, call `inv_predict!` and `acccuracy!`, respectively. To calculate relative signal, concentration or accuracy and save the result, call `update_relative_signal!`, `update_inv_predict!` (in combination, `update_quantification!`) and `update_accuracy!`, respectively. `inv_predict_accuracy!` calls `inv_predict_cal!` and `acccuracy!` subsequently. To change internal standard, call `set_isd!` with object `analyte` and `isd`.

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
│  ├──analyte_map.txt
│  └──config.txt
├──calibration
│  ├──1.mcal
│  └──2.mcal
└──data.at
   ├──0_quantity1.dt
   ├──1_quantity2.dt 
   └──2_quantity3.dt
```
Config files has the following general forms
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
The config file in `batch_name.batch` must contain property `delim` which determines the default delimiter for all `table.txt` in this directory and subdirectories.

`data.at` and `calibration` is not necessary for initializing a batch. The former can be added to the batch directly in julia, and the latter will be generated after calibration.

### *.dt
All `*.dt` files will be read as `ColumnDataTable` or `RowDataTable`. They contain `config.txt` and `table.txt`.

Config file for `ColumnDataTable` needs at least the following two properties.
```
[Type]
C

[Sample]
sample_col_name

[Analyte]
analyte_col_name_1
analyte_col_name_2
.
.
.
``` 
Config file for `RowDataTable` needs at least the following two properties.
```
[Type]
R

[Analyte]
analyte_col_name

[Sample]
sample_col_name_1
sample_col_name_2
.
.
.
``` 
User can also provide `delim` to overwrite the default one defined in `batch_name.batch/config.txt`.

### *.mt
It must contain two `*dt` files. `true_concentrstion.dt` contains true concentration for each analyte and level. The sample names must be integers.
Another `*.dt` file is signal data for each analyte and calibration point. The file name is determined by `config.txt`.

Config file for `method.mt` needs two properties, `signal` and `level_map`.
```
[signal]
area

[level_map]
level_for_1st_point
level_for_2nd_point
.
.
.
```
`signal` specify `.dt` file serving as signal data. For the above file, `method.mt/area.dt` will be `method.signaltable`.

`level_map` maps each point to level which should be integera as well.

`analyte_map.txt` contains analyte names, index of their internal standards, and index of of other analytes whose calibration curve is used. 
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
To read a batch into julia, call `QuantitativeAnalysis.read`.
```julia-repl
julia> QuantitativeAnalysis.read("batch_name.batch", T; table_type, analyte_fn)
```
`T` is the sink function for tabular data; it should create an object following `Tables.jl` interface. `table_type` is `T` parameter in the type signature of `Batch` which determines the underlying table type, and `analyte_fn` is responsible for converting `analyte` in string type into user defined type `A`.

To write project to disk, call `QuantitativeAnalysis.write`. There is a keyword argument `delim` controling whether using "\t" or "," as delim for tables.
```julia-repl
julia> QuantitativeAnalysis.write("batch_name.pjc", batch; delim = "\t")
```
A new folder `calibration` containing multiple `*.mcal` or `*.scal` folders. The former is for `MultipleCalibration` and the latter is for `SingleCalibration`.
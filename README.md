# Calibration

[![Build Status](https://github.com/yufongpeng/Calibration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yufongpeng/Calibration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yufongpeng/Calibration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yufongpeng/Calibration.jl)

`Calibration.jl` is a package for instrument calibration and analyte quantification based on tabular data.

## Tabular data wrapper
This package provides two basic wrappers, `ColumnAnalysisTable{A, T}` and `RowAnalysisTable{A, T}` which are subtypes of `AbstractAnalysisTable{A, T}`. `ColumnAnalysisTable` indicates that part of columns represent analytes, and all rows reprsent samples. `RowAnalysisTable` indicates that part of columns represent samples, and all rows represent analytes. Both types have the same properties, but the actual meanings may be different. 
|Property|`ColumnAnalysisTable{A, T}`|`RowAnalysisTable{A, T}`|
|----------|---------------------|------------------|
|`sample_name`|`Symbol`, the column name that each element is sample name.|`Vector{Symbol}`, the column names that are sample names.|
|`analyte_name`|`Vector{Symbol}` stored in field `config`, the column names that are analytes names|`Symbol`, the column name that each element is analyte name.|
|`analytes`|`Vector{A}` stored in field `config`, analytes in user-defined types. Its length and correponding analytes must matches those of `isd_map`.|same|
|`isd_map`|`Vector{Int}` stored in field `config`, each element is the index of the correponding internal standrads in `analytes` of the analyte in this position. `0` means no internal standard, and `-1` means the analyte itself is internal standard. For example, `[2, -1]` means that internal standard of the first analyte is the second analyte, and the second analyte is an internal standrad.|`Symbol` or `Nothing`, the column name that contains internal standard information.|
|`table`|Tabular data of type `T`|same|
|`samples`|`Vector{Symbol}`, the column names that are sample names.|same|
|`isds`|`Vector{A}` that each analytes are internal standards.|same|
|`nonisds`|`Vector{A}` that each analytes are not internal standards.|same|

To add new samples to `ColumnAnalysisTable{A, T}`, user can directly modify `table`; for `RowAnalysisTable{A, T}`, user have to modify `sample_name` as well. To add new analytes, user can directly modify `table`, and modify `config` for `ColumnAnalysisTable{A, T}` (`config` is a `TypedTables.Table`) and `analytes` for `RowAnalysisTable{A, T}`.

The package provides another two wrappers, `SampleWrapperTable{A, T}`, and `CalWrapperTable{A, T}` which are subtypes of `AbstractWrapperTable{A, T}`. 
### SampleWrapperTable
This type wraps a `Vector{Int}` as attribute `cal_map` and an `AbstractAnalysisTable` as attribute `analysistable`. The length and corresponding analytes of `cal_map` matches analyte in `analysistable.analytes`, and each element is the index of another analyte that its calibration curve is used for quantification of this analyte. For example, `[1, -1, 1]` means that the first analyte uses calibration curve of itself, the second analyte is internal standard, and the third analyte also uses calibration curve of the first analyte. All properties for an `AbstractAnalysisTable` are also available for `SampleWrapperTable`.
### CalWrapperTable
This type has three attribute `level_map`, `conctable`,and `signaltable`. `conctable` is an `AbstractAnalysisTable{A, T}` containing concentration for each level, and ISD information. Sample names must be symbol or string of integers for multiple levels. With only one level, `level_map` and `signaltable` are ignored. With multiple levels, `signaltable` should be an `AbstractAnalysisTable{A}` containig signal for each points, and `level_map` will be a `Vector{Int}` matching each point to level.

## Calibration
This package provides two calibration types, `MultipleCalibration{A, T}` and `SingleCalibration{A, T}` which are subtypes of `AbstractCalibration{A, T}`.

### MultipleCalibration
This type fits and stores calibration curve. It can be created from a `AbstractAnalysisTable{A, T}` containing calibration data, an analyte `A`.
|Attributes|Description|
|----------|-----------|
|`analyte`|`Int` indicates the index the analyte in `caltable.conctable.analytes`|
|`isd`|`Int` indicates the index internal standard in `caltable.conctable.analytes`. `0` indicates no internal standard.|
|`type`|`Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`).|
|`zero`|`Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`).|
|`weight`|`Float64` represents the exponential applying to each element of `x` as a weighting vector.|
|`formula`|`FormulaTerm`, the formula for fitting calibration curve.|
|`caltable`|`CalWrapperTable{A, T}`, the original calibration data.|
|`table`|`TypedTable.Table`, the clean up calibration data, containing 7 columns.|
|`model`|`GLM` object|

The columns in `table`:
|Column|Description|
|------|-----------|
|`id`|Sample name|
|`level`|The index of concentration level. The larger, the higher concentraion it represents.|
|`y`|Signal|
|`x`|Concentraion|
|`x̂`|Predicted concentration|
|`accuracy`|Accuracy, i.e. `x̂/x`.|
|`include`|Whether this point is included or not|

To fit a new model, call `calfit` or `calfit!` for inplace substitution of `model`. To predict concentration, call `inv_predict` or `inv_predict!` for inplace replacement of `table.x̂`. To calculate accuracy, call `accuracy` or `acccuracy!` for inplace replacement of `table.accuracy`. `inv_predict_accuracy!` calls `inv_predict!` and `acccuracy!` subsequently. `type`, `zero`, and `weigtht` can be modified directly. To switch internal standard, modified `isd_map` of `source` and recreate the object (It is not recommended to switch internal standard in this way, see `switch_isd!` in [`Project`](#project) section). 

### SingleCalibration
This type contains data for single pont calibration. 
|Attributes|Description|
|----------|-----------|
|`analyte`|`Int` indicates the index the analyte in `caltable.conctable.analytes`|
|`isd`|`Int` indicates the index internal standard in `caltable.conctable.analytes`. `0` indicates no internal standard.|
|`caltable`|`CalWrapperTable{A, T}`, the original calibration data.|

## Project
A type for holding all data. 
|Attributes|Description|
|----------|-----------|
|`calibration`|`Vector{MultipleCalibration{A, T}}` or `Vector{SingleCalibration{A, T}}`|
|`caltable`|`CalWrapperTable{A, T}`, the original calibration data.|
|`sampletable`|Sample data, `SampleWrapperTable{A}` or `Nothing`.|
|`resulttable`|Result data, an `SampleWrapperTable{A}` or `nothing`.|

It can be created from function `project` by providing `caltable`, and optionally `sampletable`, `resulttable`.

To fit new models for all `calibration`, call `calfit!`. To predict concentration, call `inv_predict_cal!` for `caltable` and `inv_predict_sample!` for `sampletable`. To calculate accuracy, call `accuracy` or `acccuracy!` for inplace replacement. `inv_predict_accuracy!` calls `inv_predict_cal!` and `acccuracy!` subsequently. To switch internal standard, call `switch_isd!` with object `analyte` and id of the new internal standard.

## Reading and writting data to disk
To use data on disk, user should create a directory in the following structure:
```
project_name.pjc
├──config.txt
├──cal.ctbl
│  ├──conc.tbl
│  │  ├──config.txt
│  │  └──table.txt
│  ├──signal.tbl
│  │  ├──config.txt
│  │  └──table.txt
│  └──level_map.txt
└──sample.tbl
   ├──cal_map.txt
   ├──config.txt 
   └──table.txt
```
The config file in `project_name.pjc` should be the following form:
```
[delim]
\t
```
or 
```
[delim]
,
```
It indicates that all `table.txt` in this directory and subdirectories should use "\t" or "," as delim.

Other config files for `.tbl` have the following two forms:
For `ColumnAnalysisTable`
```
[Sample]
sample_col_name

[Analyte]
analyte_col_name_1  isd1
analyte_col_name_2  isd2
.
.
.
``` 
`analyte_col_name_x` and `isdx` should be separated by tab, and `isdx` should be an integer or empty which will be parsed as the property `isd_map`. Empty string will be parsed as `0`.
For `RowAnalysisTable`
```
[ISD]
isd_col_name

[Analyte]
analyte_col_name

[Sample]
sample_col_name_1
sample_col_name_2
.
.
.
``` 
If `isd_col_name` is not in the columns, the property `isd_map` will be `nothing`, and user cannot add `isd_map` afterwards (as `RowAnalysisTable` is immutable). To be noticed, `conctable` must contain the correct internal standard information.
`[Analyte]`, `[Sample]` and `[ISD]` can be in any orders.

`level_map` matches points to levels:
```
level_for_1st_point
level_for_2nd_point
.
.
.
```
This file as well as `signal.tbl` will be ignored if there is only one level in `conc.tbl` indicating this project will use `SingleCalibration`.

`sample.tbl` is optional; user can add the sample data afterwards.
`cal_map` matches analytes to analytes whose calibration curve is used.
```
calibration_analyte_id_for_1st_analyte
calibration_analyte_id_for_2nd_analyte
.
.
.
```
All elements should be integers.
To read the project into julia, call `Calibration.read`.
```julia-repl
julia> Calibration.read("project_name.pjc", T)
```
`T` is sink for tabular data; it should create an object following `Tables.jl` interface.

To write project to disk, call `Calibration.write`. There is a keyword argument `delim` controling whether using "\t" or "," as delim for tables.
```julia-repl
julia> Calibration.write("project_name.pjc", project; delim = "\t")
```
There will be a new folder `calibration` containing multiple `.mcal` or `.scal` folders. The former is for `MultipleCalibration` and the latter is for `SingleCalibration`. Prediction result will be stored as `result.tbl`.
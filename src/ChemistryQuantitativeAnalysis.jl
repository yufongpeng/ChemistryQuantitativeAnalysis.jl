module ChemistryQuantitativeAnalysis

using GLM, CSV, TypedTables, LinearAlgebra, Dictionaries, ThreadsX, Tables, Pkg
export SampleDataTable, AnalyteDataTable, 
    AnalysisTable, analysistable, AnalysisMethod, Batch, 
    MultipleCalibration, SingleCalibration, calibration, init_calibration!, update_calibration!,
    analyteobj, sampleobj, analytename, samplename, analytecol, samplecol,
    inv_predict, set_inv_predict, set_inv_predict!, update_inv_predict!,
    relative_signal, set_relative_signal, set_relative_signal!, update_relative_signal!,
    quantification, set_quantification, set_quantification!, update_quantification!,
    accuracy, set_accuracy, set_accuracy!, update_accuracy!,
    isdof, isisd, 
    findanalyte, getanalyte, findsample, getsample, eachanalyte, eachsample,
    dynamic_range, lloq, uloq, signal_range, signal_lloq, signal_uloq, lod, loq, 
    formula_repr, weight_repr, weight_value, formula_repr_ascii, weight_repr_ascii, format_number, mkbatch, 
    typedmap

import Base: getproperty, propertynames, show, write, eltype, length, iterate, 
        getindex, setindex!, insert!, get!, delete!, get, 
        pairs, keys, values, haskey, isassigned, 
        copy
import Dictionaries: set!, unset!, isinsertable, issettable
import Tables: istable, rowaccess, rows, columnaccess, columns
    
abstract type AbstractCalibration{A, N} end
abstract type AbstractDataTable{A, S, N, T} end
# abstract type AbstractAnalysisTable{A, T, S} end
const TypeOrFn = Union{<: Type, <: Function}

"""
    SampleDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}

Tabular data wrapper indicates part of columns represent analytes, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, `N` determines numeric value type, and `T` determines table type.

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types. Use `analyteobj` to get this field.
* `sample`: `Vector{S}`, samples in user-defined types, i.e., `getproperty(table(dt), samplecol(dt))`. Use `sampleobj` to get this field.
* `samplecol`: `Symbol`, the column name that each element is sample name. Use `samplecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct SampleDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}
    analyte::Vector{A}
    sample::Vector{S}
    samplecol::Symbol
    table::T
    function SampleDataTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, samplecol::Symbol, table::T) where {A, S, T}
        allunique(sample) || throw("`sample` should have unique elemnts.")
        tp = collect(map(eltype, columns(table)))
        ps = propertynames(table)
        idp = findfirst(==(Symbol(first(analyte))), ps)
        isnothing(idp) && throw(ArgumentError("`analyte[1]` is not in the table"))
        N = tp[idp]
        for i in 2:lastindex(analyte)
            idp = findfirst(==(Symbol(analyte[i])), ps)
            isnothing(idp) && throw(ArgumentError("`analyte[$i]` is not in the table"))
            N == tp[idp] || throw(ArgumentError(string("The element type of `analyte[$i]` is `", tp[idp], "`; convert it to `", N, "`.")))
        end
        new{A, S, N, T}(convert(Vector{A}, analyte), convert(Vector{S}, sample), samplecol, table)
    end
end
"""
    SampleDataTable(analyte::AbstractVector, samplecol::Symbol, table)
    SampleDataTable(table, analytename::AbstractVector, samplecol::Symbol)

An interface equivalent to `SampleDataTable(analyte, getproperty(table, samplecol), samplecol, table)`.
"""
SampleDataTable(analyte::AbstractVector, samplecol::Symbol, table) = 
    SampleDataTable(analyte, getproperty(table, samplecol), samplecol, table)
SampleDataTable(table, analyte::AbstractVector, samplecol::Symbol) = 
    SampleDataTable(analyte, samplecol, table)
"""
    SampleDataTable(analytetype::TypeOrFn, samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))
    SampleDataTable(table, analytetype::TypeOrFn, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))
    SampleDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol]))
    SampleDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol]))
    

Higher level interfaces for `SampleDataTable{A}`.

* `analytename`: `AbstractVector`, the column names of `table` that are analyte names. It will be converted to `AbstractVector{String}` before conversion to `AbstractVector{A}` (`string.(analytename)`).
* `samplecol`: `Symbol`, the column name that each element is sample name.
* `analytetype`: type `A` or a function. After first conversion of `analytename`, it convert `analytename` to `AbstractVector{A}` using `cqamap(analytetype, analytename)`. See `cqaconvert` for the requirement of `analytetype`.
"""
SampleDataTable(analytetype::TypeOrFn, samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(cqamap(analytetype, string.(analytename)), samplecol, table)
SampleDataTable(table, analytetype::TypeOrFn, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(analytetype, samplecol, table; analytename)
SampleDataTable(samplecol::Symbol, table; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(string.(analytename), samplecol, table)
SampleDataTable(table, samplecol::Symbol; analytename = setdiff(propertynames(table), [samplecol])) = 
    SampleDataTable(samplecol, table; analytename)
"""
    AnalyteDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}

Tabular data wrapper indicates part of columns represent analyte, and all rows reprsent samples. `A` determines analyte type, `S` determines sample type, `N` determines numeric value type, and `T` determines table type.

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types, i.e., `getproperty(table(dt), analytecol(dt))`. Use `analyteobj` to get this field.
* `sample`: `Vector{S}`, samples in user-defined types. Use `sampleobj` to get this field.
* `analytecol`: `Symbol`, the column name that each element is analyte name. Use `analytecol` to get this field.
* `table`: Tabular data of type `T`. Use `table` to get this field.

# Properties
All properties of `table`.
"""
struct AnalyteDataTable{A, S, N <: Real, T} <: AbstractDataTable{A, S, N, T}
    analyte::Vector{A}
    sample::Vector{S}
    analytecol::Symbol
    table::T
    function AnalyteDataTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, analytecol::Symbol, table::T) where {A, S, T}
        allunique(analyte) || throw("`analyte` should have unique elemnts.")
        tp = collect(map(eltype, columns(table)))
        ps = propertynames(table)
        idp = findfirst(==(Symbol(first(sample))), ps)
        isnothing(idp) && throw(ArgumentError("`sample[1]` is not in the table"))
        N = tp[idp]
        for i in 2:lastindex(sample)
            idp = findfirst(==(Symbol(sample[i])), ps)
            isnothing(idp) && throw(ArgumentError("`sample[$i]` is not in the table"))
            N == tp[idp] || throw(ArgumentError(string("The element type of `sample[$i]` is `", tp[idp], "`; convert it to `", N, "`.")))
        end
        new{A, S, N, T}(convert(Vector{A}, analyte), convert(Vector{S}, sample), analytecol, table)
    end
end
"""
    AnalyteDataTable(analytecol::Symbol, samplename::AbstractVector, table)
    AnalyteDataTable(table, analytecol::Symbol, samplename::AbstractVector)
    
An interface equivalent to `AnalyteDataTable(getproperty(table, analytecol), samplename, analytecol, table)`.
"""
AnalyteDataTable(analytecol::Symbol, samplename::AbstractVector, table) = 
    AnalyteDataTable(getproperty(table, analytecol), samplename, analytecol, table)
AnalyteDataTable(table, analytecol::Symbol, samplename::AbstractVector) = 
    AnalyteDataTable(analytecol, samplename, table)
"""
    AnalyteDataTable(analytecol::Symbol, sampletype::TypeOrFn, table; samplename = setdiff(propertynames(table), [analytecol]))
    AnalyteDataTable(table, analytecol::Symbol, sampletype::TypeOrFn; samplename = setdiff(propertynames(table), [analytecol]))
    AnalyteDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol]))
    AnalyteDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol]))

Higher level interfaces for `AnalyteDataTable{A, S}`.

* `analytecol`: `Symbol`, the column name that each element is analyte name.
* `samplename`: `AbstractVector`, the column names of `table` that are sample names. It will be converted to `Vector{String}` before conversion to `AbstractVector{S}` (`string.(samplename)`).
* `sampletype`: type `S` or a function. After first conversion of `sampletype`, it converts `samplename` to `AbstractVector{S}` using `cqamap(sampletype, sampletype)`. See `cqaconvert` for the requirement of `sampletype`.
"""
AnalyteDataTable(analytecol::Symbol, sampletype::TypeOrFn, table; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, cqamap(sampletype, string.(samplename)), table)
AnalyteDataTable(table, analytecol::Symbol, sampletype::TypeOrFn; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, sampletype, table; samplename)
AnalyteDataTable(analytecol::Symbol, table; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, string.(samplename), table)
AnalyteDataTable(table, analytecol::Symbol; samplename = setdiff(propertynames(table), [analytecol])) = 
    AnalyteDataTable(analytecol, table; samplename)

"""
    AnalysisTable{A, S, T <: AbstractDataTable{A, S}}

Wrapper of multiple tables representing different types of values. `A` determines analyte type, `S` determines sample type, and `T` determines datatable type. It is implemented as a `Dictionary{Symbol, T <: AbstractDataTable{A, S}}`, but unlike dictionaries, each datatable can also be extracted using `getproperty`. 

# Fields
* `analyte`: `Vector{A}`, analytes in user-defined types. Use `analyteobj` to get this field.
* `sample`: `Vector{S}`, samples in user-defined types. Use `sampleobj` to get this field.
* `tables`: `Dictionary{Symbol, T}`, a dictionary mapping data type to datatable. Use `tables` to get this field.

# Properties
All keys of `tables`.
"""
struct AnalysisTable{A, S, T <: AbstractDataTable{A, S}}
    analyte::Vector{A}
    sample::Vector{S}
    tables::Dictionary{Symbol, T}
    AnalysisTable(analyte::AbstractVector{A}, sample::AbstractVector{S}, tables::Dictionary{Symbol, T}) where {A, S, T <: AbstractDataTable{A, S}} = new{A, S, T}(convert(Vector{A}, analyte), convert(Vector{S}, sample), tables)
end

"""
    AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S}})

A `Dictionary`-like constructor for `AnalysisTable` from two iterable inputs `keys` and `tables`. The first value of `keys` will be the index for the first value of `tables`.
"""
function AnalysisTable(keys::AbstractVector{Symbol}, tables::AbstractVector{<: AbstractDataTable{A, S}}) where {A, S}
    allequal(analyteobj(table) for table in tables) || throw(ArgumentError("Tables should have identical analyte"))
    allequal(sampleobj(table) for table in tables) || throw(ArgumentError("Tables should have identical sample"))
    AnalysisTable(deepcopy(analyteobj(first(tables))), deepcopy(sampleobj(first(tables))), Dictionary(keys, tables))
end
"""
    analysistable(iter)

A `dictionary`-like function for constructing `AnalysisTable` from an iterable iter of key-value `Pair`s (or other iterables of two elements, such as a two-tuples). Keys should be `Symbol`s, and values should be `AbstractDataTable`s.
"""
function analysistable(iter)
    allequal(analyteobj(last(kv)) for kv in iter) || throw(ArgumentError("Tables should have identical analyte"))
    allequal(sampleobj(last(kv)) for kv in iter) || throw(ArgumentError("Tables should have identical sample"))
    v = last(first(iter))
    AnalysisTable(deepcopy(analyteobj(v)), deepcopy(sampleobj(v)), dictionary(iter))
end

"""
    AnalysisMethod{A, M <: Table, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}

A type containing analytes settings, and source calibration data. `A` determines analyte type.

# Fields
* `analytetable`: `M` contaning three columns.
    * `analyte`: `AbstractVector{A}`, analytes in user-defined types.
    * `isd`: `AbstractVector{Int}`, index of internal standard. `0` means no internal standard, and `-1` means the analyte itself is a internal standard.
    * `calibration`: index of analyte for calibration curve. `-1` means the analyte itself is a internal standard, so it will not be put into any calibration curve.
* `signal`: `Symbol`, type and key name of experimental acquisition data, e.g. `:area`.
* `rel_sig`: `Symbol`, key name of relative signal.
* `est_conc`: `Symbol`, key name of estimated concentration.
* `true_conc`: `Symbol`, key name of true concentration.
* `acc`: `Symbol`, key name of accuracy.
* `pointlevel`: `Vector{Int}` matching each point to level. It can be empty if there is only one level in `conctable`.
* `conctable`: `C <: AbstractDataTable{<: A, Int}` containing concentration data for each level. Samples must integers. One level indicates using `SingleCalibration`.
* `signaltable`: `D <: AbstractDataTable{A}` containig signal for each point. It can be `nothing` if signal data is unecessary.

# Properties
* `analyte`: `AbstractVector{A}`, analytes in user-defined types, identical to `analytetable.analyte`.
* `isd`: `AbstractVector{A}`, analytes which are internal standards.
* `nonisd`: `AbstractVector{A}`, analytes which are not internal standards.
* `point`: `AbstractVector{S}`, calibration points, identical to `sampleobj(signaltable)`. If `signaltable` is `nothing`, this value is `nothing` too.
* `level`: `AbstractVector{Int}`, calibration levels, identical to `sampleobj(conctable)`.
"""
struct AnalysisMethod{A, M <: Table, C <: AbstractDataTable, D <: Union{AbstractDataTable, Nothing}}
    analytetable::M
    signal::Symbol
    rel_sig::Symbol
    est_conc::Symbol
    true_conc::Symbol
    acc::Symbol
    pointlevel::Vector{Int}
    conctable::C
    signaltable::D
    function AnalysisMethod(
                    analytetable::M,
                    signal::Symbol,
                    rel_sig::Symbol,
                    est_conc::Symbol,
                    true_conc::Symbol,
                    acc::Symbol,
                    pointlevel::AbstractVector{Int}, 
                    conctable::AbstractDataTable{B, Int}, 
                    signaltable::AbstractDataTable{C, S}
                ) where {B, M <: Table, C, S}
        cp = propertynames(analytetable)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analytetable"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analytetable"))
        :calibration in cp || throw(ArgumentError("Column `:calibration` is required in `analytetable"))
        A = eltype(analytetable.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable $B should be a subtype of $A"))
        C <: A || throw(ArgumentError("Analyte type of signaltable $C should be a subtype of $A"))
        for a in analyteobj(conctable)
            a in analytetable.analyte || throw(ArgumentError("Analyte `$a` is not in the `analytetable`."))
        end
        if length(sampleobj(conctable)) > 1
            length(pointlevel) == length(sampleobj(signaltable)) || throw(ArgumentError("The length of `pointlevel` is different from that of `sampleobj(signaltable)`."))
            for a in analyteobj(conctable)
                a in analyteobj(signaltable) || throw(ArgumentError("Analyte `$a` is not in the `signatable`."))
            end
        end
        new{A, M, typeof(conctable), typeof(signaltable)}(analytetable, signal, rel_sig, est_conc, true_conc, acc, convert(Vector{Int}, pointlevel), conctable, signaltable)
    end
    function AnalysisMethod(
                    analytetable::M,
                    signal::Symbol,
                    rel_sig::Symbol,
                    est_conc::Symbol,
                    true_conc::Symbol,
                    acc::Symbol,
                    pointlevel::AbstractVector{Int}, 
                    conctable::AbstractDataTable{B, Int}, 
                    signaltable::Nothing
                ) where {B, M <: Table}
        cp = propertynames(analytetable)
        :analyte in cp || throw(ArgumentError("Column `:analyte` is required in `analytetable"))
        :isd in cp || throw(ArgumentError("Column `:isd` is required in `analytetable"))
        :calibration in cp || throw(ArgumentError("Column `:calibration` is required in `analytetable"))
        A = eltype(analytetable.analyte)
        B <: A || throw(ArgumentError("Analyte type of conctable should be a subtype of $A"))
        for a in analyteobj(conctable)
            a in analytetable.analyte || throw(ArgumentError("Analyte `$a` is not in the `analytetable`."))
        end
        new{A, M, typeof(conctable), Nothing}(analytetable, signal, rel_sig, est_conc, true_conc, acc, convert(Vector{Int}, pointlevel), conctable, signaltable)
    end
end

"""
    AnalysisMethod(
        conctable::AbstractDataTable, 
        signaltable::Union{AbstractDataTable, Nothing}, 
        signal, 
        pointlevel = Int[]; 
        signal = :area, 
        rel_sig = :relative_signal, 
        est_conc = :estimated_concentration, 
        true_conc = :true_concentration, 
        kwargs...
    )

Convenient contructors for `AnalysisMethod` which make constructing a `Table` optional.

`signal` data type for quantification. `levelname` is the column name for `pointlevel`. See `AnalysisMethod` for details of other arguments.

# Keyword arguments
* `rel_sig`: key name of relative signal. 
* `est_conc`: key name of estimated concentration. 
* `true_conc`: key name of true concentration. 
* `acc`: key name of accuracy. 
* Other keyword arguments will be columns in `analytetable`; when `analyte`, `isd` and `calibration` are not provided, it will use analyte in `conctable`. 
"""
function AnalysisMethod(
        conctable::AbstractDataTable, 
        signaltable::Union{AbstractDataTable, Nothing},
        signal,
        pointlevel = nothing; 
        rel_sig = :relative_signal,
        est_conc = :estimated_concentration,
        true_conc = :true_concentration, 
        acc = :accuracy, 
        kwargs...
    )
    signal = Symbol(signal)
    rel_sig = Symbol(rel_sig)
    est_conc = Symbol(est_conc)
    true_conc = Symbol(true_conc)
    acc = Symbol(acc)
    if !isa(pointlevel, AbstractVector{<: Integer})
        pointlevel = getproperty(signaltable, Symbol(pointlevel))
    end
    analyte = get(kwargs, :analyte, analyteobj(conctable))
    isd = get(kwargs, :isd, zeros(Int, length(analyteobj(conctable))))
    calibration = get(kwargs, :calibration, collect(eachindex(analyteobj(conctable))))
    analytetable = Table(; analyte, isd, calibration, kwargs...)    
    AnalysisMethod(analytetable, signal, rel_sig, est_conc, true_conc, acc, pointlevel, conctable, signaltable)
end

"""
    MultipleCalibration{A, N, T <: Table} <: AbstractCalibration{A, N}

A mutable type holding all data and settings for a calibration curve. `A` determines analyte type, and `N` determines numeric type. Only `Float64` is supported because of issue of `GLM`.

# Fields
* `analyte`: `Tuple{A, Any}`. First element is the analyte being quantified, and the second element is its internal standard for which `nothing` indicates no internal standard.
* `type`: `Bool` determines whether fitting a linear line (`true`) or quadratic curve (`false`).
* `zero`: `Bool` determines whether forcing the curve crossing (0, 0) (`true`) or ignoring it (`false`).
* `weight`: `Float64` represents the exponential applying to each element of `x` as a weighting vector.
* `formula`: `FormulaTerm`, the formula for fitting calibration curve.
* `table`: `TypedTable.Table`, the cleaned up calibration data, containing 7 columns.
    * `id`: Point name
    * `level`: The index of concentration level. The larger, the higher concentraion it represents.
    * `y`: Signal or relative signal
    * `x`: True concentraion
    * `x̂`: Predicted concentration
    * `accuracy`: Accuracy, i.e. `x̂/x`.
    * `include`: Whether this point is included or not
* `model`: `GLM` object.
"""
mutable struct MultipleCalibration{A, N, T <: Table} <: AbstractCalibration{A, N}
    analyte::Tuple{A, Any}
    type::Bool
    zero::Bool
    weight::N
    formula::FormulaTerm
    table::T
    model
    function MultipleCalibration(analyte::Tuple{A, B}, type, zero, weight, formula, table::T, model) where {A, T, B}
        N = eltype(table.x)
        N == eltype(table.y) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column y is `", eltype(table.y), "`.")))
        N == eltype(table.x̂) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column x̂ is `", eltype(table.x̂), "`.")))
        N == eltype(table.accuracy) || throw(ArgumentError(string("Element of column x is `", N, "`, while element of column accuracy is `", eltype(table.accuracy), "`.")))
        new{A, N, T}(analyte, type, zero, N(weight), formula, table, model)
    end
end

"""
    SingleCalibration{A, N} <: AbstractCalibration{A, N}

A mutable type holding all data for single point calibration.

# Fields
* `analyte`: `Tuple{A}` is the analyte with known concentration (internal standard).
* `conc`: `Float64`, concentration of analyte.
"""
mutable struct SingleCalibration{A, N} <: AbstractCalibration{A, N}
    analyte::Tuple{A}
    conc::N
end

"""
    Batch{A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}, D <: Union{AnalysisTable{<: A}, Nothing}}

A type representing a batch for quantitative analysis. `A` determines analyte type, and `T` determines table type.

# Fields
* `method`: `M`, method.
* `calibration`: `C`, calibration curves.
* `data`: `D`, data for analysis, .

# Properties
* `analyte`: `AbstractVector{A}`, analytes in user-defined types, identical to `method.analyte`.
* `isd`: `AbstractVector{<: A}`, analytes which are internal standards, identical to `method.isd`.
* `nonisd`: `AbstractVector{<: A}`, analytes which are not internal standards, identical to `method.nonisd`.
* `point`: `AbstractVector` or `Nothing`, calibration points, identical to `method.point`.
* `level`: `AbstractVector{Int}`, calibration levels, identical to `method.level`.

# Constructors
* `Batch(method, calibration, data = nothing)`
"""
struct Batch{A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}, D <: Union{AnalysisTable{<: A}, Nothing}}
    method::M
    calibration::C
    data::D
    Batch(method::M, calibration::C, data::D) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}, D <: AnalysisTable{<: A}} = new{A, M, C, D}(method, calibration, data)
    Batch(method::M, calibration::C, data::Nothing) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}} = new{A, M, C, Nothing}(method, calibration, data)
    Batch(method::M, calibration::C) where {A, M <: AnalysisMethod{A}, C <: AbstractVector{<: AbstractCalibration{<: A}}} = new{A, M, C, Nothing}(method, calibration, nothing)
end

"""
    Batch(method::AnalysisMethod, data = nothing; type = true, zero = false, weight = 0)

Construct a `Batch` from `method`, and optionally `data` with specified calibration parameters. See "MultipleCalibration" for detail description of keyword arguments.
"""
function Batch(method::AnalysisMethod{A}, data = nothing;
                type = true, 
                zero = false, 
                weight = 0
            ) where A
    Batch(
        method,
        if length(sampleobj(method.conctable)) > 1 
            map(analyteobj(method.conctable)) do analyte
                calibration(method, analyte; type, zero, weight)
            end
        else
            map(analyteobj(method.conctable)) do analyte
                SingleCalibration((analyte, ), first(getanalyte(method.conctable, analyte)))
            end
        end
        ,
        data
    )
end
"""
    Batch(batch::Batch, at::AnalysisTable)

Construct a new batch from an old batch and a `AnalysisTable` as new data.
"""
Batch(batch::Batch, at::AnalysisTable) = Batch(batch.method, batch.calibration, at)
"""
    Batch(dt::AbstractDataTable; 
        signal = :area, 
        rel_sig = :relative_signal, 
        est_conc = :estimated_concentration, 
        true_conc = :true_concentration, 
        acc = :accuracy, 
        calid = r"Cal_(\\d)_(\\d*-*\\d*)", 
        order = "LR", 
        f2c = 1,
        parse_decimal = x -> replace(x, "-" => "."))

Construct a batch from data. All analytes are considered as normal analytes, so calibration curves are not contructed immediately; call `init_calibration!` after confirming isd and calibration settings in `batch.method.analytetable`.

# Arguments
* `dt`: data containing both calibration curves and samples.

# Keyword Arguments
* `signal`: type and key name of experimental acquisition data, e.g. `:area`.
* `rel_sig`: key name of relative signal.
* `est_conc`: key name of estimated concentration.
* `true_conc`: key name of true concentration.
* `acc`: key name of accuracy.
* `calid`: `Regex`, identifier for calibration data; level and concentration related factors (ratio of concentration or dilution factor) should be captured. 
The former should be able to be parsed as integer directly; the latter should be able to be parsed as floating number after applying `parse_decimal`.
* `order`: `String`, represents the order and identity of captured string; `L` is level, `R` is ratio of concentration, `D` is dilution factor (df).
* `f2c`: `Number` or vector of numbers; concentration equals to f2c * ratio or f2c / df. When a vector is provided, each element represents `f2c` value of each analyte.
* `parse_decimal`: `Function`, converts a string into another string which can be parsed as floating number.
"""
function Batch(dt::AbstractDataTable; 
                signal = :area, 
                rel_sig = :relative_signal, 
                est_conc = :estimated_concentration, 
                true_conc = :true_concentration, 
                acc = :accuracy, 
                calid = r"Cal_(\d)_(\d*-*\d*)", 
                order = "LR", 
                f2c = 1,
                parse_decimal = x -> replace(x, "-" => "."))
    signal = Symbol(signal)
    rel_sig = Symbol(rel_sig)
    est_conc = Symbol(est_conc)
    true_conc = Symbol(true_conc)
    acc = Symbol(acc)
    dilution = occursin("D", order)
    levelid = findfirst(==("L"), split(order, ""))
    tbl = Table(table(dt))
    calname = map(samplename(dt)) do s
        parse_calibration_name(s; calid, levelid, dilution, f2c, parse_decimal)
    end
    aj = analyteobj(dt)
    id = idcol(dt)
    sj = sampleobj(dt)
    if dt isa SampleDataTable
        sampledata = SampleDataTable(aj, sj[isnothing.(calname)], id, tbl[isnothing.(calname)])
        signaltable = SampleDataTable(aj, sj[(!isnothing).(calname)], id, tbl[(!isnothing).(calname)])
    else
        sampledata = AnalyteDataTable(aj, sj[isnothing.(calname)], id, getproperties(tbl, tuple(samplename(dt)[isnothing.(calname)]...)))
        signaltable = AnalyteDataTable(aj, sj[(!isnothing).(calname)], id, getproperties(tbl, tuple(samplename(dt)[(!isnothing).(calname)]...)))
    end
    filter!(!isnothing, calname)
    pointlevel = map(first, calname)
    unique!(calname)
    levels = map(first, calname)
    concs = map(last, calname)
    if dt isa SampleDataTable
        conctable = SampleDataTable(Table(; Level = levels, (analytename(dt) .=> (collect(i) for i in zip(concs...)))...), :Level)
    else
        conctable = AnalyteDataTable(Table(; (id => aj, )..., (Symbol.(levels) .=> map(x -> repeat([0], length(aj)) .+ x, concs))...), id, Int)
    end
    analytetable = Table(; analyte = aj, isd = zeros(Int, length(aj)), calibration = collect(1:length(aj)))
    Batch(AnalysisMethod(analytetable, signal, rel_sig, est_conc, true_conc, acc, pointlevel, conctable, signaltable), MultipleCalibration{eltype(analyteobj(dt))}[], analysistable((signal => sampledata, )))
end

include("interface.jl")
include("utils.jl")
include("cal.jl")
include("quant.jl")
include("io.jl")

"""
    SampleDataTable(tbl::AnalyteDataTable{A, S, N, T}, samplecol::Symbol, tablesink::TypeOrFn = T)
    SampleDataTable(samplecol::Symbol, tablesink::TypeOrFn, tbl::AnalyteDataTable)
    SampleDataTable(samplecol::Symbol, tbl::AnalyteDataTable)

Convert `AnalyteDataTable` to `SampleDataTable` with `samplecol` as the column name of sample.
"""
SampleDataTable(tbl::AnalyteDataTable{A, S, N, T}, samplecol::Symbol, tablesink::TypeOrFn = T) where {A, S, N, T} = 
    SampleDataTable(analyteobj(tbl), samplecol, tablesink((; (samplecol => sampleobj(tbl), (analytename(tbl) .=> eachanalyte(tbl))...)...)))
SampleDataTable(samplecol::Symbol, tablesink::TypeOrFn, tbl::AnalyteDataTable) = 
    SampleDataTable(tbl, samplecol, tablesink)
SampleDataTable(samplecol::Symbol, tbl::AnalyteDataTable) = 
    SampleDataTable(tbl, samplecol)

"""
    AnalyteDataTable(tbl::SampleDataTable{A, S, N, T}, analytecol::Symbol, tablesink::TypeOrFn = T)
    AnalyteDataTable(analytecol::Symbol, tablesink::TypeOrFn, tbl::SampleDataTable)
    AnalyteDataTable(analytecol::Symbol, tbl::SampleDataTable)

Convert `SampleDataTable` to `AnalyteDataTable` with `analytecol` as the column name of analyte.
"""
AnalyteDataTable(tbl::SampleDataTable{A, S, N, T}, analytecol::Symbol, tablesink::TypeOrFn = T) where {A, S, N, T} = 
    AnalyteDataTable(analytecol, sampleobj(tbl), tablesink((; (analytecol => analyteobj(tbl), (samplename(tbl) .=> eachsample(tbl))...)...)))
AnalyteDataTable(analytecol::Symbol, tablesink::TypeOrFn, tbl::SampleDataTable) = 
    AnalyteDataTable(tbl, analytecol, tablesink)
AnalyteDataTable(analytecol::Symbol, tbl::SampleDataTable) = 
    AnalyteDataTable(tbl, analytecol)

end
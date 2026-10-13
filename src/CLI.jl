module CLI

using ArgCheck
using ArgParse
using Dates
using FASTX
using Logging
using Random
using Statistics
using TOML
using ..PottsEvolver
using ..PottsEvolver: TreeTools

#================================#
######### Parameter file #########
#================================#

"""
    load_params(path) -> SamplingParameters

Read a TOML file at `path` into a `SamplingParameters`.
Any field of `SamplingParameters` may be given as a top-level key; fields that are omitted
fall back to `SamplingParameters`'s own defaults. `sampling_type`, `step_type` and
`step_meaning` are given as plain strings (*e.g.* `sampling_type = "continuous"`) and are
converted to `Symbol`s here.

`branchlength_meaning`, if given, is a subtable:
```toml
[branchlength_meaning]
type = "step"
length = "exact"
```

The `[run]` and `[schedule]` tables, as written by [`write_params`](@ref), are ignored here:
they hold metadata about a run and the sampling schedule, neither of which is a sampling
parameter, so a `<prefix>.params.toml` written by the CLI can be given back to it unchanged.
Read the schedule with [`load_schedule`](@ref).

`mutation_matrix`, if given, is a nested array **of rows**: `mutation_matrix[a]` is the
`a`-th row of the matrix, *i.e.* the relative mutation rates *out of* state `a` towards
every other state `b` (`mutation_matrix[a][b]` becomes `μ[a, b]`).
```toml
mutation_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
```
"""
function load_params(path::AbstractString)
    # Unknown keys surface as `SamplingParameters`'s own `ArgumentError`.
    raw = TOML.parsefile(path)
    kwargs = Dict{Symbol,Any}()

    symbol_fields = ("sampling_type", "step_type", "step_meaning")
    for (k, v) in raw
        if k in ("run", "schedule")
            # run metadata and sampling schedule, not sampling parameters: see the docstring
            continue
        elseif k in ("Teq", "burnin")
            throw(
                ArgumentError(
                    """
                    `$k` is no longer a sampling parameter: the schedule moved to the \
                    `[schedule]` table (and to -M/--spacing/--burnin). Replace

                        $k = $v

                    with, for instance

                        [schedule]
                        spacing = <time between samples>
                        burnin = <time before the first sample>
                        n_samples = <number of samples>
                    """
                ),
            )
        elseif k == "branchlength_meaning"
            kwargs[Symbol(k)] = BranchLengthMeaning(Symbol(v["type"]), Symbol(v["length"]))
        elseif k == "mutation_matrix"
            kwargs[Symbol(k)] = permutedims(reduce(hcat, v))
        elseif k in symbol_fields
            kwargs[Symbol(k)] = Symbol(v)
        else
            kwargs[Symbol(k)] = v
        end
    end

    return SamplingParameters(; kwargs...)
end

#================================#
######### Init sequence #########
#================================#

const _INIT_KEYWORDS = ("random_aa", "random_codon", "random_rna", "random_num")

# Symbol sets of the three fasta-representable alphabets. `CodonSequence`'s nucleotides are
# written with `T` (not `U`, see `PottsEvolver.nucleotides`), so `T` is shared with the amino
# acid alphabet (Threonine) and `U` is exclusive to `RNASequence` - the two are genuinely
# different character sets, not "DNA vs RNA" spelling of the same thing.
const _CODON_NT_SYMBOLS = Set(vcat(PottsEvolver.nucleotides, ['-']))
const _RNA_SYMBOLS = Set(PottsEvolver.symbols(PottsEvolver.rna_alphabet))
const _AA_SYMBOLS = Set(PottsEvolver.symbols(PottsEvolver.aa_alphabet))
const _AA_EXCLUSIVE_SYMBOLS = setdiff(_AA_SYMBOLS, _CODON_NT_SYMBOLS, _RNA_SYMBOLS)
const _RNA_EXCLUSIVE_SYMBOLS = setdiff(_RNA_SYMBOLS, _CODON_NT_SYMBOLS, _AA_SYMBOLS)

"""
    load_init(init, init_type, g::PottsGraph; rng=Random.default_rng()) -> AbstractSequence

Resolve the `--init` CLI value into an initial/root sequence.

`init` is either:
- one of `"random_aa"`, `"random_codon"`, `"random_rna"`, `"random_num"`: forwarded to
  `get_init_sequence`;
- a path to a fasta file, optionally followed by a sequence name, mirroring IQ-TREE's
  alisim `--root-seq` syntax: `"file.fasta"` (must contain exactly one record) or
  `"file.fasta,seqname"` (look up the record whose identifier matches `seqname` exactly).

For the fasta forms, the sequence type (`AASequence`/`CodonSequence`/`RNASequence`) is
inferred from the sequence's characters and, when that alone is ambiguous, from its length
relative to `g`'s `L`/`q`. Pass `init_type` (one of `"aa"`, `"codon"`, `"rna"`) to skip
inference and force a type; it is otherwise only needed when inference genuinely fails.
`NumSequence` cannot be built from a fasta file (it has no symbol table) - use the
`"random_num"` keyword instead.
"""
function load_init(
    init::AbstractString,
    init_type::Union{Nothing,AbstractString},
    g::PottsGraph;
    rng=Random.default_rng(),
)
    if init in _INIT_KEYWORDS
        return PottsEvolver.get_init_sequence(Symbol(init), g; rng)
    end

    file, seqname = _split_init_path(init) # split at first comma
    @argcheck isfile(file) "No such file: `$file`"
    s = _read_fasta_sequence(file, seqname)

    T = if isnothing(init_type)
        _infer_sequence_type(s, g)
    else
        _sequence_type_from_name(init_type)
    end
    return T(s)
end

function _split_init_path(init::AbstractString)
    parts = split(init, ","; limit=2)
    return if length(parts) == 2
        (String(parts[1]), String(parts[2]))
    else
        (String(parts[1]), nothing)
    end
end

function _read_fasta_sequence(file::AbstractString, seqname::Nothing)
    records = collect(FASTAReader(open(file)))
    @argcheck length(records) == 1 """
        Expected exactly one sequence in `$file` (got $(length(records))). \
        Use `$file,<seqname>` to select one explicitly.
        """
    return FASTX.sequence(only(records))
end
function _read_fasta_sequence(file::AbstractString, seqname::AbstractString)
    records = collect(FASTAReader(open(file)))
    idx = findfirst(rec -> FASTX.identifier(rec) == seqname, records)
    @argcheck !isnothing(idx) """
        No sequence named `$seqname` in `$file`. \
        Available: $(join(FASTX.identifier.(records), ", "))
        """
    return FASTX.sequence(records[idx])
end

function _sequence_type_from_name(init_type::AbstractString)
    return if init_type == "aa"
        AASequence
    elseif init_type == "codon"
        CodonSequence
    elseif init_type == "rna"
        RNASequence
    else
        throw(
            ArgumentError(
                "Unknown `--init-type` `$init_type`. Expected `aa`, `codon` or `rna`."
            ),
        )
    end
end

function _infer_sequence_type(s::AbstractString, g::PottsGraph)
    (; q, L) = size(g)
    chars = Set(s)

    if any(in(_AA_EXCLUSIVE_SYMBOLS), chars)
        @argcheck q == PottsEvolver.Q_AA && length(s) == L """
            Sequence contains symbol(s) exclusive to amino acids \
            ($(join(intersect(chars, _AA_EXCLUSIVE_SYMBOLS), ", "))), but the graph is \
            incompatible with `AASequence` (q=$q, L=$L, sequence length=$(length(s))). \
            Pass `--init-type` explicitly if this is intentional.
            """
        return AASequence
    elseif any(in(_RNA_EXCLUSIVE_SYMBOLS), chars)
        @argcheck q == PottsEvolver.Q_RNA && length(s) == L """
            Sequence contains symbol(s) exclusive to RNA \
            ($(join(intersect(chars, _RNA_EXCLUSIVE_SYMBOLS), ", "))), but the graph is \
            incompatible with `RNASequence` (q=$q, L=$L, sequence length=$(length(s))). \
            Pass `--init-type` explicitly if this is intentional.
            """
        return RNASequence
    end

    # Content alone doesn't disambiguate (only `-`, `A`, `C`, `G`, and possibly `T`): fall
    # back to matching length and `q` against each candidate compatible with the characters.
    has_t = 'T' in chars
    matches = Type[]
    q == PottsEvolver.Q_AA && length(s) == L && push!(matches, AASequence)
    q == PottsEvolver.Q_AA && length(s) == 3L && push!(matches, CodonSequence)
    !has_t && q == PottsEvolver.Q_RNA && length(s) == L && push!(matches, RNASequence)

    if length(matches) == 1
        return only(matches)
    elseif isempty(matches)
        throw(
            ArgumentError("""
          Cannot infer sequence type: content is ambiguous (only `-ACGT`) and its length \
          ($(length(s))) does not uniquely match the graph (q=$q, L=$L). Pass \
          `--init-type` explicitly (`aa`, `codon` or `rna`).
          """)
        )
    else
        throw(
            ArgumentError(
                """
Cannot infer sequence type: length $(length(s)) is compatible with more than one \
candidate for this graph (q=$q, L=$L). Pass `--init-type` explicitly (`aa`, \
`codon` or `rna`).
"""
            ),
        )
    end
end

#================================#
############ Output ##############
#================================#

# Output files are all named `<prefix>.<suffix>` in `--outdir`, in the style of IQ-TREE.

"""
    default_prefix(prefix, model)

Resolve the `--prefix` value: if not given, use `model`'s basename without its extension.
"""
default_prefix(prefix::AbstractString, model) = prefix
default_prefix(::Nothing, model) = first(splitext(basename(model)))

"""
    output_file(outdir, prefix, suffix)

Path to `outdir/prefix.suffix`, creating `outdir` if needed.
Warn if the file exists: it is going to be overwritten.
"""
function output_file(outdir, prefix, suffix)
    mkpath(outdir)
    path = joinpath(outdir, "$prefix.$suffix")
    isfile(path) && @warn "Overwriting existing output file $path"
    return path
end

"""
    setup_logfile(outdir, prefix, command, args) -> (path, io)

Create the run's log file, write a header to it, and return its path along with the **open**
stream, which the caller owns and has to close at the end of the run.

The stream is meant to be wrapped with `PottsEvolver.file_sink` and passed to `mcmc_sample` as
`extra_sinks`, so that the CLI and the sampling functions log to the same file through a
*single* handle. Two handles would not do: Julia's append mode seeks to the end of the file
when opening it and each handle then keeps its own position, so their writes would land on top
of each other.
"""
function setup_logfile(outdir, prefix, command, args)
    path = output_file(outdir, prefix, "log")
    io = open(path, "w")
    println(io, "# pottsevolver $(pkgversion(PottsEvolver)) - $command - $(Dates.now())")
    println(io, "# args: $(join(args, ' '))")
    flush(io)
    return path, io
end

#=========== Parameters ===========#

_params_pairs(p::SamplingParameters) = (f => getproperty(p, f) for f in propertynames(p))
_params_pairs(d::AbstractDict) = (Symbol(k) => v for (k, v) in d)

function _params_toml_dict(params; run=nothing, schedule=nothing)
    out = Dict{String,Any}()
    run_table = if isnothing(run)
        Dict{String,Any}()
    else
        Dict{String,Any}(string(k) => v for (k, v) in run)
    end
    for (field, value) in _params_pairs(params)
        # TOML has no null; an omitted key falls back to `SamplingParameters`'s default, which
        # is exactly what an unset field means
        isnothing(value) && continue
        if field == :sequence_type
            # not a `SamplingParameters` field: belongs to the metadata table
            run_table["sequence_type"] = string(value)
        elseif field == :branchlength_meaning
            out[string(field)] = Dict(
                "type" => string(value.type), "length" => string(value.length)
            )
        elseif field == :mutation_matrix
            out[string(field)] = map(collect, eachrow(value)) # `load_params` reads rows
        else
            out[string(field)] = value
        end
    end
    isempty(run_table) || (out["run"] = run_table)
    isnothing(schedule) || (out["schedule"] = schedule)
    return out
end

"""
    write_params(path, params; run=nothing)

Write the parameters of a run to `path` as TOML.
`params` is either a `SamplingParameters` or the `params` field of an `mcmc_sample` output.

The file can be read back by [`load_params`](@ref): fields that are `nothing` are omitted
(`load_params` falls back to the same defaults), and anything that is not a
`SamplingParameters` field goes to a `[run]` table that `load_params` ignores. `run` holds
extra metadata to add there, *e.g.* the input files and the rng seed.

`schedule` is the `[schedule]` table, as produced by [`write_tvals`](@ref): the sampling times
are not sampling parameters, but they belong in the file so that a run can be reproduced from
it with `-p` alone.
"""
function write_params(path, params; run=nothing, schedule=nothing)
    open(path, "w") do io
        # `Symbol` is not a valid TOML type: `string` converts it, and `sequence_type` too
        return TOML.print(string, io, _params_toml_dict(params; run, schedule); sorted=true)
    end
    return path
end

#=========== Substitutions ===========#

"""
    write_substitutions(path, info, tvals)

Write the substitutions tracked during a continuous chain run to `path`, as a csv file with
columns `sample,time,position,old,new`.

`sample` is the time value of the sample ending the interval the substitution happened in, and
`time` the absolute time of the substitution, counted from the start of the run (*i.e.*
including `burnin`). `position`, `old` and `new` index states of the sampled sequence type,
which is recorded in the `[run]` table of the parameter file.
"""
function write_substitutions(path, info, tvals)
    open(path, "w") do io
        println(io, "sample,time,position,old,new")
        # `info[m]` covers the interval ending at `tvals[m+1]`: the first sample is the initial
        # sequence, with nothing sampled before it
        for (m, entry) in enumerate(info), (mutation, t) in entry.substitutions
            println(
                io, join((tvals[m + 1], t, mutation.pos, mutation.old, mutation.new), ',')
            )
        end
    end
    return path
end

_tracks_substitutions(params) = get(params, :track_substitutions, false)

function log_info_summary(info, sampling_type)
    isempty(info) && return nothing
    if sampling_type == :discrete
        @info "Average fraction of accepted steps: $(mean(x -> x.ratio, info))"
    elseif sampling_type == :continuous
        @info "Total number of substitutions: $(sum(x -> x.number_substitutions, info))"
    end
    return nothing
end

#=========== Writing a run ===========#

"""
    write_chain_output(outdir, prefix, result; run=nothing)

Write the output of a chain run (`mcmc_sample(g, M, params)`) to `outdir`, and return the paths
written:
- `<prefix>.fasta`: the sample, labelled by sampling time;
- `<prefix>.params.toml`: parameters of the run, see [`write_params`](@ref);
- `<prefix>.substitutions.csv`: only if substitutions were tracked, see
  [`write_substitutions`](@ref).
"""
function write_chain_output(outdir, prefix, result; run=nothing, schedule=nothing)
    files = String[]
    schedule_table, sidecar = if isnothing(schedule)
        nothing, nothing
    else
        write_tvals(schedule; prefix_path=joinpath(outdir, prefix))
    end

    path = output_file(outdir, prefix, "fasta")
    write_fasta(path, result.sequences)
    push!(files, path)
    @info "Wrote $(length(result.sequences)) sequences to $path"

    push!(
        files,
        write_params(
            output_file(outdir, prefix, "params.toml"),
            result.params;
            run,
            schedule=schedule_table,
        ),
    )

    if !isnothing(sidecar)
        (dir, name), times = sidecar
        path = output_file(dir, first(splitext(name)), "txt")
        open(path, "w") do io
            return foreach(t -> println(io, t), times)
        end
        push!(files, path)
        @info "Wrote the sampling times to $path"
    end

    if _tracks_substitutions(result.params)
        path = write_substitutions(
            output_file(outdir, prefix, "substitutions.csv"), result.info, result.tvals
        )
        push!(files, path)
        @info "Wrote tracked substitutions to $path"
    end

    log_info_summary(result.info, get(result.params, :sampling_type, nothing))
    return files
end

"""
    write_tree_output(outdir, prefix, result; run=nothing, internals=false)

Write the output of a tree run (`mcmc_sample(g, tree, params)`) to `outdir`, and return the
paths written:
- `<prefix>.fasta`: the leaf sequences, labelled by node;
- `<prefix>.params.toml`: parameters of the run, see [`write_params`](@ref);
- `<prefix>.internals.fasta` and `<prefix>.nwk`: only if `internals`.

The tree is written along with the internal sequences because `read_tree` labels the internal
nodes that the input tree leaves unnamed: the output tree is the one whose labels the sequences
refer to.
"""
function write_tree_output(outdir, prefix, result; run=nothing, internals=false)
    files = String[]

    path = output_file(outdir, prefix, "fasta")
    write_fasta(path, result.leaf_sequences)
    push!(files, path)
    @info "Wrote $(length(result.leaf_sequences)) leaf sequences to $path"

    if internals
        path = output_file(outdir, prefix, "internals.fasta")
        write_fasta(path, result.internal_sequences)
        push!(files, path)
        @info "Wrote $(length(result.internal_sequences)) internal sequences to $path"

        path = output_file(outdir, prefix, "nwk")
        TreeTools.write_newick(path, result.tree)
        push!(files, path)
        @info "Wrote the sampled tree to $path"
    end

    push!(
        files, write_params(output_file(outdir, prefix, "params.toml"), result.params; run)
    )

    if _tracks_substitutions(result.params)
        @warn """
        `track_substitutions` is set, but tree sampling does not return substitutions:
        no substitution file is written.
        """
    end
    return files
end

#================================#
######### Running a command ######
#================================#

"""
    load_tvals(path)

Read an explicit sampling schedule from `path`: one number per line.
Blank lines, and lines whose first non-space character is `#`, are skipped.
"""
function load_tvals(path::AbstractString)
    tvals = Float64[]
    for (n, line) in enumerate(eachline(path))
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#')) && continue
        value = tryparse(Float64, stripped)
        isnothing(value) &&
            throw(ArgumentError("$path, line $n: expected a number, got `$stripped`"))
        push!(tvals, value)
    end
    return tvals
end

#=========== Sampling schedule ===========#

# A schedule is carried around as a small spec rather than a vector of times, so that it can
# be written back to a parameter file without storing every value: a 100_000-sample log
# schedule is four numbers. `kind` says how to expand it.
const _SCHEDULE_KINDS = ("linear", "log")

"""
    schedule_times(spec) -> Vector

Expand a schedule spec into the vector of sampling times passed to `mcmc_sample`.

- `(; kind=:linear, burnin, spacing, n_samples)`: `burnin`, then `n_samples` samples `spacing`
  apart;
- `(; kind=:linear, start, stop, n_samples)`: `n_samples` times from `start` to `stop`, evenly
  spaced, *i.e.* `range`;
- `(; kind=:log, start, stop, n_samples)`: the same between `start` and `stop`, but spaced
  geometrically, *i.e.* `logrange`. `start` must be `> 0`;
- `(; kind=:explicit, tvals)`: the times themselves.

The two `linear` shapes are the same schedule said two ways: `spacing` is the step between
samples, `stop` the time of the last one.
"""
function schedule_times(spec)
    spec.kind == :explicit && return collect(spec.tvals)
    @argcheck spec.kind in (:linear, :log) "Unknown schedule kind `$(spec.kind)`"

    # `range`/`logrange` reject a single point unless start == stop, so take it directly
    if spec.kind == :linear && hasproperty(spec, :spacing)
        return collect(
            spec.burnin .+
            range(0, (spec.n_samples - 1) * spec.spacing; length=spec.n_samples),
        )
    end
    spec.n_samples == 1 && return [spec.start]
    return if spec.kind == :linear
        collect(range(spec.start, spec.stop, spec.n_samples))
    else
        @argcheck spec.start > 0 """
        A log-spaced schedule needs `start > 0`. Instead $(spec.start).
        """
        collect(logrange(spec.start, spec.stop, spec.n_samples))
    end
end

"""
    parse_tvals(table; dir=".") -> spec

Build a schedule spec from the `[schedule]` table of a parameter file. See
[`schedule_times`](@ref) for the shapes, and [`write_tvals`](@ref) for the inverse.

`dir` is where a `tvals_file` entry is resolved from, so that a parameter file written next to
its sidecar keeps working when the pair is moved together.
"""
function parse_tvals(table; dir=".")
    haskey(table, "tvals") && return (; kind=:explicit, tvals=collect(table["tvals"]))
    if haskey(table, "tvals_file")
        return (; kind=:explicit, tvals=load_tvals(joinpath(dir, table["tvals_file"])))
    end

    kind = get(table, "kind", "linear")
    @argcheck kind in _SCHEDULE_KINDS """
    `[schedule] kind` should be one of $(join(_SCHEDULE_KINDS, ", ")). Instead `$kind`.
    """
    # `start`/`stop` shapes must be whole: no flag completes them, unlike `spacing`/`burnin`
    if kind == "log" || haskey(table, "stop")
        return (;
            kind=Symbol(kind),
            start=_required(table, "start"),
            stop=_required(table, "stop"),
            n_samples=_required(table, "n_samples"),
        )
    end
    return (;
        kind=:linear,
        burnin=get(table, "burnin", nothing),
        spacing=get(table, "spacing", nothing),
        n_samples=get(table, "n_samples", nothing),
    )
end

"""
    normalize_schedule(spec) -> spec

Fill in the entries a `linear` schedule leaves unset: `n_samples` defaults to 1, `spacing` to
1 and `burnin` to `5 * spacing`, matching `mcmc_sample`'s own defaults.
"""
function normalize_schedule(spec)
    (spec.kind == :linear && hasproperty(spec, :spacing)) || return spec
    spacing = something(spec.spacing, 1)
    return (;
        kind=:linear,
        burnin=something(spec.burnin, 5 * spacing),
        spacing,
        n_samples=something(spec.n_samples, 1),
    )
end

function _required(table, key)
    haskey(table, key) || throw(
        ArgumentError(
            "`[schedule]` is missing the `$key` entry (kind=$(get(table, "kind", "linear")))",
        ),
    )
    return table[key]
end

"""
    write_tvals(spec; prefix_path=nothing) -> (table, sidecar)

Render a schedule spec as the `[schedule]` table of a parameter file, in the same shape it was
given: a structured schedule stays four numbers, and an explicit one keeps its values.

`sidecar` is `nothing`, or a `path => times` pair to be written alongside: an explicit schedule
that came from a file is stored as a file again rather than inlined, so that a parameter file
stays readable no matter how many samples were asked for. It needs `prefix_path`, the output
path stem.
"""
function write_tvals(spec; prefix_path=nothing)
    if spec.kind == :linear && hasproperty(spec, :spacing)
        return Dict{String,Any}(
            "kind" => "linear",
            "burnin" => spec.burnin,
            "spacing" => spec.spacing,
            "n_samples" => spec.n_samples,
        ),
        nothing
    elseif spec.kind in (:linear, :log)
        return Dict{String,Any}(
            "kind" => string(spec.kind),
            "start" => spec.start,
            "stop" => spec.stop,
            "n_samples" => spec.n_samples,
        ),
        nothing
    end

    # explicit: inline a schedule that was typed out, keep a file one as a file
    inline = !get(spec, :from_file, false) || isnothing(prefix_path)
    inline && return Dict{String,Any}("tvals" => collect(spec.tvals)), nothing
    sidecar = "$(basename(prefix_path)).tvals.txt"
    return Dict{String,Any}("tvals_file" => sidecar),
    (dirname(prefix_path), sidecar) => spec.tvals
end

"""
    resolve_schedule(args, params_path) -> spec

Resolve the sampling schedule for `sample-chain`, from the command line and the `[schedule]`
table of the parameter file.

Flags win over the file. `--tvals`/`--tvals-file` give the times directly; `-M`, `--spacing`
and `--burnin` describe an equispaced schedule and are filled in from a `linear` table in the
file for whatever they leave unset. Asking for times directly *and* describing them is a
conflict, not a precedence: the file's table is then ignored with a warning, since the two
cannot be combined meaningfully.
"""
function resolve_schedule(args, params_path)
    from_file = load_schedule(params_path)

    explicit_flags = filter(
        !isnothing,
        [
            isempty(args["tvals"]) ? nothing : "--tvals",
            isnothing(args["tvals-file"]) ? nothing : "--tvals-file",
        ],
    )
    @argcheck length(explicit_flags) <= 1 """
    --tvals and --tvals-file both give the sampling times: keep only one.
    """
    linear_flags = filter(
        !isnothing,
        [
            isnothing(args["n-samples"]) ? nothing : "-M",
            isnothing(args["spacing"]) ? nothing : "--spacing",
            isnothing(args["burnin"]) ? nothing : "--burnin",
        ],
    )
    @argcheck isempty(explicit_flags) || isempty(linear_flags) """
    $(join(vcat(explicit_flags, linear_flags), ", ")) cannot be combined: the first gives the \
    sampling times, the others describe them. Keep one or the other.
    """

    # times given directly on the command line
    if !isempty(explicit_flags)
        _warn_ignored_table(from_file, join(explicit_flags, ", "))
        return if isempty(args["tvals"])
            (; kind=:explicit, tvals=load_tvals(args["tvals-file"]), from_file=true)
        else
            (; kind=:explicit, tvals=args["tvals"])
        end
    end

    # an equispaced schedule, with the file supplying anything the flags leave out
    if !isempty(linear_flags)
        base = if !isnothing(from_file) && from_file.kind == :linear
            # a start/stop table is restated as burnin/spacing so a flag can override one entry
            as_spacing_form(from_file)
        else
            _warn_ignored_table(from_file, join(linear_flags, ", "))
            (; kind=:linear, burnin=nothing, spacing=nothing, n_samples=nothing)
        end
        # `something` would throw when both are unset, hence `_flag_or_file`
        return normalize_schedule((;
            kind=:linear,
            burnin=_flag_or_file(args["burnin"], base.burnin),
            spacing=_flag_or_file(args["spacing"], base.spacing),
            n_samples=_flag_or_file(args["n-samples"], base.n_samples),
        ))
    end

    # nothing on the command line: the file, or a single sample
    return normalize_schedule(
        something(from_file, (; kind=:linear, burnin=nothing, spacing=nothing, n_samples=1))
    )
end

"""
    as_spacing_form(spec)

Restate a `linear` schedule given as `start`/`stop` in terms of `burnin`/`spacing`, so that a
flag can override one entry of it. The two describe the same times.
"""
function as_spacing_form(spec)
    hasproperty(spec, :spacing) && return spec
    spacing = spec.n_samples > 1 ? (spec.stop - spec.start) / (spec.n_samples - 1) : 0
    return (; kind=:linear, burnin=spec.start, spacing, n_samples=spec.n_samples)
end

_flag_or_file(flag, from_file) = isnothing(flag) ? from_file : flag

function _warn_ignored_table(from_file, flags)
    isnothing(from_file) && return nothing
    @warn """
    The `[schedule]` table of the parameter file is ignored: $flags was given on the command \
    line and describes the schedule differently.
    """
    return nothing
end

"""
    load_schedule(path) -> Union{Nothing,spec}

Read the `[schedule]` table of the parameter file at `path`, or return `nothing` if it has
none. See [`parse_tvals`](@ref).
"""
function load_schedule(path::AbstractString)
    table = get(TOML.parsefile(path), "schedule", nothing)
    return isnothing(table) ? nothing : parse_tvals(table; dir=dirname(path))
end

"""
    resolve_seed(seed) -> (seed, rng)

Return the rng seed to use and an rng seeded with it. If `seed` is `nothing`, one is drawn
from `RandomDevice`, so that a run is reproducible from the seed recorded in its output even
when `--seed` was not given. A local rng is returned rather than seeding the global one.
"""
function resolve_seed(seed::Union{Nothing,Integer})
    seed = isnothing(seed) ? rand(Random.RandomDevice(), UInt32) : seed
    return seed, Random.Xoshiro(seed)
end

"""
    compute_omega!(params, g, s0; rng, progress_meter)

Set `params.substitution_rate` to `average_transition_rate(g, params.step_type, s0)`, for
`--compute-omega`.

Does nothing if the parameter file already provides a rate: the file stays authoritative and
the flag is only a convenience for not having to compute `Ω` beforehand.
"""
function compute_omega!(params, g, s0; rng=Random.default_rng(), progress_meter=true)
    if !isnothing(params.substitution_rate)
        @warn """
        --compute-omega ignored: `substitution_rate` is already set to \
        $(params.substitution_rate) in the parameter file.
        """
        return params
    end
    @info "Computing the average substitution rate Ω..."
    params.substitution_rate = PottsEvolver.average_transition_rate(
        g, params.step_type, s0; rng, progress_meter, params=nothing
    )
    @info "Ω = $(params.substitution_rate)"
    return params
end

# Metadata recorded in the `[run]` table of `<prefix>.params.toml`, alongside the parameters.
function run_metadata(args, seed)
    return Dict{Symbol,Any}(
        :version => string(pkgversion(PottsEvolver)),
        :model => args["model"],
        :init => args["init"],
        :seed => Int(seed),
    )
end

# Shared by both subcommands: everything between parsing and calling `mcmc_sample`.
function setup_run(args)
    seed, rng = resolve_seed(args["seed"])
    @info "Using rng seed $seed"

    @info "Reading Potts model from $(args["model"])"
    g = read_graph(args["model"])

    @info "Reading parameters from $(args["params"])"
    params = load_params(args["params"])

    @info "Building the initial sequence from `$(args["init"])`"
    s0 = load_init(args["init"], args["init-type"], g; rng)

    args["compute-omega"] &&
        compute_omega!(params, g, s0; rng, progress_meter=args["verbose"] >= 1)

    return (; g, params, s0, rng, seed)
end

# `--translate` is a flag, so it is only forwarded when set: leaving it out keeps
# `mcmc_sample`'s own default, which is the same (`false`) for chains and trees.
_translate_kwargs(args) = args["translate"] ? (; translate_output=true) : (;)

function run_chain(args, sink)
    (; g, params, s0, rng, seed) = setup_run(args)
    schedule = resolve_schedule(args, args["params"])

    result = mcmc_sample(
        g,
        schedule_times(schedule),
        s0,
        params;
        rng,
        verbose=args["verbose"],
        extra_sinks=(sink,),
        progress_meter=args["verbose"] >= 1,
        _translate_kwargs(args)...,
    )

    return write_chain_output(
        args["outdir"],
        default_prefix(args["prefix"], args["model"]),
        result;
        run=run_metadata(args, seed),
        schedule,
    )
end

function run_tree(args, sink)
    (; g, params, s0, rng, seed) = setup_run(args)

    @info "Reading tree from $(args["tree"])"
    result = mcmc_sample(
        g,
        read_tree(args["tree"]),
        params;
        init=s0,
        burnin=something(args["burnin"], 0),
        rng,
        verbose=args["verbose"],
        extra_sinks=(sink,),
        _translate_kwargs(args)...,
    )

    metadata = run_metadata(args, seed)
    metadata[:tree] = args["tree"]
    return write_tree_output(
        args["outdir"],
        default_prefix(args["prefix"], args["model"]),
        result;
        run=metadata,
        internals=args["internals"],
    )
end

# Errors a user can cause: reported as a message, without a stacktrace. Anything else is a bug
# and is rethrown.
const _USER_ERRORS = Union{ArgumentError,SystemError,ArgParse.ArgParseError,KeyError}

# `showerror` on an `ArgParseError` prints its constructor, not its message, which is no use to
# someone who simply mistyped a flag.
_error_message(err) = sprint(showerror, err)
_error_message(err::ArgParse.ArgParseError) = err.text

"""
    (@main)(ARGS)

Entry point of the `pottsevolver` app. Return value is the process exit code: `0` on success,
`1` for an error caused by the input (bad file, invalid parameters, conflicting flags).
"""
function (@main)(ARGS)
    parsed = try
        parse_cli(ARGS)
    catch err
        err isa _USER_ERRORS || rethrow()
        # no logger yet: the output directory is not even known at this point
        println(stderr, "pottsevolver: $(_error_message(err))")
        println(stderr, ArgParse.usage_string(cli_settings()))
        println(stderr, "Run `pottsevolver --help` for the list of commands and options.")
        return 1
    end

    command = parsed["%COMMAND%"]
    args = parsed[command]

    prefix = default_prefix(args["prefix"], args["model"])
    logpath, io = setup_logfile(args["outdir"], prefix, command, ARGS)
    # one sink for the whole run: the CLI and `mcmc_sample` share this single handle
    sink = PottsEvolver.file_sink(io, args["log-verbose"])
    logger = PottsEvolver.get_logger(
        args["verbose"], nothing, args["log-verbose"]; extra_sinks=(sink,)
    )

    return try
        with_logger(logger) do
            try
                if command == "sample-chain"
                    run_chain(args, sink)
                else
                    run_tree(args, sink)
                end
                @info "Done. Log written to $logpath"
                return 0
            catch err
                if !(err isa _USER_ERRORS)
                    # not the user's fault: keep the stacktrace
                    @error "pottsevolver failed - this is likely a bug" exception = (
                        err, catch_backtrace()
                    )
                    rethrow()
                end
                # `@error` reaches both the console and the log file, so one message is enough
                @error "pottsevolver failed: $(_error_message(err))"
                return 1
            end
        end
    finally
        close(io)
    end
end

#================================#
########### CLI parsing ##########
#================================#

# Flags common to both subcommands: input files, output location, init sequence, seed,
# verbosity. Model parameters themselves live in the --params TOML file, not here.
function _add_shared_args!(settings::ArgParseSettings)
    @add_arg_table! settings begin
        "model"
        help = "Path to a Potts graph file (see `read_graph`)"
        required = true
        "--params", "-p"
        help = "Path to a TOML file with `SamplingParameters` fields (see `load_params`)"
        required = true
        "--outdir", "-o"
        help = "Output directory"
        default = pwd()
        "--prefix"
        help = "Output filename prefix (defaults to `model`'s basename without extension)"
        default = nothing
        "--init"
        help = "One of random_aa/random_codon/random_rna/random_num, or a fasta path \
                optionally suffixed with `,seqname` (see `load_init`)"
        default = "random_aa"
        "--init-type"
        help = "One of aa/codon/rna: force the sequence type read from an --init fasta \
                path, skipping inference"
        default = nothing
        "--burnin"
        help = "For sample-chain, time before the first sample (defaults to 5 * --spacing). \
                For sample-tree, time spent equilibrating the root sequence before sampling"
        arg_type = Float64
        default = nothing
        "--seed"
        help = "Integer RNG seed"
        arg_type = Int
        default = nothing
        "--verbose", "-v"
        help = "Verbosity: <0 error, 0 warn, 1 info, >=2 debug"
        arg_type = Int
        default = 0
        "--log-verbose"
        help = "Verbosity for the log file always written alongside the other outputs"
        arg_type = Int
        default = 1
        "--compute-omega"
        help = "If `substitution_rate` is missing from --params (continuous sampling \
                only), compute it via `average_transition_rate` before sampling"
        action = :store_true
        "--translate"
        help = "Translate the output to amino acids; only relevant when the initial/root \
                sequence is a `CodonSequence`. Off by default: the output keeps the type \
                that was sampled"
        action = :store_true
    end
    return settings
end

"""
    cli_settings() -> ArgParseSettings

The `ArgParseSettings` describing `pottsevolver`'s whole command line.

Separate from [`parse_cli`](@ref) so that `(@main)` can build the usage message to show
alongside a parse error: `exc_handler=debug_handler` makes `parse_args` throw instead of
printing usage and exiting, which is what lets the error be caught, but it also means nobody
prints the usage line unless we do.
"""
function cli_settings()
    settings = ArgParseSettings(; prog="pottsevolver", exc_handler=ArgParse.debug_handler)
    @add_arg_table! settings begin
        "sample-chain"
        action = :command
        help = "Sample sequences along a Markov chain"
        "sample-tree"
        action = :command
        help = "Sample sequences along the branches of a phylogenetic tree"
    end

    _add_shared_args!(settings["sample-chain"])
    @add_arg_table! settings["sample-chain"] begin
        "-M", "--n-samples"
        help = "Number of samples taken along the chain. Combines with --spacing/--burnin \
                and with the `[schedule]` table of --params; defaults to 1"
        arg_type = Int
        default = nothing # `nothing` means 1: it tells an explicit `-M 1` from no `-M` at all
        "--tvals"
        help = "Explicit sampling schedule (space-separated numbers), instead of \
                -M/--spacing/--burnin. The first value is the burn-in time; later values \
                are absolute times/steps, not deltas"
        nargs = '*'
        arg_type = Float64
        default = Float64[]
        "--spacing"
        help = "Time between samples (default 1, which is seldom enough to decorrelate). \
                Combines with -M/--burnin, and overrides the `[schedule]` table of --params"
        arg_type = Float64
        default = nothing
        "--tvals-file"
        help = "Path to a file with one time value per line, as an alternative to --tvals"
        default = nothing
    end

    _add_shared_args!(settings["sample-tree"])
    @add_arg_table! settings["sample-tree"] begin
        "tree"
        help = "Path to a Newick tree file"
        required = true
        "--internals"
        help = "Also write internal-node sequences and the (possibly relabeled) sampled tree"
        action = :store_true
    end

    return settings
end

"""
    parse_cli(args) -> Dict

Parse `pottsevolver`'s command-line arguments for the `sample-chain` and `sample-tree`
subcommands. Returns the nested `Dict` produced by `ArgParse.parse_args`: top-level
`"%COMMAND%"` names the chosen subcommand, and `d[d["%COMMAND%"]]` holds that subcommand's
own flags. Run `pottsevolver --help`, `pottsevolver sample-chain --help` or
`pottsevolver sample-tree --help` for the full, authoritative flag list.

A parse error raises an `ArgParse.ArgParseError` (rather than exiting the process), so it
can be caught and reported by the caller.
"""
parse_cli(args) = parse_args(args, cli_settings())
end # module CLI

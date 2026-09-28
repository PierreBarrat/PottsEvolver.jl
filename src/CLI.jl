module CLI

using ArgCheck
using ArgParse
using Dates
using FASTX
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

A `[run]` table, as written by [`write_params`](@ref), is ignored: it holds metadata about a
run (input files, seed, sequence type) rather than sampling parameters, so a
`<prefix>.params.toml` written by the CLI can be given back to it unchanged.

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
        if k == "run"
            # run metadata, not sampling parameters: see the docstring
            continue
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

    T = isnothing(init_type) ? _infer_sequence_type(s, g) : _sequence_type_from_name(init_type)
    return T(s)
end

function _split_init_path(init::AbstractString)
    parts = split(init, ","; limit=2)
    return length(parts) == 2 ? (String(parts[1]), String(parts[2])) : (String(parts[1]), nothing)
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
        throw(ArgumentError("Unknown `--init-type` `$init_type`. Expected `aa`, `codon` or `rna`."))
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
        throw(ArgumentError("""
            Cannot infer sequence type: content is ambiguous (only `-ACGT`) and its length \
            ($(length(s))) does not uniquely match the graph (q=$q, L=$L). Pass \
            `--init-type` explicitly (`aa`, `codon` or `rna`).
            """))
    else
        throw(ArgumentError("""
            Cannot infer sequence type: length $(length(s)) is compatible with more than one \
            candidate for this graph (q=$q, L=$L). Pass `--init-type` explicitly (`aa`, \
            `codon` or `rna`).
            """))
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

function _params_toml_dict(params; run=nothing)
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
"""
function write_params(path, params; run=nothing)
    open(path, "w") do io
        # `Symbol` is not a valid TOML type: `string` converts it, and `sequence_type` too
        TOML.print(string, io, _params_toml_dict(params; run); sorted=true)
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
function write_chain_output(outdir, prefix, result; run=nothing)
    files = String[]

    path = output_file(outdir, prefix, "fasta")
    write_fasta(path, result.sequences)
    push!(files, path)
    @info "Wrote $(length(result.sequences)) sequences to $path"

    push!(files, write_params(output_file(outdir, prefix, "params.toml"), result.params; run))

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

    push!(files, write_params(output_file(outdir, prefix, "params.toml"), result.params; run))

    if _tracks_substitutions(result.params)
        @warn """
        `track_substitutions` is set, but tree sampling does not return substitutions:
        no substitution file is written.
        """
    end
    return files
end

function (@main)(ARGS)
    parsed_args = parse_cli(ARGS)
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    return nothing
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
    parse_cli(args) -> Dict

Parse `pottsevolver`'s command-line arguments for the `sample-chain` and `sample-tree`
subcommands. Returns the nested `Dict` produced by `ArgParse.parse_args`: top-level
`"%COMMAND%"` names the chosen subcommand, and `d[d["%COMMAND%"]]` holds that subcommand's
own flags. Run `pottsevolver --help`, `pottsevolver sample-chain --help` or
`pottsevolver sample-tree --help` for the full, authoritative flag list.

A parse error raises an `ArgParse.ArgParseError` (rather than exiting the process), so it
can be caught and reported by the caller.
"""
function parse_cli(args)
    settings = ArgParseSettings(; exc_handler=ArgParse.debug_handler)
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
        help = "Number of samples taken along the chain (spaced by `Teq`, per --params)"
        arg_type = Int
        default = 1
        "--tvals"
        help = "Explicit sampling schedule (space-separated numbers), bypassing -M/Teq. \
                The first value is the burn-in time; later values are absolute \
                times/steps, not deltas"
        nargs = '*'
        arg_type = Float64
        default = Float64[]
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

    return parse_args(args, settings)
end
end # module CLI

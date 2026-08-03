module CLI

using ArgCheck
using ArgParse
using FASTX
using Random
using TOML
using ..PottsEvolver

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
        if k == "branchlength_meaning"
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

    file, seqname = _split_init_path(init)
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

    if !isempty(intersect(chars, _AA_EXCLUSIVE_SYMBOLS))
        @argcheck q == PottsEvolver.Q_AA && length(s) == L """
            Sequence contains symbol(s) exclusive to amino acids \
            ($(join(intersect(chars, _AA_EXCLUSIVE_SYMBOLS), ", "))), but the graph is \
            incompatible with `AASequence` (q=$q, L=$L, sequence length=$(length(s))). \
            Pass `--init-type` explicitly if this is intentional.
            """
        return AASequence
    elseif !isempty(intersect(chars, _RNA_EXCLUSIVE_SYMBOLS))
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

function (@main)(ARGS)
    parsed_args = parse_cli(ARGS)
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    return nothing
end

function parse_cli(args)
    tab = ArgParseSettings()
    @add_arg_table tab begin
        "sample-tree"
            help="Sample along branches of a tree"
            action=:command
        "sample-chain"
            help="Sample a chain"
            action=:command
    end
    return parse_args(args, tab)
end
end # module CLI

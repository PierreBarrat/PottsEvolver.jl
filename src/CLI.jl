module CLI

using ArgParse
using FASTX
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

`mutation_matrix`, if given, is a nested array **of rows**: `mutation_matrix[i]` is the
`i`-th row of the matrix, *i.e.* the relative mutation rates *out of* state `i` towards
every other state `j`. This matches `SamplingParameters`'s own `(i, a) -> (i, b)` rate
convention: `mutation_matrix[i][j]` becomes `μ[i, j]`.
```toml
mutation_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
```

Unknown keys or invalid field combinations surface as `SamplingParameters`'s own
`ArgumentError`.
"""
function load_params(path::AbstractString)
    raw = TOML.parsefile(path)
    kwargs = Dict{Symbol,Any}()

    symbol_fields = ("sampling_type", "step_type", "step_meaning")
    for (k, v) in raw
        k == "branchlength_meaning" && continue
        kwargs[Symbol(k)] = if k in symbol_fields
            Symbol(v)
        elseif k == "mutation_matrix"
            # `v` is a vector of rows; `hcat` stacks rows as columns, so `permutedims`
            # brings row `i` of the TOML array back to row `i` of the matrix.
            permutedims(reduce(hcat, v))
        else
            v
        end
    end

    if haskey(raw, "branchlength_meaning")
        blm = raw["branchlength_meaning"]
        kwargs[:branchlength_meaning] = BranchLengthMeaning(
            Symbol(blm["type"]), Symbol(blm["length"])
        )
    end

    return SamplingParameters(; kwargs...)
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

# just a wrapper to have <:TreeNodeData
@kwdef struct Sequence{T<:AbstractSequence} <: TreeNodeData
    seq::T = T(0)
end

Base.copy(S::Sequence) = Sequence(copy(S.seq))

"""
    mcmc_sample(g, tree, M=1, params; alignment_output, translate_output, kwargs...)

Sample `g` along branches of `tree`.
Repeat the process `M` times, returning an array of named tuples of the form
  `(; tree, leaf_sequences, internal_sequences)`.
If `M` is omitted, the output is just a named tuple (no array).

The sequence to be used as the root should be provided using the `init` kwarg,
  see `?PottsEvolver.get_init_sequence`.

If `alignment_output`, the sequences will be wrapped into an `Alignment` strucutre.
Otherwise, they are in a dictionary indexed by node label.
If `translate_output` and if the root sequence was a `CodonSequence`, the output alignment
will contain the amino acid sequence and not the codons.

## Warning
The `Teq` field of `params` is not used in the sampling.
However, the `burnin` field will be used to set the root sequence: `burnin` mcmc steps
will be performed starting from the input sequence, and the result is set at the root.
If you want a precise root sequence to be used, set `burnin=0` in `params`.
"""
function mcmc_sample(
    g, tree, params;
    alignment_output=true, translate_output=true, kwargs...
)
    # one sequence per node --> two alignments as output (+ tree)
    sampled_tree = mcmc_sample_tree(g, tree, params; kwargs...)
    leaf_sequences = map(n -> data(n).seq, leaves(sampled_tree))
    internal_sequences = map(n -> data(n).seq, internals(sampled_tree))

    return (;
        tree = sampled_tree,
        leaf_sequences = fmt_output(
            leaf_sequences, alignment_output, translate_output;
            names = map(label, leaves(sampled_tree)), dict=true,
        ),
        internal_sequences = fmt_output(
            internal_sequences, alignment_output, translate_output;
            names = map(label, internals(sampled_tree)), dict=true,
        )
    )
end
function mcmc_sample(g, tree,  M::Int, params; kwargs...)
    # M sequences per node --> [(tree, leaf, internals)] of length `M`
    return [mcmc_sample(g, tree, params; kwargs...) for _ in 1:M]
end


#=
Functions below just return a tree.
=#
"""
    mcmc_sample_tree(g::PottsGraph, tree::Tree, rootseq::AbstractSequence, params; kwargs...)
    mcmc_sample_tree(g::PottsGraph, tree::Tree, params::SamplingParameters; init, kwargs...)

Sample `g` along `tree`, using branch lengths to set the number of steps.
Return a sampled copy `tree`.
- If `rootseq` is supplied, use it as the root of the tree.
- If not, then the `init` kwarg will be used, see `PottsEvolver.get_init_sequence` to
  pick a root sequence. The code then falls back on the first form.
  The field `params.burnin` is then used to perform initial equilibration of the
    root sequence (useful if *e.g.* `init=:random_codon`)
"""
function mcmc_sample_tree(
    g::PottsGraph, tree::Tree, params::SamplingParameters;
    init=:random_num, verbose=0, kwargs...,
)
    # Pick initial sequence
    s0 = get_init_sequence(init, g; verbose)
    # If burnin, then equilibrate the picked sequence first
    @unpack burnin = params
    if burnin > 0
        verbose > 0 && println("Initializing with $(burnin) burnin iterations... ")
        gibbs_holder = get_gibbs_holder(s0)
        mcmc_steps!(s0, g, burnin, params; gibbs_holder, kwargs...)
    end
    # Sample and return the tree
    return mcmc_sample_tree(g, tree, s0, params; verbose, kwargs...)
end
function mcmc_sample_tree(
    g::PottsGraph, tree::Tree, rootseq::AbstractSequence, params::SamplingParameters;
    kwargs...
)
    tree_copy = prepare_tree(tree, rootseq)
    return mcmc_sample_tree_main!(g, tree_copy, params; kwargs...) # returns tree_copy
end
"""
    mcmc_sample_main!(
        g::PottsGraph, tree::Tree{<:Sequence{S}}, params::SamplingParameters;
        rng=Random.GLOBAL_RNG, verbose=0,
    ) where S <: AbstractSequence

Sample one sequence per node of `tree`.
Expects that `tree` is already "ready": it has a non empty root sequence, and branch
lengths are integers.
"""
function mcmc_sample_tree_main!(
    g::PottsGraph,
    tree::Tree{<:Sequence{S}},
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
    verbose=0,
) where {S<:AbstractSequence}
    # checks
    @argcheck all(n -> is_approx_integer(branch_length(n)), nodes(tree; skiproot=true)) """
        Branches of tree should be integers.\
        Instead $(map(branch_length, nodes(tree; skiproot=true)))
    """
    tmp_check_alphabet_consistency(g, data(root(tree)).seq)

    # logging
    if verbose > 0
        @info """
            Sampling $M sequences using the following settings:
                - Type of sequence = $(supertype(s0))
                - Steps between samples = $(Teq)
                - burnin = $(burnin)
                - Step style = $(params.step_type)
                - Branch length style = $(params.step_meaning)
                - fraction of gap steps (if codon) = $(params.fraction_gap_step)
            """
    end
    if params.Teq > 0 || params.burnin > 0
        @warn "Teq and burnin values from parameters will be ignored" params.Teq params.burnin
    end

    # some settings
    gibbs_holder = get_gibbs_holder(data(root(tree)).seq)

    # Sampling
    time = @elapsed sample_children!(root(tree), g, params, gibbs_holder; rng, verbose)
    verbose > 0 && @info "Sampling done in $time seconds"

    return tree
end

function sample_children!(node::TreeNode{<:Sequence}, g, params, gibbs_holder; kwargs...)
    for c in children(node)
        # copy the ancestral sequence
        s0 = copy(data(node).seq)
        # mcmc using it as an init
        nsteps = round(Int, branch_length(c))
        mcmc_steps!(s0, g, nsteps, params; gibbs_holder, kwargs...) # from sampling.jl
        # copy result to child
        data!(c, Sequence(s0))
        # recursive
        sample_children!(c, g, params, gibbs_holder; kwargs...)
    end
    return nothing
end

#=====================#
######## Utils ########
#=====================#

is_approx_integer(x) = isapprox(x, round(x))

function prepare_tree(tree::Tree, rootseq::S) where S <: AbstractSequence
    # convert to right type -- this makes a copy
    tree_copy = convert(Tree{Sequence{S}}, tree)
    # set root sequence
    data!(root(tree_copy), Sequence(copy(rootseq)))
    #
    round_tree_branch_length!(tree_copy)
    return tree_copy
end

function round_tree_branch_length!(tree::Tree; rtol=1e-3)
    for node in nodes(tree; skiproot=true)
        x = branch_length(node)
        if ismissing(x) || !isapprox(x, round(Int, x); rtol)
            throw(ArgumentError(
                "Expected tree with (approximately) integer branch lengths. Instead $x"
            ))
        else
            branch_length!(node, round(Int, x))
        end
    end
    return tree
end


"""
    pernode_alignment(data)

Transform data from an array of dictionaries to a dictionary of alignments.
If the trees are indexed by `m`, then this will transform `data[m].leaf_sequences[label]`
  to `new_data[label][m]`.

Use on the output of `mcmc_sample(graph, tree, M, params)`.
The latter function returns an array `data = [d1, d2, ...]` where each `d` is a named tuple
  of the form `d.tree, d.leaf_sequences, d.internal_sequences`.
Each `d` corresponds to one sampling run along `tree`, with `tree` being **shared** accross
    runs.

Transform this to a dictionary `label => alignment`, where `label` corresponds to tree nodes
  and `alignment` to the set of sequences obtained for this node in `data`.

**Note**: assume that all elements `d.tree` represent the same underlying tree (same labels).
"""
function pernode_alignment(data::AbstractVector{<:NamedTuple})
    if isempty(data)
        return (;tree=nothing, leaf_sequences=[], internal_sequences=nothing)
    end
    tree = first(data).tree
    leaf_sequences = _pernode_alignment([d.leaf_sequences for d in data])
    internal_sequences = _pernode_alignment([d.internal_sequences for d in data])
    sequences = merge(leaf_sequences, internal_sequences)
    return (; tree, sequences)
end

function _pernode_alignment(data::Vector{Dict{String,T}}) where T<:AbstractSequence
    out = Dict{String,Vector{T}}()
    for (m, tree_data) in enumerate(data), (label, seq) in tree_data
        push!(get!(out, label, T[]), seq)
    end
    return out
end

function _pernode_alignment(data::Vector{T}) where T<:Alignment
    @argcheck allequal(A -> A.alphabet, data)
    alphabet = first(data).alphabet
    labels = first(data).names
    L, N_nodes = size(first(data))
    M = length(data) # number of trees

    out = Dict{String,T}()
    for (i, label) in enumerate(labels)
        # construct data matrix for this node
        D = zeros(Int, L, M)
        for m in 1:M
            D[:, m] .= find_sequence(label, data[m])[2] # data[m] is an Alignment
        end
        out[label] = T(; data=D, alphabet, names=1:M)
    end
    return out
end

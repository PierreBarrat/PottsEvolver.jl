# just a wrapper to have <:TreeNodeData
@kwdef struct Sequence{T<:AbstractSequence} <: TreeNodeData
    seq::T = T(0) # sequence of null length as init.
end

Base.copy(S::Sequence) = Sequence(copy(S.seq))

#================================#
########### CONTINUOUS ###########
#================================#
"""
    mcmc_sample_continuous_tree(
        g::PottsGraph, tree::Tree, rootseq::AbstractSequence, params; kwargs...
    )
    mcmc_sample_continuous_tree(
        g::PottsGraph, tree::Tree, params::SamplingParameters; init, kwargs...
    )

Sample `g` along `tree`, using branch lengths to set the number of steps.
Return a sampled copy `tree`.
- If `rootseq` is supplied, use it as the root of the tree.
- If not, then the `init` kwarg will be used, see `PottsEvolver.get_init_sequence` to
  pick a root sequence.
  The field `params.burnin` is used to perform initial equilibration of the
  root sequence (useful if *e.g.* `init=:random_codon`).
  The code then falls back on the first form.
"""
function mcmc_sample_continuous_tree(
    g::PottsGraph, tree::Tree, params::SamplingParameters; init=:random_num, kwargs...
)
    # Pick initial sequence
    s0 = get_init_sequence(init, g)
    # If burnin, then equilibrate the picked sequence first
    @unpack burnin = params
    if burnin > 0
        @info "Equilibrating root sequence with $(burnin) burnin steps... "
        gibbs_holder = get_gibbs_holder(s0)
        time = @elapsed mcmc_steps!(s0, g, Float64(burnin), params; kwargs...) # Float(burnin) to ensure this goes to the continuous version
        @info "done in $time seconds"
    end

    # Sample and return the tree
    return mcmc_sample_continuous_tree(g, tree, s0, params; kwargs...)
end
function mcmc_sample_continuous_tree(
    g::PottsGraph,
    tree::Tree,
    rootseq::AbstractSequence,
    params::SamplingParameters;
    kwargs...,
)
    tree_copy = prepare_tree(tree, rootseq)
    return mcmc_sample_continuous_tree!(g, tree_copy, params; kwargs...) # returns tree_copy
end

"""
    mcmc_sample_continuous_tree!(g::PottsGraph, tree::Tree, params; kwargs...)

Sample one sequence per node of `tree`.
Expects `tree` to have a root sequence.
"""
function mcmc_sample_continuous_tree!(
    g::PottsGraph,
    tree::Tree{<:Sequence{S}},
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
) where {S<:AbstractSequence}
    tmp_check_alphabet_consistency(g, data(root(tree)).seq)
    @argcheck params.sampling_type == :continuous

    # logging & warnings
    let
        M = length(leaves(tree))
        s0 = data(root(tree)).seq
        @info """
        Sampling sequences of tree with $M leaves using the following settings:
            - type of sampling = `:continuous`
            - type of sequence = $(typeof(s0))
            - step style = $(params.step_type)
        """
    end

    @unpack Teq, burnin = params
    Teq > 0 && @info "Sampling on a tree: `Teq` field in parameters is ignored" Teq

    # Setting the CTMCState
    rootseq = data(root(tree)).seq
    state = CTMCState(rootseq)
    state.R = if !isnothing(params.substitution_rate)
        @info "Using provided average model substitution rate $(params.substitution_rate)"
        params.substitution_rate
    else
        @info "Computing average substitution rate for the model using discrete sampling..."
        (; value, time) = @timed average_transition_rate(g, params.step_type, rootseq; rng)
        @info "Done in $time seconds"
        @info "Average substitution rate: $value"
        value
    end

    # Sampling
    time = @elapsed sample_children_continuous!(root(tree), g, params, state; rng)
    @info "Sampling done in $time seconds"

    return tree
end

function sample_children_continuous!(
    node::TreeNode{<:Sequence}, g, params, state::CTMCState; kwargs...
)
    for c in children(node)
        # copy the ancestral sequence
        s0 = copy(data(node).seq)
        # mcmc using it as an init
        state.seq = s0 # state is uninitialised, this just avoids allocation
        state.previous_seq = nothing # to make sure not to take state seriously
        mcmc_steps!(state, g, branch_length(c), params.step_type; kwargs...)
        # copy result to child
        data!(c, Sequence(state.seq))
        # recursive
        sample_children_continuous!(c, g, params, state; kwargs...)
    end
    return nothing
end

#==========================#
######### DISCRETE #########
#==========================#
"""
    mcmc_sample_tree(g::PottsGraph, tree::Tree, rootseq::AbstractSequence, params; kwargs...)
    mcmc_sample_tree(g::PottsGraph, tree::Tree, params::SamplingParameters; init, kwargs...)

Sample `g` along `tree`, using branch lengths to set the number of steps.
Return a sampled copy `tree`.
- If `rootseq` is supplied, use it as the root of the tree.
- If not, then the `init` kwarg will be used, see `PottsEvolver.get_init_sequence` to
  pick a root sequence.
  The field `params.burnin` is used to perform initial equilibration of the
  root sequence (useful if *e.g.* `init=:random_codon`).
  The code then falls back on the first form.
"""
function mcmc_sample_tree(
    g::PottsGraph, tree::Tree, params::SamplingParameters; init=:random_num, kwargs...
)
    # Pick initial sequence
    s0 = get_init_sequence(init, g)
    # If burnin, then equilibrate the picked sequence first
    @unpack burnin = params
    if burnin > 0
        @info "Equilibrating root sequence with $(burnin) burnin steps... "
        gibbs_holder = get_gibbs_holder(s0)
        time = @elapsed mcmc_steps!(s0, g, burnin, params; gibbs_holder, kwargs...)
        @info "done in $time seconds"
    end

    # Sample and return the tree
    return mcmc_sample_tree(g, tree, s0, params; kwargs...)
end
function mcmc_sample_tree(
    g::PottsGraph,
    tree::Tree,
    rootseq::AbstractSequence,
    params::SamplingParameters;
    kwargs...,
)
    tree_copy = prepare_tree(tree, rootseq)
    return mcmc_sample_tree!(g, tree_copy, params; kwargs...) # returns tree_copy
end

"""
    mcmc_sample_tree!(
        g::PottsGraph, tree::Tree{<:Sequence{S}}, params::SamplingParameters;
        rng=Random.GLOBAL_RNG,
    ) where S <: AbstractSequence

Sample one sequence per node of `tree`, and return `tree`.
Expects a "prepared" tree with a non empty root sequence.
"""
function mcmc_sample_tree!(
    g::PottsGraph,
    tree::Tree{<:Sequence{S}},
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
) where {S<:AbstractSequence}
    tmp_check_alphabet_consistency(g, data(root(tree)).seq)
    @argcheck params.sampling_type == :discrete

    # logging & warnings
    let
        M = length(leaves(tree))
        s0 = data(root(tree)).seq
        @info """
        Sampling sequences of tree with $M leaves using the following settings:
            - type of sampling = `:discrete`
            - type of sequence = $(typeof(s0))
            - step style = $(params.step_type)
            - step meaning = $(params.step_meaning)
            - fraction of gap steps (if codon) = $(params.fraction_gap_step)
            - branch length meaning = $(params.branchlength_meaning)
        """
    end

    @unpack Teq, burnin = params
    Teq > 0 && @info "Sampling on a tree: `Teq` field in parameters is ignored" Teq

    # some settings
    gibbs_holder = get_gibbs_holder(data(root(tree)).seq)

    # Sampling
    time = @elapsed sample_children!(root(tree), g, params, gibbs_holder; rng)
    @info "Sampling done in $time seconds"

    return tree
end

function sample_children!(node::TreeNode{<:Sequence}, g, params, gibbs_holder; kwargs...)
    L = size(g).L
    for c in children(node)
        # copy the ancestral sequence
        s0 = copy(data(node).seq)
        # mcmc using it as an init
        nsteps = steps_from_branchlength(branch_length(c), params.branchlength_meaning, L)
        @debug """
        branchlen=$(branch_length(c)) - L=$L --> nsteps=$nsteps
        """
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

"""
    prepare_tree(tree::Tree, rootseq::S) where {S<:AbstractSequence}

Convert `tree` to an appropriate type for sampling, by attaching the proper datastructure to each node, to contain the sampled sequence.
"""
function prepare_tree(tree::Tree, rootseq::S) where {S<:AbstractSequence}
    # convert to right type -- this makes a copy
    tree_copy = convert(Tree{Sequence{S}}, tree)
    # set root sequence
    data!(root(tree_copy), Sequence(copy(rootseq)))

    return tree_copy
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
        return (; tree=nothing, leaf_sequences=[], internal_sequences=nothing)
    end
    tree = first(data).tree
    leaf_sequences = _pernode_alignment([d.leaf_sequences for d in data])
    internal_sequences = _pernode_alignment([d.internal_sequences for d in data])
    sequences = merge(leaf_sequences, internal_sequences)
    return (; tree, sequences)
end

function _pernode_alignment(data::Vector{Dict{String,T}}) where {T<:AbstractSequence}
    out = Dict{String,Vector{T}}()
    for (m, tree_data) in enumerate(data), (label, seq) in tree_data
        push!(get!(out, label, T[]), seq)
    end
    return out
end

function _pernode_alignment(data::Vector{T}) where {T<:Alignment}
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

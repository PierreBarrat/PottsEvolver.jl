#====================================================================================#
############################# mcmc_sample: chain version #############################
#====================================================================================#

"""
    mcmc_sample(
        g::PottsGraph, M::Integer, s0::AbstractSequence, params::SamplingParameters;
        rng=Random.GLOBAL_RNG, verbose=0, progress_meter=true, alignment_output=true,
    )

Sample `g` for `M` steps starting from `s0`, using parameters in `params`.
Return value: named tuple with fields
    - `sequences`: alignment (or vector) of sequences
    - `tvals`: vector with the number of steps at each sample
    - `info`: information about the run
    - `params`: parameters of the run.


*Note*: this function is not very efficient if `M` is small.

Sampling details are determined by `parameters`, see `?SamplingParameters`.
Whether to use the genetic code is determined by the type of the `init` sequence:
it is used if `init::CodonSequence`, otherwise not.
"""
function mcmc_sample end

function mcmc_sample(
    g::PottsGraph,
    M::Int,
    s0::AbstractSequence,
    params::SamplingParameters;
    verbose=0,
    logfile=nothing,
    logfile_verbose=1,
    kwargs...
)
    logger = get_logger(verbose, logfile, logfile_verbose)
    with_logger(logger) do
        mcmc_sample_chain(g, M, s0, params; kwargs...)
    end
end

"""
    mcmc_sample(
        g::PottsGraph, M::Integer, params::SamplingParameters; init=:random_num, kwargs...)
    )

Sample `g` for `M` steps, using parameters in `params`.
Choose the initial sequence based on `init` and `g`.
See `PottsEvolver.get_init_sequence` for more information.
"""
function mcmc_sample(
    g::PottsGraph, M::Integer, params::SamplingParameters;
    init=:random_num, verbose=0, kwargs...
)
    s0 = get_init_sequence(init, g; verbose)
    return mcmc_sample(g, M, s0, params; verbose, kwargs...)
end

#=================================================================================#
############################ mcmc_sample: tree version ############################
#=================================================================================#

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
    verbose=0,
    logfile=nothing,
    logfile_verbose=1,
    alignment_output=true,
    translate_output=true,
    kwargs...
)    # one sequence per node --> two alignments as output (+ tree)
    logger = get_logger(verbose, logfile, logfile_verbose)
    with_logger(logger) do
    # Actual MCMC
        sampled_tree = mcmc_sample_tree(g, tree, params; kwargs...)

        # Constructing output
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
end
function mcmc_sample(g, tree,  M::Int, params; kwargs...)
    # M sequences per node --> [(tree, leaf, internals)] of length `M`
    return [mcmc_sample(g, tree, params; kwargs...) for _ in 1:M]
end
function mcmc_sample(g::PottsGraph, tree::Tree, s0::AbstractSequence, params; kwargs...)
    return mcmc_sample(g, tree, params; init=s0, kwargs...)
end

#======================================================#
################### Initial sequence ###################
#======================================================#

"""
    get_init_sequence(s0, g::PottsGraph; kwargs...)

Try to guess a reasonable init sequence from `s0`:
- if `s0::AbstractSequence`, use a **copy** of it;
- if `s0::Symbol`, then it should be among `[:random_codon, :random_aa, :random_num]`;
  a random sequence of the corresponding type is created, using the length of `g`;
- if `s0` is a vector of integers, convert it to `AASequence`, `CodonSequence` or `NumSequence`;
  the conversion type depends on the alphabet of `g` and on the maximum element of `s0`.
"""
function get_init_sequence(s0::Symbol, g::PottsGraph; kwargs...)
    (; L, q) = size(g)
    return if s0 == :random_codon
        @argcheck q == length(aa_alphabet) """
            For sampling from `CodonSequence`, graph alphabet size must be $(length(aa_alphabet)).
            Instead $q.
            """
        CodonSequence(L)
    elseif s0 == :random_aa
        @argcheck q == length(aa_alphabet) """
            For sampling from `AASequence`, graph alphabet size must be $(length(aa_alphabet)).
            Instead $q.
            """
        AASequence(L)
    elseif s0 == :random_num
        NumSequence(L, q)
    else
        error(
            "Invalid symbol `init = $s0`. Options: `[:random_codon, :random_aa, :random_num]`",
        )
    end
end
get_init_sequence(s0::AbstractSequence, g; kwargs...) = copy(s0)
function get_init_sequence(s0::AbstractVector{<:Integer}, g; verbose=true)
    return if g.alphabet == aa_alphabet
        if maximum(s0) <= 21
            AASequence(s0)
        elseif 21 < maximum(s0) <= 65
            CodonSequence(s0)
        else
            error("Sequence $s0 incompatible with graph of size $(size(g))")
        end
    else
        NumSequence(s0)
    end
end

#=============================================#
################ Logging utils ################
#=============================================#

function get_logger(verbose, logfile, logfile_verbose)
    return if verbose < 0
        MinLevelLogger(ConsoleLogger(), Logging.Error)
    elseif verbose == 0
        MinLevelLogger(ConsoleLogger(), Logging.Warn)
    elseif verbose == 1
        MinLevelLogger(ConsoleLogger(), Logging.Info)
    else
        MinLevelLogger(ConsoleLogger(), Logging.Debug)
    end
end

#=====================#
######## Utils ########
#=====================#

function tmp_check_alphabet_consistency(g::PottsGraph, s0::CodonSequence)
    if isnothing(g.alphabet) || symbols(g.alphabet) != symbols(aa_alphabet)
        @warn """
            For now, sampling is only possible for graphs with the default alphabet $(aa_alphabet)
            Instead $(g.alphabet)
            """
    end
    return false
end
function tmp_check_alphabet_consistency(g::PottsGraph, s0::AASequence)
    if isnothing(g.alphabet) || symbols(g.alphabet) != symbols(aa_alphabet)
        @warn """
            For now, sampling is only possible for graphs with the default alphabet $(aa_alphabet)
            Instead $(g.alphabet)
            """
    end
    return false
end
tmp_check_alphabet_consistency(g::PottsGraph, s0::AbstractSequence) = true

function return_params(p::SamplingParameters, s::T) where {T<:AbstractSequence}
    d = Dict()
    for field in propertynames(p)
        d[field] = getproperty(p, field)
    end
    d[:sequence_type] = T
    return d
end

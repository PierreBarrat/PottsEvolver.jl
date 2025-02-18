#====================================================================================#
############################# mcmc_sample: chain version #############################
#====================================================================================#

"""
    mcmc_sample(
        g::PottsGraph, M::Integer, s0::AbstractSequence, params::SamplingParameters;
        rng=Random.GLOBAL_RNG, verbose=0, progress_meter=true, alignment_output=true,
    )
    mcmc_sample(
        g::PottsGraph, M::Integer, params::SamplingParameters; init=:random_num, kwargs...)
    )

First form: sample `g` for `M` steps starting from `s0`, using parameters in `params`.
Return value: named tuple with fields
- `sequences`: alignment (or vector) of sequences
- `tvals`: vector with the number of steps at each sample
- `info`: information about the run
- `params`: parameters of the run.

Second form: same, but initial sequence is provided through the `init` kwarg.


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
    kwargs...,
)
    logger = get_logger(verbose, logfile, logfile_verbose)
    with_logger(logger) do
        return if params.sampling_type == :continuous
            params = convert(SamplingParameters{FloatType}, params)
            mcmc_sample_continuous_chain(g, M, s0, params; kwargs...)
        elseif params.sampling_type == :discrete
            params = convert(SamplingParameters{Int}, params)
            mcmc_sample_chain(g, M, s0, params; kwargs...)
        else
            throw(ArgumentError("Invalid sampling type: $(params.sampling_type)"))
        end
    end
end
function mcmc_sample(
    g::PottsGraph,
    M::Integer,
    params::SamplingParameters;
    init=:random_num,
    verbose=0,
    kwargs...,
)
    s0 = get_init_sequence(init, g; verbose)
    return mcmc_sample(g, M, s0, params; verbose, kwargs...)
end

#=================================================================================#
############################ mcmc_sample: tree version ############################
#=================================================================================#

"""
    mcmc_sample(g, tree, M=1, params; alignment_output, translate_output, init, kwargs...)

Sample `g` along branches of `tree`.
Repeat the process `M` times, returning an array of named tuples of the form
  `(; tree, leaf_sequences, internal_sequences)`.
If `M` is omitted, the output is just a named tuple (no array).
Sequences in `leaf_sequences` and `internal_sequences` are sorted in post-order traversal.

The sequence to be used as the root should be provided using the `init` kwarg,
  see `?PottsEvolver.get_init_sequence`.

If `alignment_output`, the sequences will be wrapped into an `Alignment` structure.
Otherwise, they are in a dictionary indexed by node label.
If `translate_output` and if the root sequence was a `CodonSequence`, the output alignment
will contain the amino acid sequence and not the codons.

## Warning
The `Teq` field of `params` is not used in the sampling.
However, the `burnin` field will be used to set the root sequence: `burnin` mcmc steps
will be performed starting from the input sequence, and the result is placed at the root.
If you want a precise root sequence to be used, set `burnin=0` in `params`.
"""
function mcmc_sample(
    g,
    tree::Tree,
    params;
    verbose=0,
    logfile=nothing,
    logfile_verbose=1,
    alignment_output=true,
    translate_output=true,
    kwargs..., # init=get_init_sequence(...) here: passed to mcmc_sample_tree
)    # one sequence per node --> two alignments as output (+ tree)
    logger = get_logger(verbose, logfile, logfile_verbose)
    with_logger(logger) do
        # Actual MCMC
        sampled_tree = if params.sampling_type == :continuous
            params = convert(SamplingParameters{FloatType}, params)
            mcmc_sample_continuous_tree(g, tree, params; kwargs...)
        elseif params.sampling_type == :discrete
            params = convert(SamplingParameters{Int}, params)
            mcmc_sample_tree(g, tree, params; kwargs...)
        else
            throw(ArgumentError("Invalid sampling type: $(params.sampling_type)"))
        end

        leaf_names = map(label, traversal(sampled_tree, :postorder; internals=false))
        internal_names = map(label, traversal(sampled_tree, :postorder; leaves=false))

        # Constructing output
        leaf_sequences = map(n -> data(sampled_tree[n]).seq, leaf_names)
        internal_sequences = map(n -> data(sampled_tree[n]).seq, internal_names)

        return (;
            tree=sampled_tree,
            leaf_sequences=fmt_output(
                leaf_sequences,
                alignment_output,
                translate_output;
                names=leaf_names,
                dict=true,
            ),
            internal_sequences=fmt_output(
                internal_sequences,
                alignment_output,
                translate_output;
                names=internal_names,
                dict=true,
            ),
        )
    end
end
function mcmc_sample(g, tree::Tree, M::Int, params; kwargs...)
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
        NumSequence(s0, size(g).q)
    end
end

#=============================================#
################ Logging utils ################
#=============================================#

function get_logger(verbose, logfile, logfile_verbose)
    function min_lvl(val)
        return if val < 0
            Logging.Error
        elseif val == 0
            Logging.Warn
        elseif val == 1
            Logging.Info
        else
            Logging.Debug
        end
    end
    loggers = []
    # console logger
    console = MinLevelLogger(ConsoleLogger(), min_lvl(verbose))
    push!(loggers, console)
    # file logger
    file = if !isnothing(logfile) && !isempty(logfile)
        push!(loggers, MinLevelLogger(FileLogger(logfile), min_lvl(logfile_verbose)))
    end

    return TeeLogger(loggers...)
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

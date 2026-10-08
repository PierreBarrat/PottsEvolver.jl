#====================================================================================#
############################# mcmc_sample: chain version #############################
#====================================================================================#

"""
    mcmc_sample(
        g::PottsGraph, M::Integer, s0::AbstractSequence, params::SamplingParameters;
        spacing=0, burnin=5*spacing,
        rng=Random.GLOBAL_RNG, verbose=0, progress_meter=true, pack_output=true,
        translate_output=false,
    )
    mcmc_sample(
        g::PottsGraph, M::Integer, params::SamplingParameters; init=:random_num, kwargs...)
    )
    mcmc_sample(g, tvals::AbstractVector, s0, params::SamplingParameters; kwargs...)

First form: take `M` samples of `g` starting from sequence `s0`, using parameters in `params`.
Samples are `spacing` apart, and the first one is taken after `burnin`.
Return value: named tuple with fields
- `sequences`: a [`SequenceSample`](@ref), or a vector of sequences if `pack_output=false`
- `tvals`: vector with the time (number of steps, if discrete) at each sample
- `info`: information about the run
- `params`: parameters of the run.

Second form: same, but initial sequence is provided through the `init` kwarg.
See `?get_init_sequence` for details on how the initial sequence is determined from `init`.

Third form: give the sampling times explicitly instead of `M`/`spacing`/`burnin`. The first
value of `tvals` acts as the burn-in, and the following ones are absolute times, not gaps:
any `AbstractVector` will do, so `range`, `logrange` or a hand-written vector all work. Can
also be used with the `init` kwarg.

The schedule is *not* part of `params`: `SamplingParameters` describes how a step is taken,
while `spacing`/`burnin`/`tvals` say when samples are taken. For `:discrete` sampling, times
are numbers of mcmc steps and must be whole numbers.


If `pack_output`, the sequences are wrapped into a [`SequenceSample`](@ref), otherwise
they are returned as a vector.
If `translate_output` and if `s0` was a `CodonSequence`, the output sample will contain the
amino acid sequences and not the codons. It defaults to `false`: the output keeps the type
that was sampled, and `PottsEvolver.translate` can be applied to it afterwards.

*Note*: this function is not very efficient if `M` is small.

Sampling details are determined by `parameters`, see `?SamplingParameters`.
Whether to use the genetic code is determined by the type of the `init` sequence:
it is used if `init::CodonSequence`, otherwise not.
"""
function mcmc_sample end

function mcmc_sample(
    g::PottsGraph,
    M::Integer,
    s0::AbstractSequence,
    params;
    spacing=0,
    burnin=5 * spacing,
    kwargs...,
)
    @argcheck M > 0 "Number of samples `M` must be >0. Instead $M"

    @argcheck spacing >= 0 && burnin >= 0
    if params.sampling_type == :discrete
        @argcheck isinteger(spacing) && isinteger(burnin) """
        For :discrete mcmc, `spacing` and `burnin` are numbers of steps and must be integers.
        Instead $spacing and $burnin.
        """
        spacing, burnin = Int(spacing), Int(burnin)
    end
    tvals = if spacing > 0
        burnin .+ range(0, (M - 1) * spacing; step=spacing)
    else
        burnin .+ zeros(typeof(burnin), M)
    end
    if params.sampling_type == :continuous
        tvals = Float64.(tvals)
    end

    return mcmc_sample(g, tvals, s0, params; kwargs...)
end
function mcmc_sample(
    g::PottsGraph,
    M::Integer,
    params::SamplingParameters;
    init=:random_num,
    rng=Random.default_rng(),
    verbose=0,
    kwargs...,
)
    s0 = get_init_sequence(init, g; rng)
    return mcmc_sample(g, M, s0, params; verbose, rng, kwargs...)
end

function mcmc_sample(
    g::PottsGraph,
    tvals::AbstractVector,
    s0::AbstractSequence,
    params::SamplingParameters;
    verbose=0,
    logfile=nothing,
    logfile_verbose=1,
    extra_sinks=(),
    kwargs...,
)
    logger = get_logger(verbose, logfile, logfile_verbose; extra_sinks)
    with_logger(logger) do
        warn_repeated_times(tvals)
        return if params.sampling_type == :continuous
            mcmc_sample_continuous_chain(g, Float64.(tvals), s0, params; kwargs...)
        elseif params.sampling_type == :discrete
            mcmc_sample_chain(g, discrete_times(tvals), s0, params; kwargs...)
        else
            throw(ArgumentError("Invalid sampling type: $(params.sampling_type)"))
        end
    end
end
function mcmc_sample(
    g::PottsGraph,
    tvals::AbstractVector,
    params::SamplingParameters;
    init=:random_num,
    rng=Random.default_rng(),
    verbose=0,
    kwargs...,
)
    s0 = get_init_sequence(init, g; verbose, rng)
    return mcmc_sample(g, tvals, s0, params; verbose, kwargs...)
end

#=================================================================================#
############################ mcmc_sample: tree version ############################
#=================================================================================#

"""
    mcmc_sample(g, tree, M=1, params; pack_output, translate_output, init, kwargs...)

Sample `g` along branches of `tree`.
Repeat the process `M` times, returning an array of named tuples of the form
  `(; tree, leaf_sequences, internal_sequences)`.
If `M` is omitted, the output is just a named tuple (no array).
Sequences in `leaf_sequences` and `internal_sequences` are sorted in post-order traversal.

The sequence to be used as the root should be provided using the `init` kwarg,
  see `?PottsEvolver.get_init_sequence`.

If `pack_output`, the sequences will be wrapped into a [`SequenceSample`](@ref).
Otherwise, they are in a dictionary indexed by node label.
If `translate_output` and if the root sequence was a `CodonSequence`, the output sample
will contain the amino acid sequence and not the codons.
`translate_output` defaults to `false`: the output keeps the type that was sampled, and
`PottsEvolver.translate` can be applied to it afterwards.

Sampling times come from the branch lengths of the tree, so there is no schedule to give.
The `burnin` kwarg sets the root sequence: `burnin` mcmc steps are performed starting from the
input sequence, and the result is placed at the root. It defaults to 0, which places the input
sequence at the root exactly as given.
"""
function mcmc_sample(
    g,
    tree::Tree,
    params;
    verbose=0,
    logfile=nothing,
    logfile_verbose=1,
    extra_sinks=(),
    pack_output=true,
    translate_output=false,
    kwargs..., # init=get_init_sequence(...) here: passed to mcmc_sample_tree
)    # one sequence per node --> two alignments as output (+ tree)
    logger = get_logger(verbose, logfile, logfile_verbose; extra_sinks)
    with_logger(logger) do
        # Actual MCMC
        sampled_tree = if params.sampling_type == :continuous
            mcmc_sample_continuous_tree(g, tree, params; kwargs...)
        elseif params.sampling_type == :discrete
            mcmc_sample_tree(g, tree, params; kwargs...)
        else
            throw(ArgumentError("Invalid sampling type: $(params.sampling_type)"))
        end

        leaf_names = map(label, traversal(sampled_tree, :postorder; internals=false))
        internal_names = map(label, traversal(sampled_tree, :postorder; leaves=false))

        # Constructing output
        leaf_sequences = map(n -> data(sampled_tree[n]).seq, leaf_names)
        internal_sequences = map(n -> data(sampled_tree[n]).seq, internal_names)
        params = return_params(params, root(sampled_tree).data.seq)
        return (;
            tree=sampled_tree,
            leaf_sequences=fmt_output(
                leaf_sequences,
                pack_output,
                translate_output;
                names=leaf_names,
                dict=true,
            ),
            internal_sequences=fmt_output(
                internal_sequences,
                pack_output,
                translate_output;
                names=internal_names,
                dict=true,
            ),
            params,
        )
    end
end
function mcmc_sample(g, tree::AbstractString, params; kwargs...)
    # read tree from a file
    return mcmc_sample(g, read_tree(tree), params; kwargs...)
end
function mcmc_sample(g, tree::Tree, M::Int, params; kwargs...)
    # M sequences per node --> [(tree, leaf, internals)] of length `M`
    return [mcmc_sample(g, tree, params; kwargs...) for _ in 1:M]
end
function mcmc_sample(g::PottsGraph, tree::Tree, s0::AbstractSequence, params; kwargs...)
    # init sequence as positional arg
    return mcmc_sample(g, tree, params; init=s0, kwargs...)
end
function mcmc_sample(g, tree::AbstractString, x::Any, params; kwargs...)
    tree = read_tree(tree)
    return mcmc_sample(g, tree, x, params; kwargs...)
end

#======================================================#
################### Initial sequence ###################
#======================================================#

"""
    get_init_sequence(s0, g::PottsGraph; kwargs...)

Try to guess a reasonable init sequence from `s0`:
- if `s0::AbstractSequence`, use a **copy** of it;
- if `s0::Symbol`, then it should be among `[:random_codon, :random_aa, :random_rna, :random_num]`;
  a random sequence of the corresponding type is created, using the length of `g`;
- if `s0` is a vector of integers, convert it to `AASequence`, `CodonSequence` or `NumSequence`;
  the conversion type depends on the number of states `q` of `g` and on the maximum element
  of `s0`. Graphs with `q == 21` are assumed to represent amino acids.

The graph `g` is only used to determine the length of the sequence, and the alphabet size in the case of a numerical sequence.
"""
function get_init_sequence(s0::Symbol, g::PottsGraph; rng=Random.default_rng(), kwargs...)
    (; L, q) = size(g)
    return if s0 == :random_codon
        @argcheck q == Q_AA """
            For sampling from `CodonSequence`, graph alphabet size must be $(Q_AA).
            Instead $q.
            """
        CodonSequence(rng, L)
    elseif s0 == :random_aa
        @argcheck q == Q_AA """
            For sampling from `AASequence`, graph alphabet size must be $(Q_AA).
            Instead $q.
            """
        AASequence(rng, L)
    elseif s0 == :random_rna
        @argcheck q == Q_RNA """
            For sampling from `RNASequence`, graph alphabet size must be $(Q_RNA).
            Instead $q.
            """
        RNASequence(rng, L)
    elseif s0 == :random_num
        NumSequence(rng, L, q)
    else
        error(
            "Invalid symbol `init = $s0`. Options: `[:random_codon, :random_aa, :random_rna, :random_num]`",
        )
    end
end
get_init_sequence(s0::AbstractSequence, g; kwargs...) = copy(s0)
function get_init_sequence(
    s0::AbstractVector{<:Integer}, g; rng=Random.default_rng(), kwargs...
)
    return if size(g).q == Q_AA
        if maximum(s0) <= 21
            AASequence(s0)
        elseif 21 < maximum(s0) <= 65
            CodonSequence(s0; source=:codon, rng)
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

"""
    log_level(verbose)

`LogLevel` corresponding to a verbosity: `<0` error, `0` warn, `1` info, `>=2` debug.
"""
function log_level(verbose)
    return if verbose < 0
        Logging.Error
    elseif verbose == 0
        Logging.Warn
    elseif verbose == 1
        Logging.Info
    else
        Logging.Debug
    end
end

"""
    file_sink(io::IO, verbose)

A logger writing to the already open stream `io` at verbosity `verbose`.

Pass it in `extra_sinks` to log to a file that the *caller* owns and closes, instead of having
[`get_logger`](@ref) open one itself through its `logfile` argument. This is the only safe way
for several loggers to write to the same file: Julia's append mode seeks to the end of the file
when *opening* it, and each handle then keeps its own position, so two handles onto one file
overwrite each other's output. Sharing a single stream leaves a single position.
"""
file_sink(io::IO, verbose) = MinLevelLogger(FileLogger(io), log_level(verbose))

"""
    get_logger(verbose, logfile, logfile_verbose; extra_sinks=())

Build a `TeeLogger` writing to the console at verbosity `verbose`, and to `logfile` (if given,
as a path) at verbosity `logfile_verbose`. Loggers in `extra_sinks` are added to it, see
[`file_sink`](@ref).
"""
function get_logger(verbose, logfile, logfile_verbose; extra_sinks=())
    loggers = []
    # console logger
    console = MinLevelLogger(ConsoleLogger(Logging.Debug), log_level(verbose))
    push!(loggers, console)
    # file logger - `logfile` is a path, opened (and never closed) here
    file = if !isnothing(logfile) && !isempty(logfile)
        push!(loggers, MinLevelLogger(FileLogger(logfile), log_level(logfile_verbose)))
    end
    append!(loggers, extra_sinks)

    return TeeLogger(loggers...)
end

#=====================#
######## Utils ########
#=====================#

"""
    discrete_times(tvals)

Convert a sampling schedule to integers, for discrete sampling.

Discrete time is a number of mcmc steps, so the schedule has to be whole numbers - and
`mcmc_sample_chain` dispatches on `AbstractVector{<:Integer}`, which would otherwise turn a
`Float64` schedule into a `MethodError` rather than a usable message.
"""
function discrete_times(tvals)
    fractional = filter(!isinteger, tvals)
    @argcheck isempty(fractional) """
    For :discrete mcmc, sampling times are numbers of steps and must be whole numbers.
    Instead got $(fractional).
    """
    return map(Int, tvals)
end

"""
    warn_repeated_times(tvals)

Warn if a sampling schedule asks for several samples at the same time.

No mcmc step happens between two samples taken at the same time, so they come out as identical
sequences. The usual cause is a `spacing` of 0, which is still the default: it is a valid
request (*e.g.* `M` copies of the initial sequence) but rarely what is wanted.
"""
function warn_repeated_times(tvals)
    repeats = length(tvals) - length(unique(tvals))
    repeats > 0 && @warn """
    The sampling schedule asks for the same time more than once ($repeats repeated value(s)).
    No mcmc step happens between samples taken at the same time, so the corresponding
    sequences will be identical. Use a non-zero `spacing`, or pass distinct sampling times.
    """
    return nothing
end

function return_params(p::SamplingParameters, ::T) where {T<:AbstractSequence}
    d = Dict()
    for field in propertynames(p)
        d[field] = getproperty(p, field)
    end
    d[:sequence_type] = T
    return d
end

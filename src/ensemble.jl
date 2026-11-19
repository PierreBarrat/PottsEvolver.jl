#====================================================================================#
################################## Parallel chains ###################################
#====================================================================================#

"""
    parallel_chains(
        g::PottsGraph,
        schedule::Union{Integer,AbstractVector},
        inits::AbstractVector{<:AbstractSequence},
        params::SamplingParameters;
        seed=rand(UInt32), verbose=0, kwargs...,
    )
    parallel_chains(g, schedule, params; n_chains, init=:random_num, kwargs...)

Run one independent chain per element of `inits`, in parallel over the available threads.

Return a `Vector{SequenceSample}` with one entry per chain, in the order of `inits`. Each
entry is the trajectory of that chain, labelled by sampling time, exactly as the `sequences`
field of `mcmc_sample` — which is what each chain is run through.

`schedule` is passed on to `mcmc_sample`: either a number `M` of samples, with the `spacing`
and `burnin` kwargs, or a vector of sampling times. Remaining `kwargs` go to `mcmc_sample`
too, except that `pack_output` is fixed and `progress_meter` defaults to `false`: with
several chains running at once the bars are of little use.

Chain `i` is run with `Xoshiro(seed + i)`, so a run is reproducible from `seed` alone
and gives the same result whether or not threads are available. A seed is drawn if none
is given.

The second form draws `n_chains` initial sequences from `init`, see
`PottsEvolver.get_init_sequence`; they follow from `seed` as well, so it reproduces the
whole run. Pass `inits` instead to control them yourself.

To regroup the output by time rather than by chain, pass it to
`PottsEvolver.regroup_samples`. The special case of one sample per chain is
`equilibrium_sample`.
"""
function parallel_chains(
    g::PottsGraph,
    schedule::Union{Integer,AbstractVector},
    inits::AbstractVector{<:AbstractSequence},
    params::SamplingParameters;
    seed=rand(UInt32),
    verbose=0,
    progress_meter=false,
    kwargs...,
)
    @argcheck !isempty(inits) "Need at least one initial sequence"

    with_logger(get_logger(verbose, nothing, 0)) do
        if Threads.nthreads() == 1
            @info "Only one thread available: running $(length(inits)) chains sequentially."
        else
            @info "Running $(length(inits)) chains on $(Threads.nthreads()) threads."
        end
    end

    out = Vector{SequenceSample}(undef, length(inits))
    # `g` and `params` are only read: sampling copies the initial sequence and never writes
    # to the parameters, so the chains share them safely.
    if Threads.nthreads() == 1
        for (i, s0) in enumerate(inits)
            out[i] = _one_chain(
                g, schedule, s0, params, seed + i; verbose, progress_meter, kwargs...
            )
        end
    else
        Threads.@threads for i in eachindex(inits)
            out[i] = _one_chain(
                g, schedule, inits[i], params, seed + i; verbose, progress_meter, kwargs...
            )
        end
    end
    return identity.(out) # narrow the element type now that all chains are in
end

function parallel_chains(
    g::PottsGraph,
    schedule::Union{Integer,AbstractVector},
    params::SamplingParameters;
    n_chains,
    init=:random_num,
    seed=rand(UInt32),
    kwargs...,
)
    @argcheck n_chains > 0 "Number of chains must be >0. Instead $n_chains"
    # one rng advanced across the draws, so the initial sequences differ from each other but
    # still follow from `seed`; the chains themselves use `seed + i`, a different stream.
    rng = Random.Xoshiro(seed)
    inits = [get_init_sequence(init, g; rng) for _ in 1:n_chains]
    return parallel_chains(g, schedule, inits, params; seed, kwargs...)
end

function _one_chain(g, schedule, s0, params, seed; kwargs...)
    return mcmc_sample(
        g, schedule, s0, params; rng=Random.Xoshiro(seed), pack_output=true, kwargs...
    ).sequences
end

"""
    equilibrium_sample(
        g::PottsGraph, n_chains::Integer, params::SamplingParameters;
        burnin, init=:random_num, kwargs...,
    )

Draw `n_chains` sequences from `g`, one per independent chain.

Each chain starts from its own sequence drawn from `init` and is sampled once, after
`burnin` steps. Return a single `SequenceSample` of `n_chains` sequences, labelled by
chain index.

Taking one sequence per chain buys independence between samples outright, rather than
relying on a `spacing` large enough to decorrelate a single chain — and starting from
several random sequences avoids reading one trapped chain as the equilibrium. That makes
it the way to get the sample that `average_transition_rate` needs.

`burnin` has no default: how long a chain needs to equilibrate is a property of the
model, and the whole point of the sample is that it is drawn at equilibrium. Other
`kwargs` are passed to `parallel_chains`.
"""
function equilibrium_sample(
    g::PottsGraph,
    n_chains::Integer,
    params::SamplingParameters;
    burnin,
    init=:random_num,
    kwargs...,
)
    @argcheck burnin >= 0 "`burnin` must be >=0. Instead $burnin"
    chains = parallel_chains(g, [burnin], params; n_chains, init, kwargs...)
    return only(regroup_samples(chains))
end

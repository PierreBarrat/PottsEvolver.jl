"""
    parallel_chains(
        g::PottsGraph,
        time_steps::AbstractVector{<:Integer},
        inits::AbstractVector{<:AbstractSequence},
        params::SamplingParameters;
        verbose=0,
        kwargs...,
    )
"""
function _sample_chains_sequential(
    g::PottsGraph,
    time_steps::AbstractVector{<:Integer},
    inits::AbstractVector{T},
    params::SamplingParameters,
    logger;
    kwargs...,
) where {T<:AbstractSequence}
    n_chains = length(inits)
    S = Vector{Vector{T}}(undef, n_chains)
    with_logger(logger) do
        for i in eachindex(inits)
            rng = Random.default_rng()
            (; sequences) = PottsEvolver.mcmc_sample_chain(
                g,
                time_steps,
                inits[i],
                params;
                rng,
                alignment_output=false,
                progress_meter=false,
                store_info=false,
                kwargs...,
            )
            S[i] = sequences
        end
    end
    return S
end

function _sample_chains_parallel(
    g::PottsGraph,
    time_steps::AbstractVector{<:Integer},
    inits::AbstractVector{T},
    params::SamplingParameters,
    logger;
    kwargs...,
) where {T<:AbstractSequence}
    n_chains = length(inits)
    S = Vector{Vector{T}}(undef, n_chains)
    with_logger(logger) do
        Threads.@threads for i in eachindex(inits)
            rng = Random.default_rng()
            (; sequences) = PottsEvolver.mcmc_sample_chain(
                g,
                time_steps,
                inits[i],
                params;
                rng,
                alignment_output=false,
                progress_meter=false,
                store_info=false,
                kwargs...,
            )
            S[i] = sequences
        end
    end
    return S
end

function parallel_chains(
    g::PottsGraph,
    time_steps::AbstractVector{<:Integer},
    inits::T,
    params::SamplingParameters;
    verbose=0,
    kwargs...,
) where {T<:AbstractVector{<:AbstractSequence}}
    logger = PottsEvolver.get_logger(verbose, nothing, 0)

    if Threads.nthreads() == 1
        verbose > 0 && @info "Only one thread available. Running sequentially."
        return _sample_chains_sequential(g, time_steps, inits, params, logger; kwargs...)
    else
        verbose > 0 && @info "Running in parallel with $(Threads.nthreads()) threads."
        return _sample_chains_parallel(g, time_steps, inits, params, logger; kwargs...)
    end
end

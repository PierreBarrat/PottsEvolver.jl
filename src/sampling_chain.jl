#==========================#
######### DISCRETE #########
#==========================#

function mcmc_sample_chain(
    g::PottsGraph,
    time_steps::AbstractVector{<:Integer},
    s0::AbstractSequence,
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
    progress_meter=true,
    alignment_output=true,
    translate_output=false,
)
    # Argument checks
    @argcheck params.sampling_type == :discrete
    M = length(time_steps)
    @argcheck M > 0 "Number of samples must be >0. Instead $(M)"
    @argcheck issorted(time_steps) && all(>=(0), time_steps) """
    Time steps must be positive and sorted in ascending order. Instead $time_steps
    """
    tmp_check_alphabet_consistency(g, s0)
    if !alignment_output && translate_output
        error("I have to implement this case")
    end

    @info """
    Sampling $M sequences using the following settings:
        - Type of sampling: $(params.sampling_type)
        - Type of sequence = $(typeof(s0))
        - Step style = $(params.step_type)
        - Step meaning = $(params.step_meaning)
        - Fraction of gap steps (if codon) = $(params.fraction_gap_step)
        - Initial/last time step: $(first(time_steps))/$(last(time_steps))
    """
    @info "Initial sequence: $s0"

    # vector for return value
    conf = copy(s0)
    S = similar([conf], M) # the sample
    # tvals = Vector{Int}(undef, M) # number of steps at each sample

    # Holder if gibbs steps are used
    # the size of the holder depends on the symbols of the sequence: amino acids or codons
    gibbs_holder = get_gibbs_holder(s0)

    # Burnin
    burnin = first(time_steps)
    @info "Initializing with $(burnin) burnin iterations... "
    burnin > 0 && mcmc_steps!(conf, g, burnin, params; rng, gibbs_holder)
    S[1] = copy(conf)

    # Sampling
    log_info = []
    progress = Progress(
        M - 1;
        barglyphs=BarGlyphs("[=> ]"),
        desc="Sampling: ",
        showspeed=true,
        enabled=progress_meter,
        dt=2,
    )
    @info "Sampling..."
    t_previous = first(time_steps)
    time = @elapsed for (m, t_now) in enumerate(Iterators.drop(time_steps, 1))
        T = t_now - t_previous
        # doing T steps on the current configuration
        @debug "Doing $T discrete steps"
        _, proposed, performed = mcmc_steps!(conf, g, T, params; rng, gibbs_holder)
        # storing the result in S
        S[m+1] = copy(conf)
        # misc.
        push!(
            log_info, (proposed=proposed, performed=performed, ratio=performed / proposed)
        )
        t_previous = t_now
        next!(progress; showvalues=[("steps", m + 1), ("total", M)])
    end
    @info "Sampling done in $time seconds"

    sequences = fmt_output(S, alignment_output, translate_output; names=time_steps)
    return (;
        sequences,
        tvals=collect(time_steps),
        info=log_info,
        params=return_params(params, s0)
    )
end

#=================================#
########### CONTINUOUS ############
#=================================#

function mcmc_sample_continuous_chain(
    g::PottsGraph,
    time_steps::AbstractVector{<:AbstractFloat},
    s0::AbstractSequence,
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
    progress_meter=true,
    alignment_output=true,
    translate_output=false,
)
    # Argument checks
    @argcheck params.sampling_type == :continuous
    M = length(time_steps)
    @argcheck M > 0 "Number of samples must be >0. Instead $(M)"
    @argcheck issorted(time_steps) && all(>=(0), time_steps) """
    Time steps must be positive and sorted in ascending order. Instead $time_steps
    """
    tmp_check_alphabet_consistency(g, s0)
    if !alignment_output && translate_output
        error("I have to implement this case")
    end

    @unpack step_type = params

    @info """
       Sampling $M sequences using the following settings:
        - Type of sampling: $(params.sampling_type)
        - Type of sequence = $(typeof(s0))
        - Step style = $(params.step_type)
        - Step meaning = $(params.step_meaning)
        - Fraction of gap steps (if codon) = $(params.fraction_gap_step)
        - Initial/last time step: $(first(time_steps))/$(last(time_steps))
    """

    @info "Initial sequence: $s0"

    # vector for return value
    conf = copy(s0)
    S = similar([conf], M) # the sample

    R = if !isnothing(params.substitution_rate)
        @info "Using provided average model substitution rate $(params.substitution_rate)"
        params.substitution_rate
    else
        @info "Computing average substitution rate for the model using discrete sampling..."
        (; value, time) = @timed average_transition_rate(
            g, step_type, s0; rng, progress_meter
        )
        @info "Done in $time seconds"
        @info "Average substitution rate: $value"
        value
    end

    # State structure for the sampling
    state = CTMCState(conf)
    state.R = R

    # Burnin
    burnin = first(time_steps)
    @info "Initializing with $(burnin) burnin iterations... "
    burnin > 0 && mcmc_steps!(state, g, burnin, step_type; rng)
    S[1] = copy(state.seq)

    # Sampling
    log_info = []
    progress = Progress(
        M - 1;
        barglyphs=BarGlyphs("[=> ]"),
        desc="Sampling: ",
        showspeed=true,
        enabled=progress_meter,
        dt=2,
    )
    @info "Sampling..."
    t_previous = first(time_steps)
    time = @elapsed for (m, t_now) in enumerate(Iterators.drop(time_steps, 1))
        T = t_now - t_previous
        # doing T steps on the current configuration
        @debug "Sampling for time $T"
        _, number_substitutions, substitutions = mcmc_steps!(
            state, g, T, step_type; rng, params.track_substitutions,
        )
        # storing the result in S
        S[m+1] = copy(state.seq)
        # misc.
        substitutions = if params.track_substitutions
            # add the substitutions of the last T steps
            substitutions[(end-number_substitutions+1):end]
        else
            # not tracking substitutions
            Tuple{Mutation, Float64}[]
        end
        new_info = (;
            number_substitutions,
            substitutions,
        )
        push!(log_info, new_info)
        t_previous = t_now
        next!(progress; showvalues=[("steps", m + 1), ("total", M)])
    end
    @info "Sampling done in $time seconds"
    sequences = fmt_output(S, alignment_output, translate_output; names=time_steps)
    return (;
        sequences,
        tvals=collect(time_steps),
        info=log_info,
        params=return_params(params, s0)
    )
end

#=====================#
######## Utils ########
#=====================#

function fmt_output(
    sequences::AbstractVector{T}, alignment, translate; names=nothing, dict=false
) where {T<:CodonSequence}
    return if alignment
        A = Alignment(sequences; names)
        translate ? genetic_code(A) : A
    elseif dict
        Dict{String,T}(name => seq for (name, seq) in zip(names, sequences))
    else
        sequences
    end
end
function fmt_output(
    sequences::AbstractVector{T}, alignment, translate; names=nothing, dict=false
) where {T<:AbstractSequence}
    return if alignment
        Alignment(sequences; names)
    elseif dict
        Dict{String,T}(name => seq for (name, seq) in zip(names, sequences))
    else
        sequences
    end
end

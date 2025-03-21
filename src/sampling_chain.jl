#==========================#
######### DISCRETE #########
#==========================#

function mcmc_sample_chain(
    g::PottsGraph,
    M::Integer,
    s0::AbstractSequence,
    params::SamplingParameters{<:Integer};
    rng=Random.GLOBAL_RNG,
    progress_meter=true,
    alignment_output=true,
    translate_output=false,
)
    # Argument checks
    @argcheck params.sampling_type == :discrete
    @argcheck M > 0 "Number of samples `M` must be >0. Instead $M"
    tmp_check_alphabet_consistency(g, s0)
    if !alignment_output && translate_output
        error("I have to implement this case")
    end

    @unpack Teq, burnin = params

    @info """
    Sampling $M sequences using the following settings:
        - Type of sampling: $(params.sampling_type)
        - Type of sequence = $(typeof(s0))
        - Steps between samples = $(Teq)
        - burnin = $(burnin)
        - Step style = $(params.step_type)
        - Step meaning = $(params.step_meaning)
        - fraction of gap steps (if codon) = $(params.fraction_gap_step)
    """

    @info "Initial sequence: $s0"

    # vector for return value
    conf = copy(s0)
    S = similar([conf], M) # the sample
    tvals = Vector{Int}(undef, M) # number of steps at each sample

    # Holder if gibbs steps are used
    # the size of the holder depends on the symbols of the sequence: amino acids or codons
    gibbs_holder = get_gibbs_holder(s0)

    # Burnin
    @info "Initializing with $(burnin) burnin iterations... "
    burnin > 0 && mcmc_steps!(conf, g, burnin, params; rng, gibbs_holder)
    S[1] = copy(conf)
    tvals[1] = burnin

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
    time = @elapsed for m in 2:M
        # doing Teq steps on the current configuration
        _, proposed, performed = mcmc_steps!(conf, g, Teq, params; rng, gibbs_holder)
        # storing the result in S
        S[m] = copy(conf)
        tvals[m] = tvals[m - 1] + Teq
        # misc.
        push!(
            log_info, (proposed=proposed, performed=performed, ratio=performed / proposed)
        )
        next!(progress; showvalues=[("steps", m + 1), ("total", M)])
    end
    @info "Sampling done in $time seconds"

    sequences = fmt_output(S, alignment_output, translate_output; names=tvals)
    return (; sequences, tvals, info=log_info, params=return_params(params, s0))
end

#=================================#
########### CONTINUOUS ############
#=================================#

function mcmc_sample_continuous_chain(
    g::PottsGraph,
    M::Integer,
    s0::AbstractSequence,
    params::SamplingParameters{<:AbstractFloat};
    rng=Random.GLOBAL_RNG,
    progress_meter=true,
    alignment_output=true,
    translate_output=false,
)
    # Argument checks
    @argcheck params.sampling_type == :continuous
    @argcheck M > 0 "Number of samples `M` must be >0. Instead $M"
    tmp_check_alphabet_consistency(g, s0)
    if !alignment_output && translate_output
        error("I have to implement this case")
    end

    @unpack Teq, burnin, step_type = params

    @info """
    Sampling $M sequences using the following settings:
        - Type of sampling: $(params.sampling_type)
        - Type of sequence = $(typeof(s0))
        - Time between samples = $(FloatType(Teq))
        - burnin = $(burnin)
        - Type of sampling = continuous
        - Step style = $(params.step_type)
        - Step meaning = $(params.step_meaning)
        - fraction of gap steps (if codon) = $(params.fraction_gap_step)
    """

    @info "Initial sequence: $s0"

    # vector for return value
    conf = copy(s0)
    S = similar([conf], M) # the sample
    tvals = Vector{FloatType}(undef, M) # number of steps at each sample
    # substitution_times = Vector{FloatType}(undef, 0) # time of each substitution

    R = if !isnothing(params.substitution_rate)
        @info "Using provided average model substitution rate $(params.substitution_rate)"
        params.substitution_rate
    else
        @info "Computing average substitution rate for the model using discrete sampling..."
        (;value, time) = @timed average_transition_rate(g, step_type, s0; rng, progress_meter)
        @info "Done in $time seconds"
        @info "Average substitution rate: $value"
        value
    end

    # State structure for the sampling
    state = CTMCState(conf)
    state.R = R

    # Burnin
    @info "Initializing with burnin time $(burnin)... "
    burnin > 0 && mcmc_steps!(state, g, Float64(params.burnin), step_type; rng) # Float to ensure this goes to the continuous version - improve it later
    S[1] = copy(state.seq)
    tvals[1] = burnin

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
    time = @elapsed for m in 2:M
        # doing Teq steps on the current configuration
        _, number_substitutions = mcmc_steps!(state, g, Float64(Teq), step_type; rng)
        # storing the result in S
        S[m] = copy(state.seq)
        tvals[m] = tvals[m - 1] + Float64(Teq)
        # misc.
        push!(
            log_info, (;number_substitutions)
        )
        next!(progress; showvalues=[("steps", m + 1), ("total", M)])
    end
    @info "Sampling done in $time seconds"

    sequences = fmt_output(S, alignment_output, translate_output; names=tvals)
    return (; sequences, tvals, info=log_info, params=return_params(params, s0))
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

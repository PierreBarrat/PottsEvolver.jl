@kwdef mutable struct CTMCState{S<:AbstractSequence}
    seq::S # current working sequence
    previous_seq::Union{Nothing, S} # used for efficient calculation of ΔE
    ΔE::Matrix{FloatType} # between seq and all other sequences at one mutation
    accessibility_mask::Matrix{Bool} = similar(ΔE, Bool)
    qL_buffer::Matrix{FloatType} = similar(ΔE) # useful for i.e. transition rates
    i::Union{Nothing,Integer} # position of the next mutation
    x::Union{Nothing,Integer} # state of the next mutation
end

function CTMCState(sequence::CodonSequence)
    q = length(codon_alphabet)
    L = length(sequence)
    return CTMCState{CodonSequence}(
        sequence,
        nothing,
        Matrix{FloatType}(undef, q, L),
        Matrix{FloatType}(undef, q, L),
        Matrix{FloatType}(undef, q, L),
        nothing,
        nothing,
    )
end
function CTMCState(sequence::AASequence)
    q = 21
    L = length(sequence)
    return CTMCState{NumSequence{q}}(
        sequence,
        nothing,
        Matrix{FloatType}(undef, q, L),
        Matrix{FloatType}(undef, q, L),
        Matrix{FloatType}(undef, q, L),
        nothing,
        nothing,
    )
end
function CTMCState(sequence::NumSequence{q}) where {q}
    L = length(sequence)
    return CTMCState{NumSequence{q}}(
        sequence,
        nothing,
        Matrix{FloatType}(undef, q, L),
        Matrix{FloatType}(undef, q, L),
        Matrix{FloatType}(undef, q, L),
        nothing,
        nothing,
    )
end

#=============================================#
################ Main function ################
#=============================================#

function mcmc_sample_continuous_chain(
    g::PottsGraph,
    M::Integer,
    s0::AbstractSequence,
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
    progress_meter=true,
    alignment_output=true,
    translate_output=false,
)
    # Argument checks
    if !alignment_output && translate_output
        error("I have to implement this case")
    end
    @unpack Teq, burnin, step_type = params
    @argcheck Teq > 0 "Time between samples `Teq` should be >0. Instead $(Teq)"
    @argcheck M > 0 "Number of samples `M` must be >0. Instead $M"
    tmp_check_alphabet_consistency(g, s0)

    @info """
    Sampling $M sequences using the following settings:
        - Type of sequence = $(supertype(typeof(s0)))
        - Time between samples = $(Teq)
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

    # State structure for the sampling
    state = CTMCState(conf)

    # Burnin
    @info "Initializing with burnin time $(burnin)... "
    burnin > 0 && mcmc_steps!(state, g, Teq, step_type; rng)
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
        substitution_times = mcmc_steps!(state, g, Teq, step_type; rng)
        # storing the result in S
        S[m] = copy(state.seq)
        tvals[m] = tvals[m - 1] + Teq
        # misc.
        push!(
            log_info, (;substitution_times)
        )
        next!(progress; showvalues=[("steps", m + 1), ("total", M)])
    end
    @info "Sampling done in $time seconds"

    sequences = fmt_output(S, alignment_output, translate_output; names=tvals)
    return (; sequences, tvals, info=log_info, params=return_params(params, s0))
end


#=================================#
############ Gillespie ############
#=================================#

function mcmc_steps!(
    sequence::AbstractSequence, g, Tmax::AbstractFloat, step_type; kwargs...
)
    state = CTMCState(sequence)
    return mcmc_steps!(state, g, Tmax, step_type; kwargs...)
end

"""
    mcmc_steps!(
        s::AbstractSequence, g, Tmax::Float64, p::SamplingParameters;
        gibbs_holder, kwargs...
    )

Sample `g` during time `Tmax` from sequence `s`.
The step type (`:glauber`, `:metropolis`) is set using `p` (see `?SamplingParameters`).
Modifies the input sequence `s` and returns it.
"""
function mcmc_steps!(
    state::CTMCState,
    g::PottsGraph,
    Tmax::AbstractFloat,
    step_type::Symbol;
    kwargs...
)
    @argcheck length(state.seq) == size(g).L """
    Size of sequence and model do not match: $(length(state.seq)) and $(size(g).L)
    """

    substitution_times = gillespie!(state, g, Tmax, step_type; kwargs...)
    return state.seq, substitution_times # all substitutions are accepted
end

function gillespie!(state, g, Tmax, step_type; rng=Random.GLOBAL_RNG)
    times = Float64[]

    t = 0
    while true
        Q, R = transition_rates!(state, g, step_type) # R ~ sum(Q): total rate
        @assert R > 0 "Something went wrong in gillespie: null total rate"
        # pick time of next substitution
        Δt = rand(rng, Exponential(1/R))
        t += Δt
        t > Tmax && break

        push!(times, t)

        # pick next substitution
        r = rand(rng)*R
        h = 0
        a = 0
        i = 0
        for ii in 1:L, aa in 1:q
            h += Q[aa,ii]
            if h > r
                a = aa
                i = ii
                break
            end
        end
        @assert a > 0 && i > 0 "Something went wrong in gillespie: no substitution chosen"

        # do the substitution
        gillespie_substitute!(state, i, a)
    end

    return times
end

function gillespie_substitute!(state, i, a)
    # copy sequence to state.previous_sequence, without allocation
    copy!(state.previous_seq, state.seq)
    # change sequence
    _update_sequence!(state.seq, i, a)
    # set i, a
    state.i = i
    state.x = a

    return state
end
function _update_sequence!(sequence::AbstractSequence, i, a)
    sequence.seq[i] = a
end
function _update_sequence!(sequence::CodonSequence, i, a)
    sequence.seq[i] = a
    sequence.aaseq[i] = genetic_code(a)
end

#========================================================================#
############################ Transition rates ############################
#========================================================================#

function transition_rates!(state::CTMCState, g::PottsGraph, step_type)
    # Compute energy differences
    state.ΔE .= if isnothing(state.previous_seq)
        compute_energy_differences(state.seq, g) # do the calculation from scratch
    else
        compute_energy_differences(state, g) # use the previous sequence and ΔE to help
    end

    # filter out inaccessible sequences (because of genetic code)
    set_accessibility_mask!(state.accessibility_mask, state.seq)


    return if step_type == :glauber
        transition_rates_glauber!(state, g)
    end
end

function transition_rates_glauber!(state::CTMCState, g::PottsGraph)
    # Assume that `state` has `ΔE` already computed and filtered
    (;L, q) = size(g)
    Q = state.qL_buffer
    for i in 1:L, a in 1:q
        Q[a, i] = state.accessibility_mask[a,i] ? 1/(1 + exp(state.ΔE[a, i])) : 0
    end
    Qtot = sum(Q)
    return Q, Qtot
end

"""
    set_accessibility_mask!(mask::Matrix{Bool}, sequence)

`mask` is a `q x L` matrix. For each position `i`, set `mask[c, i]` to `true` if and only if state `c` is accessible from `seq[i]`.
For non codon sequences, `mask` is trivially all true.
For `CodonSequence`, many entries are false.
"""
function set_accessibility_mask!(mask::Matrix{Bool}, seq::CodonSequence)
    for i in 1:length(seq)
        mask[:, i] .= false
        for a in accessible_codons[seq[i]][1]
            mask[a, i] = true
        end
    end
    return mask
end
function set_accessibility_mask!(mask::Matrix{Bool}, seq)
    mask .= false
    return mask
end

#============================================================#
##################### Energy differences #####################
#============================================================#


function compute_energy_differences(state::CTMCState, g::PottsGraph)
    return compute_energy_differences(
        state.ΔE, state.previous_seq, state.i, state.x, g; ΔE=state.qL_buffer
    )
end

## General case
"""
    compute_energy_differences(
        ΔE::Vector{<:AbstractFloat},
        refseq::AbstractSequence,
        i::Integer,
        x::Integer,
        g::PottsGraph;
        rng=Random.GLOBAL_RNG,
        ΔE,
    )
    compute_energy_differences(refseq::AbstractSequence, g::PottsGraph)

In the first form: define the new reference `b` sequence is equal to `refseq` at positions different from `i`, and to `x` at position `i`.
For all mutations `(j,z)`, giving rise to sequence `c`, compute `E(c)-E(b)` and store the result in output matrix `ΔE` of dimensions `qxL`.
Keyword argument `ΔE` can be provided to avoid allocation.


In the second form: compute the energy difference between any sequence `b` at distance one from `a` and `a` itself.
The output is a matrix of dimensions `qxL`, always allocated in the function.

*Note*: For `CodonSequence`, the output is of size `65xL` where `L` is the length.
Stop lead to infinite energy. Providing a reference sequence that contains a stop codon is undefined.
"""
function compute_energy_differences(
    ref_ΔE::Matrix{FloatType},
    refseq::AbstractSequence,
    i::Integer,
    x::Integer,
    g::PottsGraph;
    ΔE=similar(ref_ΔE),
)
    q, L = size(ref_ΔE)
    ai = sequence(refseq)[i]
    for j in 1:L
        if j == i
            for y in 1:q
                ΔE[y,i] = ref_ΔE[y,i] - ref_ΔE[x,i]
            end
        else
            aj = sequence(refseq)[j]
            for y in 1:q
                ΔE[y,j] = ref_ΔE[y,j]
                ΔE[y,j] += g.J[aj, x, j, i] - g.J[aj, ai, j, i]
                ΔE[y,j] -= g.J[y, x, j, i] - g.J[y, ai, j, i]

            end
        end
    end
    return ΔE
end

function compute_energy_differences(refseq::AASequence, g::PottsGraph)
    ΔE = zeros(FloatType, 21, length(refseq))
    _compute_energy_differences!(ΔE, refseq, g)
    return ΔE
end
function compute_energy_differences(refseq::NumSequence{T,q}, g::PottsGraph) where {T,q}
    ΔE = zeros(FloatType, q, length(refseq))
    _compute_energy_differences!(ΔE, refseq, g)
    return ΔE
end
function _compute_energy_differences!(
    ΔE_buffer::Matrix{FloatType}, refseq::AbstractSequence, g::PottsGraph;
)
    q, L = size(ΔE_buffer)
    seq = copy(refseq)
    for i in 1:L, a in 1:q
        seq.seq[i] = a
        ΔE_buffer[a, i] = energy(seq, g) - energy(refseq, g) # can be improved
        seq.seq[i] = refseq.seq[i]
    end
    return ΔE_buffer
end



## Codon case
function compute_energy_differences(
    ref_ΔE::Matrix{FloatType},
    refseq::CodonSequence,
    i::Integer,
    x::Integer,
    g::PottsGraph;
    ΔE=similar(ref_ΔE),
)
    L = length(refseq)
    q = 21
    ai = refseq.aaseq[i]
    x_aa = genetic_code(x)
    for j in 1:L
        if j == i
            for y_aa in 1:q
                codons = reverse_code(y_aa)
                new_ΔE = ref_ΔE[first(codons),i] - ref_ΔE[x,i]
                for y in codons
                    ΔE[y,i] = new_ΔE
                end
            end
        else
            aj = refseq.aaseq[j]
            for y_aa in 1:q
                codons = reverse_code(y_aa)
                y = first(codons)
                new_ΔE = ref_ΔE[y,j]
                new_ΔE += g.J[aj, x_aa, j, i] - g.J[aj, ai, j, i]
                new_ΔE -= g.J[y_aa, x_aa, j, i] - g.J[y_aa, ai, j, i]
                for y in codons
                    ΔE[y,j] = new_ΔE
                end
            end
        end
    end
    return ΔE
end

function compute_energy_differences(refseq::CodonSequence, g::PottsGraph)
    q = length(codon_alphabet)
    L = length(refseq)
    ΔE = zeros(Float64, q, L)

    seq = copy(refseq)
    for i in 1:L, c in 1:q
        if isstop(c)
            ΔE[c, i] = Inf
        else
            seq.aaseq[i] = genetic_code(c)
            ΔE[c, i] = energy(seq, g) - energy(refseq, g)
            seq.aaseq[i] = refseq.aaseq[i]
        end
    end

    return ΔE
end

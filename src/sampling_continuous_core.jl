@kwdef mutable struct CTMCState{S<:AbstractSequence}
    seq::S # current working sequence
    previous_seq::Union{Nothing, S} # used for efficient calculation of ΔE
    ΔE::Matrix{FloatType} # between seq and all other sequences at one mutation
    accessibility_mask::Matrix{Bool} = similar(ΔE, Bool)
    qL_buffer::Matrix{FloatType} = similar(ΔE) # useful for i.e. transition rates
    ΔE_copy::Matrix{FloatType} = copy(ΔE)
    i::Union{Nothing,Integer} # position of the next mutation
    x::Union{Nothing,Integer} # state of the next mutation
    R::Union{Nothing,FloatType} = nothing # average transition rate, used for scaling
    function CTMCState{S}(sequence::S, q::Integer, L::Integer) where {S<:AbstractSequence}
        return new{S}(
            sequence,
            nothing,
            Matrix{FloatType}(undef, q, L), # ΔE
            Matrix{Bool}(undef, q, L), # accessibility mask
            Matrix{FloatType}(undef, q, L), # qxL buffer
            Matrix{FloatType}(undef, q, L), # ΔE copy
            nothing,
            nothing,
            nothing,
        )
    end
end

function CTMCState(sequence::CodonSequence)
    q = length(codon_alphabet)
    L = length(sequence)
    return CTMCState{CodonSequence}(sequence, q, L)
end
function CTMCState(sequence::AASequence)
    q = 21
    L = length(sequence)
    CTMCState{AASequence}(sequence, q, L)
end
function CTMCState(sequence::NumSequence{T,q}) where {T,q}
    L = length(sequence)
    CTMCState{NumSequence{T,q}}(sequence, q, L)
end

Base.size(state::CTMCState) = size(state.ΔE)

#=============================================#
################ Main function ################
#=============================================#

"""
    mcmc_steps!(
        sequence::AbstractSequence, g, Tmax::AbstractFloat, parameters::SamplingParameters; 
        kwargs...
    )
    mcmc_steps!(
        state::CTMCState, g, Tmax::AbstractFloat, parameters::SamplingParameters; kwargs...
    )

Sample `g` during time `Tmax` (continuous) from sequence `seq`.
The step type is `p.step_type` (see `?SamplingParameters`).
Modifies the input sequence `s` and returns it.
Allocates a `CTMCState` (three matrices of order `q*L`).

## Notes
- The form with `state::CTMCState` will use a pre-allocated `CTMCState` for calculations. 
- Expects `parameters.sampling_type` to be `:continuous`. Fails if otherwise. 
"""
function mcmc_steps!(
    sequence::AbstractSequence, g, Tmax::AbstractFloat, parameters::SamplingParameters; 
    kwargs...
)
    @argcheck parameters.sampling_type == :continuous """
    `Tmax` is `AbstractFloat`: expected sampling type to be :continuous. 
    Instead :$(parameters.sampling_type).
    """
    state = CTMCState(sequence)
    return mcmc_steps!(state, g, Tmax, parameters.step_type; kwargs...)
end
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
    n_substitutions = gillespie!(state, g, Tmax, step_type; kwargs...)
    return state.seq, n_substitutions # all substitutions are accepted
end


#=================================#
############ Gillespie ############
#=================================#

function gillespie!(state, g, Tmax, step_type; rng=Random.GLOBAL_RNG)
    Tmax == 0 && return Float64[]

    q, L = size(state)
    # times = Float64[]
    # sizehint!(times, Tmax*L)
    n_substitutions = 0

    t = 0
    while true
        state.ΔE_copy .= state.ΔE # memorize the old ΔE in case no move is made

        # compute transition rates and scaled substitution rate
        Q = transition_rates!(state, g, step_type)
        R = sum(Q)
        @assert R > 0 "Something went wrong in gillespie: null total rate"
        R_scaled = isnothing(state.R) ? R : R/state.R*L

        # pick time of next substitution
        Δt = rand(rng, Exponential(1/R_scaled))
        t += Δt

        # If time is too large, we stop
        # This does not really influence the distribution of time to the next substitution
        # as an exponential process is memory-less
        if t > Tmax
            state.ΔE .= state.ΔE_copy # reset `state` to its old state since no move is made
            # state.seq and other important fields have not been modified yet
            # state.qL_buffer and state.accessibility_mask were changed, but this is ok
            break
        end

        # Time is in range (0, Tmax], we choose a substitution
        # push!(times, t)

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
        n_substitutions += 1
    end

    # return times
    return n_substitutions
end

function gillespie_substitute!(state, i, a)
    # copy sequence to state.previous_sequence, without allocation
    if isnothing(state.previous_seq)
        state.previous_seq = copy(state.seq)
    else
        copy!(state.previous_seq, state.seq)
        # state.previous_seq = copy(state.seq)
    end
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


#==============================================================================#
########################### Average transition rates ###########################
#==============================================================================#

"""
    average_transition_rate(g::PottsGraph, step_type, s0::AbstractSequence; kwargs...)

Compute the average transition rate for the continuous time Markov chain based on `g`. 
`s0` is only used to initialize the `CTMCState` and provide a type. 
"""
function average_transition_rate(
    g::PottsGraph, step_type, s0::AbstractSequence
    ; rng=Random.GLOBAL_RNG, progress_meter=true
)
    (;L) = size(g)

    M = 100
    params = SamplingParameters(Teq = L*10) # default continuous sampling parameters
    sample_discrete = mcmc_sample(
        g, M, params; 
        alignment_output=false, init=s0, rng=rng, progress_meter=progress_meter
    ).sequences
    state = CTMCState(sample_discrete[1])
    Rmean = mean(sample_discrete) do sequence
        copy!(state.seq, sequence)
        sum(transition_rates!(state, g, step_type))
    end

    return Rmean
end

#========================================================================#
############################ Transition rates ############################
#========================================================================#

"""
    transition_rates(sequence::AbstractSequence, g::PottsGraph, step_type)

Compute the transition rate matrix for the sequence `sequence`.
Return the a `q` by `L` matrix `Q` and the total rate `R`.
Allocates a `CTMCState`. 
"""
function transition_rates(sequence::AbstractSequence, g::PottsGraph, step_type)
    state = CTMCState(sequence)
    return transition_rates!(state, g, step_type)
end

"""
    transition_rates!(state::CTMCState, g::PottsGraph, step_type; from_scratch=false)

Compute the transition rate matrix for the sequence `state.seq`.
Return the a `q` by `L` matrix `Q` and the total rate `R`. 
"""
function transition_rates!(state::CTMCState, g::PottsGraph, step_type; from_scratch=false)
    # Compute energy differences
    state.ΔE .= compute_energy_differences!(state, g; from_scratch)

    # filter out inaccessible sequences (because of genetic code)
    set_accessibility_mask!(state.accessibility_mask, state.seq)

    return if step_type == :glauber
        transition_rates_glauber!(state, g)
    elseif step_type == :metropolis
        transition_rates_metropolis!(state, g)
    elseif step_type == :sqrt
        transition_rates_sqrt!(state, g)
    else
        throw(ArgumentError("Step type $step_type not implemented"))
    end
end

function transition_rates_glauber!(state::CTMCState, g::PottsGraph)
    # Assume that `state` has `ΔE` already computed and filtered
    (q, L) = size(state)
    @inbounds for i in 1:L, a in 1:q
        if state.accessibility_mask[a,i]
            state.qL_buffer[a, i] = 1/(1. + exp(state.ΔE[a, i]))
        else
            state.qL_buffer[a, i] = 0.
        end
    end
    return state.qL_buffer
end

function transition_rates_metropolis!(state::CTMCState, g::PottsGraph)
    # Assume that `state` has `ΔE` already computed and filtered
    (q, L) = size(state)
    @inbounds for i in 1:L, a in 1:q
        if state.accessibility_mask[a,i]
            if state.ΔE[a, i] < 0
                state.qL_buffer[a, i] = 1.
            else
                state.qL_buffer[a, i] = exp(-state.ΔE[a, i])
            end
        else
            state.qL_buffer[a, i] = 0.
        end
    end
    return state.qL_buffer
end

function transition_rates_sqrt!(state::CTMCState, g::PottsGraph)
    # Assume that `state` has `ΔE` already computed and filtered
    (q, L) = size(state)
    @inbounds for i in 1:L, a in 1:q
        if state.accessibility_mask[a,i]
            state.qL_buffer[a, i] = exp(-0.5*state.ΔE[a, i])
        else
            state.qL_buffer[a, i] = 0
        end
    end
    return state.qL_buffer
end


"""
    set_accessibility_mask!(mask::Matrix{Bool}, sequence)

Set the accessibility mask for a sequence.

# Arguments
- `mask::Matrix{Bool}`: The mask to set.
- `sequence`: The sequence to determine accessibility.

# Returns
The updated mask.
"""
function set_accessibility_mask!(mask::Matrix{Bool}, seq::CodonSequence)
    for i in 1:length(seq)
        mask[:, i] .= false
        for a in accessible_codons(seq[i])[1]
            mask[a, i] = true
        end
    end
    return mask
end

"""
    set_accessibility_mask!(mask::Matrix{Bool}, sequence)

Set the accessibility mask for a sequence.

# Arguments
- `mask::Matrix{Bool}`: The mask to set.
- `sequence`: The sequence to determine accessibility.

# Returns
The updated mask.
"""
function set_accessibility_mask!(mask::Matrix{Bool}, seq::AbstractSequence)
    q, L = size(mask)
    for i in 1:L, a in 1:q
        mask[a, i] = !(sequence(seq)[i] == a) # cannot mutate to current state
    end
    return mask
end

#============================================================#
##################### Energy differences #####################
#============================================================#




# dispatch to two cases
function compute_energy_differences!(state::CTMCState, g::PottsGraph; from_scratch=false)
    # the buffer of state is changed
    return if isnothing(state.previous_seq) || from_scratch
        # do the calculation from scratch, using the buffer to avoid allocation
        compute_energy_differences!(state.qL_buffer, state.seq, g) 
    else
        compute_energy_differences!(
            state.qL_buffer, state.ΔE, state.previous_seq, state.i, state.x, g,
        ) 
    end    
end

# Case with using the previous sequence ~ update
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
function compute_energy_differences!(
    ΔE::Matrix{FloatType},
    ref_ΔE::Matrix{FloatType},
    refseq::AbstractSequence,
    i::Integer,
    x::Integer,
    g::PottsGraph,
)
    q, L = size(ref_ΔE)
    ai = sequence(refseq)[i]
    ΔE .= ref_ΔE
    @inbounds for j in 1:L
        if j == i
            for y in 1:q
                ΔE[y,i] -= ref_ΔE[x,i]
            end
        else
            aj = sequence(refseq)[j]
            dJ = g.J[aj, x, j, i] - g.J[aj, ai, j, i]
            for y in 1:q
                ΔE[y,j] += dJ
                ΔE[y,j] -= g.J[y, x, j, i] - g.J[y, ai, j, i]
            end
        end
    end
    return ΔE
end
function compute_energy_differences!(
    ΔE::Matrix{FloatType},
    ref_ΔE::Matrix{FloatType},
    refseq::CodonSequence,
    i::Integer,
    x::Integer, # codon
    g::PottsGraph,
)
    L = length(refseq)
    q = IntType(21)
    ai = refseq.aaseq[i] # previous aa at pos. i
    x_aa = genetic_code(x) # current aa at pos. i
    ΔE .= ref_ΔE
    for j in 1:L
        # The calculation below considers all amino acids y_aa
        # This is inefficient, as some amino-acids are not accessible from the current codon x
        # However, it allows me to keep an updated ΔE matrix
        # If considering only accessible codons, I believe one would have to eventually
        # recompute ΔE from scratch at each step, without being able to use the ΔE from
        # the previous sequence
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
            dJ = g.J[aj, x_aa, j, i] - g.J[aj, ai, j, i]
            for y_aa in 1:q
                codons = reverse_code(y_aa)
                y = first(codons)
                new_ΔE = 0.
                new_ΔE += dJ
                new_ΔE -= g.J[y_aa, x_aa, j, i] - g.J[y_aa, ai, j, i]
                for y in codons
                    ΔE[y,j] += new_ΔE
                end
            end
        end
    end
    return ΔE
end

# Case from scratch, without using a pre-existing sequence/ΔE pairs
"""
    compute_energy_differences(refseq::AbstractSequence, g::PottsGraph)

Compute the energy differences between `refseq` and neighbouring sequences.
Return a matrix of dimensions `q` by `L`.
"""
function compute_energy_differences(refseq::AASequence, g::PottsGraph)
    ΔE = zeros(FloatType, 21, length(refseq))
    return compute_energy_differences!(ΔE, refseq, g)
end
function compute_energy_differences(refseq::NumSequence{T,q}, g::PottsGraph) where {T,q}
    ΔE = zeros(FloatType, q, length(refseq))
    return compute_energy_differences!(ΔE, refseq, g)
end
function compute_energy_differences(refseq::CodonSequence, g::PottsGraph)
    q = length(codon_alphabet)
    L = length(refseq)
    ΔE = zeros(Float64, q, L)
    return compute_energy_differences!(ΔE, refseq, g)
end

"""
    compute_energy_differences!(
        ΔE::Matrix{FloatType}, refseq::AbstractSequence, g::PottsGraph;
    )

Compute the energy differences between `refseq` and neighbouring sequences, storing result 
in pre-allocated matrix ΔE. 
"""
function compute_energy_differences!(
    ΔE::Matrix{FloatType}, refseq::AbstractSequence, g::PottsGraph;
)
    q, L = size(ΔE)
    seq = copy(refseq)
    for i in 1:L, a in 1:q
        seq.seq[i] = a
        ΔE[a, i] = _delta_energy(seq, refseq, i, g)
        seq.seq[i] = refseq.seq[i]
    end
    return ΔE
end
function compute_energy_differences!(
    ΔE::Matrix{FloatType}, refseq::CodonSequence, g::PottsGraph;
)
    (q, L) = size(ΔE)
    seq = copy(refseq)
    for i in 1:L, c in 1:q
        if isstop(c)
            ΔE[c, i] = Inf
        else
            seq.aaseq[i] = genetic_code(c)
            ΔE[c, i] = _delta_energy(seq, refseq, i, g)
            ΔE[c, i] += sum(aa_degeneracy, seq.aaseq) - sum(aa_degeneracy, refseq.aaseq)
        end
        seq.aaseq[i] = refseq.aaseq[i]
    end

    return ΔE
end

function _delta_energy(seq::CodonSequence, refseq::CodonSequence, i, g)
# energy difference between seq and refseq, assuming that they differ only at i
    dE = -g.h[seq.aaseq[i], i] + g.h[refseq.aaseq[i], i]
    for j in 1:length(seq)
        if j != i
            @argcheck seq.seq[j] == refseq.seq[j]
            dE += -g.J[seq.aaseq[i], seq.aaseq[j], i, j]
            dE += g.J[refseq.aaseq[i], seq.aaseq[j], i, j]
        end
    end
    return dE
end
function _delta_energy(seq::AbstractSequence, refseq::AbstractSequence, i, g)
    # energy difference between seq and refseq, assuming that they differ only at i
    dE = -g.h[seq[i], i] + g.h[refseq[i], i]
    for j in 1:length(seq)
        if j != i
            @argcheck seq[j] == refseq[j]
            dE += -g.J[seq[i], seq[j], i, j] + g.J[refseq[i], seq[j], i, j]
        end
    end
    return dE
end



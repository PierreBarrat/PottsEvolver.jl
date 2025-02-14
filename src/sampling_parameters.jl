#============================================================#
##################### SamplingParameters #####################
#============================================================#

const VALID_STEP_TYPES = [:gibbs, :metropolis, :glauber, :sqrt]
const VALID_STEP_MEANINGS = [:proposed, :accepted, :changed]
const VALID_STEP_TYPES_CONTINUOUS = [:metropolis, :glauber, :sqrt]
const VALID_STEP_TYPES_DISCRETE = [:gibbs, :metropolis]

"""
```
BranchLengthMeaning

    type::Symbol
    length::Symbol
```

Define how the branch lengths should be interpreted when sampling on a tree.
This is done in two consecutive steps.
1. Convert the branch lengths to mcmc steps.
    If `type==:sweep`, multiply the branch length by the length of the sequence.
    Else, if `type==:step`, leave it as is.
2. Convert the branch length to an integer value which will be the number of mcmc steps.
  - if `length==:exact`, the branch length is expected to be an integer (`rtol=1e-6`).
  - if `length==:round`, round it to the closest `Int`.
  - if `length==:poisson`, sample a poisson distribution with the branch length as mean.
"""
@kwdef struct BranchLengthMeaning
    type::Symbol
    length::Symbol
    function BranchLengthMeaning(type, length)
        valid_types = (:step, :sweep)
        valid_lengths = (:exact, :round, :poisson)
        @argcheck type in valid_types """
        `type` field should be in `$(valid_types)`. Instead `:$type`.
        """
        @argcheck length in valid_lengths """
        `length` field should be in `$(valid_lengths)`. Instead `:$length`.
        """
        return new(type, length)
    end
end

function steps_from_branchlength(τ::Real, p::BranchLengthMeaning, L::Int)
    if p.type == :sweep
        τ *= L
    end

    rtol = 1e-6
    if p.length == :exact
        τi = round(Int, τ)
        if !isapprox(τ, τi; rtol)
            msg = """
            Expected tree with approx integer branch lengths (`rtol=$rtol`). Instead $τ.
            """
            throw(ArgumentError(msg))
        end
        return τi
    elseif p.length == :round
        return round(Int, τ)
    elseif p.length == :poisson
        return pois_rand(τ)
    else throw(ArgumentError("Invalid `length` field $(p.length)"))
    end
end

"""
    mutable struct SamplingParameters

Construct using keyword arguments:
```
sampling_type::Symbol = :discrete
step_type::Symbol = :gibbs
step_meaning::Symbol = :accepted
Teq::Int
burnin::Int = 5*Teq
fraction_gap_step::Float64 = 0.9
branchlength_meaning::BranchLengthMeaning
```

- `Teq` is measured in swaps: attempted (or accepted) change of one sequence position.
- `burnin`: number of steps starting from the initial sequence before the first `Teq` are
    made.
- `step_meaning` is only relevant for discrete sampling. It can take three values
    - `:proposed`: all mcmc steps count towards equilibration
    - `:accepted`: only accepted steps count (all steps for non-codon Gibbs)
    - `:changed`: only steps that lead to a change count (Gibbs can resample the same state)
    *Note*: Gibbs steps for codons are more complicated, since they involve the possibility
    Metropolis step for gaps, which can be rejected.

- `fraction_gap_step` is only relevant for discrete sampling with codons.
- `branchlength_meaning` is only useful if you sample along branches of a tree.
  See `?BranchLengthMeaning` for information.
"""
@kwdef mutable struct SamplingParameters{T<:Real}
    sampling_type::Symbol = :discrete
    step_type::Symbol = :gibbs
    Teq::T = Int(0)
    burnin::T = 5 * Teq
    # for discrete sampling only
    step_meaning::Symbol = :accepted
    # for continuous sampling only - average substitution rate for a given Potts model
    substitution_rate::Union{Nothing,FloatType} = nothing
    # for the discrete codon algorithm (Gibbs for aa and Metropolis for gaps)
    fraction_gap_step::Float64 = 0.9
    # when sampling on trees
    branchlength_meaning::BranchLengthMeaning = BranchLengthMeaning(:step, :exact)
    function SamplingParameters(
        sampling_type,
        step_type,
        Teq::T,
        burnin::T,
        step_meaning,
        substitution_rate,
        fraction_gap_step,
        branchlength_meaning,
    ) where {T}
        step_meaning = try
            Symbol(step_meaning)
        catch err
            @error "Invalid `step_meaning` $step_meaning"
            throw(err)
        end

        @argcheck sampling_type in (:discrete, :continuous)
        if sampling_type == :discrete
            @argcheck step_type in VALID_STEP_TYPES_DISCRETE """
            For :discrete mcmc, `step_type` should be in $VALID_STEP_TYPES_DISCRETE.
            Instead $(step_type).
            """
            if !(T <: Integer)
                try 
                    Teq = Int(Teq)
                    burnin = Int(burnin)
                catch err
                    @error """
                    For :discrete mcmc, `Teq` should be an integer. Instead $(Teq).
                    """
                    throw(err)
                end
            end
        elseif sampling_type == :continuous
            @argcheck step_type in VALID_STEP_TYPES_CONTINUOUS """
            For :continuous mcmc, `step_type` should be in $VALID_STEP_TYPES_CONTINUOUS.
            Instead $(step_type).
            """
        end

        @argcheck step_meaning in VALID_STEP_MEANINGS """
            `step_meaning` should be in $VALID_STEP_MEANINGS.
            Instead $(step_meaning).
        """

        @argcheck isnothing(substitution_rate) || substitution_rate > 0

        return new{T}(
            sampling_type,
            step_type,
            Teq,
            burnin,
            step_meaning,
            substitution_rate,
            fraction_gap_step,
            branchlength_meaning,
        )
    end
end

function Base.convert(::Type{SamplingParameters{T}}, params::SamplingParameters) where {T}
    Teq = T(params.Teq)
    burnin = T(params.burnin)
    return SamplingParameters(
        params.sampling_type,
        params.step_type,
        Teq,
        burnin,
        params.step_meaning,
        params.substitution_rate,
        params.fraction_gap_step,
        params.branchlength_meaning,
    )
end
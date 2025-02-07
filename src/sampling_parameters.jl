#============================================================#
##################### SamplingParameters #####################
#============================================================#

const VALID_STEP_TYPES = [:gibbs, :metropolis]
const VALID_STEP_MEANINGS = [:proposed, :accepted, :changed]

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
    end
    @assert false "This state should never be reached"
end

"""
    mutable struct SamplingParameters

Construct using keyword arguments:
```
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
- `step_meaning` can take three values
    - `:proposed`: all mcmc steps count towards equilibration
    - `:accepted`: only accepted steps count (all steps for non-codon Gibbs)
    - `:changed`: only steps that lead to a change count (Gibbs can resample the same state)
    *Note*: Gibbs steps for codons are more complicated, since they involve the possibility
    Metropolis step for gaps, which can be rejected.
- `branchlength_meaning` is only useful if you sample along branches of a tree.
  See `?BranchLengthMeaning` for information.
"""
@kwdef mutable struct SamplingParameters{T<:Real}
    step_type::Symbol = :gibbs
    step_meaning::Symbol = :accepted
    Teq::T
    burnin::T = 5 * Teq
    fraction_gap_step::Float64 = 0.9
    branchlength_meaning::BranchLengthMeaning = BranchLengthMeaning(:step, :exact)
    function SamplingParameters(
        step_type, step_meaning, Teq::T, burnin::T, fraction_gap_step, branchlength_meaning
    ) where {T}
        step_meaning = try
            Symbol(step_meaning)
        catch err
            @error "Invalid `step_meaning` $step_meaning"
            throw(err)
        end
        @argcheck step_meaning in VALID_STEP_MEANINGS """
                `step_meaning` should be in $VALID_STEP_MEANINGS.
                Instead $(step_meaning).
            """
        @argcheck step_type in VALID_STEP_TYPES """
                `step_type` should be in $VALID_STEP_TYPES.
                Instead $(step_type).
            """
        return new{T}(
            step_type, step_meaning, Teq, burnin, fraction_gap_step, branchlength_meaning
        )
    end
end

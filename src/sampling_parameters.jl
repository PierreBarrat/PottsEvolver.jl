#============================================================#
##################### SamplingParameters #####################
#============================================================#

const VALID_STEP_TYPES = [:gibbs, :metropolis, :glauber, :sqrt]
const VALID_STEP_MEANINGS = [:proposed, :accepted, :changed]
const VALID_STEP_TYPES_CONTINUOUS = [:metropolis, :glauber, :sqrt]
const VALID_STEP_TYPES_DISCRETE = [:gibbs]

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
Teq::Real = 0
burnin::Real = 5*Teq
step_meaning::Symbol = :accepted # relevant for :discrete
fraction_gap_step::Float64 = 0.9 # relevant for :discrete and codon sampling
branchlength_meaning::BranchLengthMeaning # relevant for :discrete and sampling on a tree
substitution_rate::Union{Nothing,Float64} = nothing # relevant for :continuous
```

`sampling_type` can be `:discrete` or `:continuous`.
The time between samples is `Teq`, and an extra time `burnin` is taken before the first sample. 
Allowed values of `step_type` differ between the continuous and discrete cases. 

## Discrete sampling
- `Teq` is measured in swaps: accepted (or attempted) changes of one sequence position.
  It must be an integer.
- `burnin`: number of steps starting from the initial sequence before the first sample. 
  Must be an integer. 
- `step_type`: `:gibbs` only. I intend to implement `:metropolis` and maybe `:glauber`.
- `step_meaning` (:discrete only): whether a swap count towards equilibration or not. 
  It can take three values
    - `:proposed`: all mcmc steps count towards equilibration
    - `:accepted`: only accepted steps count (all steps for non-codon Gibbs)
    - `:changed`: only steps that lead to a change count (Gibbs can resample the same state)
    *Note*: Gibbs steps for codons are more complicated, since they involve the possibility
    Metropolis step for gaps, which can be rejected.

- `fraction_gap_step`: fraction of `:metropolis` steps concerning gaps when sampling with 
  codons. This is only relevant for discrete sampling with codons.
- `branchlength_meaning`: how branch-lengths on a tree must be converted to swaps. 
    See `?BranchLengthMeaning` for information.
    This is only useful if you sample along branches of a tree.


## Continuous sampling

For continuous sampling, there is no notion of swap/sweep or of accepted/rejected steps. 
Since any positive real number is acceptable as a sampling time, branch lenghts of trees 
can be used directly. 

- `Teq` and `burnin` are floats. 
- `step_type` can be `:metropolis`, `:glauber` or `:sqrt`.
- `substitution_rate` is the average substitution rate for a given Potts model  
  (the average runs over sequences). 
  It is the result of `average_substitution_rate`. 
  This is computed automatically if not provided (but takes some time).

"""
@kwdef mutable struct SamplingParameters{T<:Real}
    sampling_type::Symbol = :discrete
    step_type::Symbol = (sampling_type == :discrete ? :gibbs : :glauber)
    Teq::T = Int(0)
    burnin::T = 5 * Teq
    # for discrete sampling only
    step_meaning::Symbol = :accepted
    # for the discrete codon algorithm (Gibbs for aa and Metropolis for gaps)
    fraction_gap_step::Float64 = 0.9
    # when sampling on trees
    branchlength_meaning::BranchLengthMeaning = BranchLengthMeaning(:step, :exact)
    # for continuous sampling only - average substitution rate for a given Potts model
    substitution_rate::Union{Nothing,Float64} = nothing
    function SamplingParameters(
        sampling_type,
        step_type,
        Teq::T,
        burnin::T,
        step_meaning,
        fraction_gap_step,
        branchlength_meaning,
        substitution_rate,
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
            fraction_gap_step,
            branchlength_meaning,
            substitution_rate,
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
        params.fraction_gap_step,
        params.branchlength_meaning,
        params.substitution_rate,
    )
end

function Base.show(io::IO, params::SamplingParameters)
    println(io, "SamplingParameters:")
    println(io, "  sampling_type = :$(params.sampling_type)")
    println(io, "  step_type = :$(params.step_type)")
    println(io, "  Teq = $(params.Teq)")
    println(io, "  burnin = $(params.burnin)")
    
    if params.sampling_type == :discrete
        println(io, "  step_meaning = :$(params.step_meaning)")
        println(io, "  fraction_gap_step = $(params.fraction_gap_step)")
        println(io, "  branchlength_meaning = BranchLengthMeaning(:$(params.branchlength_meaning.type), :$(params.branchlength_meaning.length))")
    elseif params.sampling_type == :continuous
        if !isnothing(params.substitution_rate)
            println(io, "  substitution_rate = $(params.substitution_rate)")
        else
            println(io, "  substitution_rate = not set")
        end
    end
end


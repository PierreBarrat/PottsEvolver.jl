abstract type AbstractSequence end

sequence(x::AbstractSequence; kwargs...) = x.seq

Base.getindex(s::AbstractSequence, i) = getindex(sequence(s), i)
Base.setindex!(s::AbstractSequence, x, i) = setindex!(sequence(s), x, i)
function Base.:(==)(x::T, y::T) where {T<:AbstractSequence}
    return all(p -> getproperty(x, p) == getproperty(y, p), propertynames(x))
end
function Base.hash(x::AbstractSequence, h::UInt)
    return hash(x.seq, h)
end

Base.iterate(s::AbstractSequence) = iterate(sequence(s))
Base.iterate(s::AbstractSequence, state) = iterate(sequence(s), state)
Base.length(s::AbstractSequence) = length(sequence(s))
Base.eltype(s::AbstractSequence) = eltype(sequence(s))

# used in Alignment
_sequence_alphabet(::Type{<:AbstractSequence}; kwargs...) = nothing

#=
Methods that a subtype should implement
- sequence: access integer vector (necessary if field is not called `seq`)
- copy !NECESSARY!
- equality and hash
- indexing
- _sequence_alphabet: return default alphabet for the type
=#

#====================================#
############# AASequence #############
#====================================#

"""
    mutable struct AASequence{T<:Integer} <: AbstractSequence

Field: `seq::Vector{T}`.
Wrapper around a vector of integers, with implied alphabet `PottsEvolver.aa_alphabet`.
"""
mutable struct AASequence{T<:Integer} <: AbstractSequence
    seq::Vector{T}
    function AASequence(x::AbstractVector{T}) where {T}
        q = length(aa_alphabet)
        @argcheck all(<=(q), x) "AA are represented by `(1..$(q))` integers. Instead, $x"
        return new{T}(x)
    end
end

Base.copy(s::AASequence) = AASequence(copy(s.seq))
function Base.copy!(dest::AASequence, source::AASequence)
    @argcheck length(dest) == length(source)
    for (i, a) in enumerate(source.seq)
        dest.seq[i] = a
    end
    return dest
end
"""
    AASequence(L; T)

Return a random `AASequence{T}` of length `L`.
"""
AASequence(L::Integer; T=IntType) = AASequence(rand(T(1):T(length(aa_alphabet)), L))
AASequence{T}(L::Integer) where {T<:Integer} = AASequence(L; T)

_sequence_alphabet(::Type{<:AASequence}; kwargs...) = aa_alphabet

#=============================================#
################ CodonSequence ################
#=============================================#

mutable struct CodonSequence{T<:Integer} <: AbstractSequence
    seq::Vector{T} # the codons
    aaseq::Vector{T} # the translation
    function CodonSequence(seq::Vector{T}, aaseq::Vector{T}) where {T}
        qc = length(codon_alphabet)
        qaa = length(aa_alphabet)
        @argcheck all(<=(qc), seq) """
            Codons are represented by `(1..$(qc))` integers. Instead $seq
        """
        @argcheck all(<=(qaa), aaseq) """
            AA are represented by `(1..$(qaa))` integers. Instead, $aaseq
        """
        any(isstop, seq) && @warn "Sequence contains stop codon"
        @argcheck all(x -> genetic_code(x[1]) == x[2], zip(seq, aaseq)) """
            Codon and amino acid sequences do not match. Got $seq and $aaseq
        """
        return new{T}(seq, aaseq)
    end
end

## Constructors

"""
    CodonSequence(seq::Vector{Integer}; source=:aa)

Build a `CodonSequence` from `seq`:
- if `source==:codon`, `seq` is interpreted as representing codons (see `codon_alphabet`);
- if `source==:aa`, `seq` is interpreted as representing amino acids (see `aa_alphabet`);
  matching codons are randomly chosen using the `PottsEvolver.reverse_code_rand` method.
"""
function CodonSequence(seq::AbstractVector{T}; source=:aa) where {T<:Integer}
    return if source == :aa
        CodonSequence(
            convert(Vector{T}, map(reverse_code_rand, seq)), convert(Vector{T}, seq)
        )
    elseif source == :codon
        aaseq = map(genetic_code, seq)
        any(isnothing, aaseq) && error("""
            Cannot build `CodonSequence` from input that contains stop codon.
            Input sequence was $seq.""")
        CodonSequence(convert(Vector{T}, seq), aaseq)
    else
        error("Unknown `source` $(source). Use `:aa` or `:codon`.")
    end
end
"""
    CodonSequence(L::Int; source=:aa, T)

Sample `L` states at random of the type of `source` (`:aa` or `:codon`):
- if `:codon`, sample codons at random
- if `:aa`, sample amino acids at random and reverse translate them randomly to matching codons

Underlying integer type is `T`.
"""
function CodonSequence(L::Int; source=:aa, T=IntType)
    return CodonSequence(rand(T(1):T(length(aa_alphabet)), L); source)
end
CodonSequence{T}(L::Int; kwargs...) where {T<:Integer} = CodonSequence(L; T, kwargs...)

## Methods

function Base.setindex!(s::CodonSequence, x::Integer, i)
    isstop(x) && @warn "Introducing stop codon in sequence"
    setindex!(s.aaseq, genetic_code(x), i)
    return setindex!(s.seq, x, i)
end
Base.copy(s::CodonSequence) = CodonSequence(copy(s.seq), copy(s.aaseq))
function Base.copy!(dest::CodonSequence, source::CodonSequence)
    @argcheck length(dest) == length(source)
    for (i, (c, aa)) in enumerate(zip(source.seq, source.aaseq))
        dest.seq[i] = c
        dest.aaseq[i] = aa
    end
    return dest
end

translate(s::CodonSequence) = AASequence(s.aaseq)
function sequence(x::CodonSequence; as_codons=true)
    return as_codons ? x.seq : x.aaseq
end

function _sequence_alphabet(::Type{<:CodonSequence}; as_codons=true)
    return as_codons ? codon_alphabet : aa_alphabet
end

#============================================================#
##################### Numerical sequence #####################
#============================================================#

"""
    NumSequence{T<:Integer, q}

A mutable struct representing a sequence of integers with a maximum value constraint `q`.
1. **Explicit Construction**: `NumSequence(seq::AbstractVector{T}, q::Integer)` or `NumSequence{T,q}(seq)`
2. **Random Construction**: `NumSequence{T,q}(L::Integer)` or `NumSequence(L::Integer, q::Integer; T=IntType)`.
    Construct a `NumSequence` of length `L` with random integers of type `T` in the range `[1, q]`.


# Examples

```julia-repl
julia> seq = [1, 2, 3, 4]
julia> num_seq = NumSequence(seq, 4)

julia> random_seq = NumSequence(10, 5; T=Int8)


julia> max_value = num_seq.q  # Returns 4


julia> copied_seq = copy(num_seq)
```
"""
@kwdef mutable struct NumSequence{T<:Integer,q} <: AbstractSequence
    seq::Vector{T}
    function NumSequence{T,q}(seq::AbstractVector) where {T<:Integer,q}
        @argcheck q isa Integer "Expect `Integer` for maximum value `q`. Instead $q"
        @argcheck all(x -> 0 < x <= q, seq) "Expect `0 < x < q=$q` for all elements."
        q_convert = convert(T, q) # can potentially fail if say T==Int8 and q very large.
        return new{T,q_convert}(seq)
    end
end

## Constructors
NumSequence(seq::AbstractVector{T}, q::Integer) where {T} = NumSequence{T,q}(seq)
function NumSequence(seq::AbstractVector)
    err = ArgumentError("Provide a maximum value `q`.")
    throw(err)
    return nothing
end

function NumSequence{T,q}(L::Integer) where {T,q}
    seq = rand(T(1):T(q), L)
    return NumSequence{T,q}(seq)
end
NumSequence(L::Integer, q::Integer; T=IntType) = NumSequence{T,q}(L)

Base.copy(x::NumSequence{T,q}) where {T,q} = NumSequence(copy(x.seq), q)
function Base.copy!(dest::NumSequence{T,q}, source::NumSequence{T,q}) where {T,q}
    @argcheck length(source) == length(dest)
    for (i,a) in enumerate(sequence(source))
        dest.seq[i] = a
    end
    return dest
end

function Base.getproperty(x::NumSequence{T,q}, sym::Symbol) where {T,q}
    if sym == :q
        return q
    elseif hasproperty(x, sym)
        return getfield(x, sym)
    end
    throw(ErrorException("type NumSequence has no field $sym"))
end
#===========================================================================#
########################## Converting to Alignment ##########################
#===========================================================================#

"""
    Alignment(sequences; alphabet, names, as_codons=true)

Construct a `BioSequenceMappings.Alignment` from a set of sequences.
If the sequences are `AASequence`, `alphabet` defaults to `aa_alphabet`.
If they are `CodonSequence`, `as_codons` can be used to decide whether the
alignment should store codons or amino acids. `alphabet` can be determined automatically
from this.
"""
function Alignment(
    S::AbstractVector{T};
    names=nothing,
    as_codons=true,
    alphabet=_sequence_alphabet(T; as_codons),
) where {T<:AbstractSequence}
    if !allequal(Iterators.map(length, S))
        error("Sequences do not have the same length")
    end
    if !isnothing(names) && length(names) != length(S)
        error("Got $(length(names)) names but $(length(S)) sequences.")
    end

    data = hcat([sequence(s; as_codons) for s in S]...)

    return Alignment(data, alphabet; names)
end

function genetic_code(A::Alignment, alphabet=aa_alphabet)
    if A.alphabet != codon_alphabet
        error("""
            Function should be called on alignment of codon sequences.
            Instead `A.alphabet`: $(A.alphabet)
        """)
    end
    S = map(x -> genetic_code.(x), A) # this translates to default aa_alphabet
    if alphabet != aa_alphabet
        # translate to requested alphabet if needed
        S = map(x -> BioSequenceMappings.translate(x, aa_alphabet, alphabet))
    end
    return Alignment(S, alphabet; A.names, A.weights)
end

#==================#
####### Misc #######
#==================#

"""
    hamming(x::AbstractSequence, y::AbstractSequence)
    hamming(x::CodonSequence, y::CodonSequence; source=:codon, kwargs...)
"""
function BioSequenceMappings.hamming(
    x::AbstractSequence, y::AbstractSequence; source=nothing, kwargs...
)
    # source kwarg to allow blind use of hamming
    return hamming(x.seq, y.seq; kwargs...)
end
function BioSequenceMappings.hamming(
    x::CodonSequence, y::CodonSequence; source=:codon, kwargs...
)
    return if source == :codon
        hamming(x.seq, y.seq; kwargs...)
    elseif source == :aa
        hamming(x.aaseq, y.aaseq; kwargs...)
    else
        error("Valid `source` values: `:codon` or `:aa`. Instead $source")
    end
end

function intvec_to_sequence(s::AbstractVector{<:Integer}; v=true)
    q = maximum(s)
    return if q < 21 || q > 65
        v && @info "Assume sequence $s is a `NumSequence`"
        NumSequence(s)
    elseif q == 21
        v && @info "Assume sequence $s is an `AASequence`"
        AASequence(s)
    else
        v && @info "Assume sequence $s is a `CodonSequence`"
        CodonSequence(s, q; source=:codon)
    end
end

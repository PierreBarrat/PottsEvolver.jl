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

#=
Methods that a subtype should implement
- sequence: access integer vector (necessary if field is not called `seq`)
- copy !NECESSARY!
- equality and hash
- indexing
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
        q = Q_AA
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
function AASequence(rng::AbstractRNG, L::Integer; T=IntType)
    return AASequence(rand(rng, T(1):T(Q_AA), L))
end
"""
    AASequence(s::AbstractString)

Build an `AASequence` from a string of amino acid symbols, *e.g.* `"AC-DE"`.
Symbols outside of `PottsEvolver.aa_alphabet` raise an error.
"""
function AASequence(s::AbstractString; T=IntType)
    seq = map(collect(s)) do c
        @argcheck in(c, symbols(aa_alphabet)) """
            Symbol '$c' is not an amino acid. Expected one of \
            "$(prod(symbols(aa_alphabet)))".
            """
        return aa_alphabet(c)
    end
    return AASequence(convert(Vector{T}, seq))
end
AASequence(L::Integer; T=IntType) = AASequence(Random.default_rng(), L; T)
AASequence{T}(rng::AbstractRNG, L::Integer) where {T<:Integer} = AASequence(rng, L; T)
AASequence{T}(L::Integer) where {T<:Integer} = AASequence(L; T)

#=============================================#
################ CodonSequence ################
#=============================================#

mutable struct CodonSequence{T<:Integer} <: AbstractSequence
    seq::Vector{T} # the codons
    aaseq::Vector{T} # the translation
    function CodonSequence(seq::Vector{T}, aaseq::Vector{T}) where {T}
        qc = Q_CODON
        qaa = Q_AA
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
function CodonSequence(
    seq::AbstractVector{T}; source=:aa, rng=Random.default_rng()
) where {T<:Integer}
    return if source == :aa
        codons = map(x -> reverse_code_rand(x; rng), seq)
        CodonSequence(convert(Vector{T}, codons), convert(Vector{T}, seq))
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
function CodonSequence(rng::AbstractRNG, L::Int; source=:aa, T=IntType)
    # Base function
    return if source == :aa
        CodonSequence(rand(rng, T(1):T(Q_AA), L); source, rng)
    elseif source == :codon
        codons = T.(rand(rng, coding_codons, L))
        CodonSequence(codons; source=:codon, rng) # rng useless in this case
    else
        error("Unknown `source` $(source). Use `:aa` or `:codon`.")
    end
end
function CodonSequence(L::Int; source=:aa, T=IntType)
    return CodonSequence(Random.default_rng(), L; source, T)
end
function CodonSequence{T}(rng::AbstractRNG, L::Int; kwargs...) where {T<:Integer}
    return CodonSequence(rng, L; T, kwargs...)
end
function CodonSequence{T}(L::Int; kwargs...) where {T<:Integer}
    return CodonSequence{T}(Random.default_rng(), L, kwargs...)
end
"""
    CodonSequence(s::AbstractString)

Build a `CodonSequence` from a string of nucleotides, *e.g.* `"ATGAAA"`.
The length of `s` must be a multiple of three. Gap codons are written `"---"`;
mixing gaps and nucleotides inside a codon is a frameshift and is rejected.
"""
function CodonSequence(s::AbstractString; T=IntType)
    @argcheck length(s) % 3 == 0 """
        Length of a nucleotide string must be a multiple of 3. Instead $(length(s)).
        """
    codons = map(Iterators.partition(s, 3)) do chunk
        codon = Codon(join(chunk))
        @argcheck isvalid(codon) """
            Invalid codon "$(join(chunk))": expected three nucleotides or three gaps.
            """
        return codon_alphabet(codon)
    end
    return CodonSequence(convert(Vector{T}, codons); source=:codon)
end

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

function NumSequence{T,q}(rng::AbstractRNG, L::Integer) where {T,q}
    seq = rand(rng, T(1):T(q), L)
    return NumSequence{T,q}(seq)
end
NumSequence{T,q}(L::Integer) where {T,q} = NumSequence{T,q}(Random.default_rng(), L)
NumSequence(rng::AbstractRNG, L::Integer, q::Integer; T=IntType) = NumSequence{T,q}(rng, L)
NumSequence(L::Integer, q::Integer; T=IntType) = NumSequence(Random.default_rng(), L, q; T)

Base.copy(x::NumSequence{T,q}) where {T,q} = NumSequence(copy(x.seq), q)
function Base.copy!(dest::NumSequence{T,q}, source::NumSequence{T,q}) where {T,q}
    @argcheck length(source) == length(dest)
    for (i, a) in enumerate(sequence(source))
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
#=========================================================================#
########################## Converting to String ##########################
#=========================================================================#

"""
    string(s::AASequence)
    string(s::CodonSequence; as_aa=false)

Return the symbolic representation of `s`.
For a `CodonSequence`, the default is the nucleotide string (three characters per position);
use `as_aa=true` to get the translated amino acid string instead.

Inverse of the `AASequence(::AbstractString)` / `CodonSequence(::AbstractString)`
constructors, which is how sequences read from a fasta file (using *e.g.* `FASTX`) enter
`PottsEvolver`.
"""
Base.string(s::AASequence) = String(map(aa_alphabet, s.seq))

function Base.string(s::CodonSequence; as_aa=false)
    return if as_aa
        String(map(aa_alphabet, s.aaseq))
    else
        join(join(bases(codon_alphabet(c))) for c in s.seq)
    end
end

function Base.string(::NumSequence)
    return throw(
        ArgumentError("""
        A `NumSequence` has no symbolic representation and cannot be converted to a `String`.
        Use `AASequence` or `CodonSequence`, or access the integers with `sequence(s)`.
        """),
    )
end

#==================#
####### Misc #######
#==================#

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

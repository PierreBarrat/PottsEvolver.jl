"""
    SequenceSample{S<:AbstractSequence,T<:Integer}

A set of sampled sequences stored as a matrix of integers, with one label per sequence.

```
data::Matrix{T}         # `L x M`: sequences are stored in *columns*
labels::Vector{String}  # one per sequence, in the same order as the columns
```

The type parameter `S` records what the integers mean — `AASequence`, `CodonSequence` or
`NumSequence` — which is what lets [`write_fasta`](@ref) and `PottsEvolver.translate`
know how to interpret them. A codon sample that has been translated is tagged `AASequence`,
since its integers are then amino acids.

## Methods

- `size(X)` returns `(L, M)`; `length(X)` returns the number of sequences `M`.
- `X[i]` returns a view of the `i`-th sequence; `X["label"]` looks a sequence up by label.
- `for s in X` iterates over sequences.
"""
struct SequenceSample{S<:AbstractSequence,T<:Integer}
    data::Matrix{T}
    labels::Vector{String}
    function SequenceSample{S,T}(data, labels) where {S<:AbstractSequence,T<:Integer}
        @argcheck size(data, 2) == length(labels) """
            Got $(length(labels)) labels for $(size(data, 2)) sequences.
            """
        return new{S,T}(Matrix{T}(data), collect(String, labels))
    end
end

#==================================#
########### Constructors ###########
#==================================#

# For a `CodonSequence`, `as_codons=false` stores the translation: the sample then holds
# amino acids and must be tagged as such.
_sample_type(::Type{S}, as_codons) where {S<:AbstractSequence} = S
function _sample_type(::Type{CodonSequence{T}}, as_codons) where {T}
    return as_codons ? CodonSequence{T} : AASequence{T}
end

"""
    SequenceSample(sequences; labels=nothing, as_codons=true)

Build a `SequenceSample` from a vector of sequences, which must all have the same length.
`labels` defaults to the index of each sequence.
For `CodonSequence`, `as_codons=false` stores the translated amino acids instead of codons.
"""
function SequenceSample(
    sequences::AbstractVector{S}; labels=nothing, as_codons=true
) where {S<:AbstractSequence}
    @argcheck !isempty(sequences) "Cannot build a `SequenceSample` from no sequences"
    @argcheck allequal(Iterators.map(length, sequences)) """
        Sequences do not have the same length
        """
    labels = isnothing(labels) ? (1:length(sequences)) : labels
    @argcheck length(labels) == length(sequences) """
        Got $(length(labels)) labels but $(length(sequences)) sequences.
        """

    T = eltype(first(sequences))
    L = length(first(sequences))
    data = Matrix{T}(undef, L, length(sequences))
    for (m, s) in enumerate(sequences)
        data[:, m] .= sequence(s; as_codons)
    end

    return SequenceSample{_sample_type(S, as_codons),T}(data, string.(labels))
end

#=========================================#
########### Iterating / Indexing ##########
#=========================================#

"""
    size(X::SequenceSample)

Return a tuple with (in order) the length of the sequences and their number.
"""
Base.size(X::SequenceSample) = size(X.data)
Base.size(X::SequenceSample, dim) = size(X.data, dim)

"""
    length(X::SequenceSample)

Return the number of sequences in `X`.
"""
Base.length(X::SequenceSample) = size(X.data, 2)

Base.getindex(X::SequenceSample, i::Integer) = view(X.data, :, i) # returns a view!
function Base.getindex(X::SequenceSample, label::AbstractString)
    i = findfirst(==(label), X.labels)
    isnothing(i) && throw(KeyError(label))
    return view(X.data, :, i)
end
Base.firstindex(X::SequenceSample) = 1
Base.lastindex(X::SequenceSample) = length(X)
Base.keys(X::SequenceSample) = LinearIndices(1:length(X))

Base.iterate(X::SequenceSample) = iterate(eachcol(X.data))
Base.iterate(X::SequenceSample, state) = iterate(eachcol(X.data), state)
Base.eltype(::SequenceSample{S,T}) where {S,T} = AbstractVector{T}

#==================#
####### Misc #######
#==================#

function Base.:(==)(X::SequenceSample, Y::SequenceSample)
    return typeof(X) == typeof(Y) && X.data == Y.data && X.labels == Y.labels
end

Base.copy(X::SequenceSample{S,T}) where {S,T} = SequenceSample{S,T}(copy(X.data), copy(X.labels))

"""
    sequence_type(X::SequenceSample)

Return the type of sequence that `X` stores, *e.g.* `CodonSequence{Int64}`.
"""
sequence_type(::SequenceSample{S,T}) where {S,T} = S

function Base.show(io::IO, X::SequenceSample{S,T}) where {S,T}
    L, M = size(X)
    return print(io, "SequenceSample{$S}: M=$M sequences of length L=$L")
end
function Base.show(io::IO, x::MIME"text/plain", X::SequenceSample{S,T}) where {S,T}
    L, M = size(X)
    println(io, "SequenceSample{$S}: M=$M sequences of length L=$L - shown as `MxL` matrix")
    return show(io, x, X.data')
end

# Rebuild a sequence object from one column, used when writing to fasta.
_sequence_from_column(::SequenceSample{S}, col) where {S<:AASequence} = AASequence(collect(col))
function _sequence_from_column(::SequenceSample{S}, col) where {S<:CodonSequence}
    return CodonSequence(collect(col); source=:codon)
end
function _sequence_from_column(::SequenceSample{S}, col) where {T,q,S<:NumSequence{T,q}}
    return NumSequence(collect(col), q)
end

"""
    translate(X::SequenceSample{<:CodonSequence})

Translate a sample of codons into the corresponding sample of amino acids.
"""
function translate(X::SequenceSample{S,T}) where {S<:CodonSequence,T}
    data = map(x -> convert(T, genetic_code(x)), X.data)
    return SequenceSample{AASequence{T},T}(data, copy(X.labels))
end

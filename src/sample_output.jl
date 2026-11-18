"""
    SequenceSample{S<:AbstractSequence,T<:Integer}

A set of sampled sequences stored as a matrix of integers, with one label per sequence.

```
data::Matrix{T}         # `L x M`: sequences are stored in *columns*
labels::Vector{String}  # one per sequence, in the same order as the columns
```

The type parameter `S` records what the integers mean — `AASequence`, `RNASequence`,
`CodonSequence` or `NumSequence` — which is what lets [`write_fasta`](@ref) and `PottsEvolver.translate`
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

function Base.copy(X::SequenceSample{S,T}) where {S,T}
    return SequenceSample{S,T}(copy(X.data), copy(X.labels))
end

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

"""
    sequences(X::SequenceSample)

Rebuild the sequences of `X` as a `Vector` of `AbstractSequence`, of the type that `X` is
tagged with: the inverse of `SequenceSample(sequences)`, dropping the labels.

Needed by the functions that take a plain vector of sequences, *e.g.*
`average_transition_rate`.
"""
function sequences(X::SequenceSample)
    return map(m -> _sequence_from_column(X, X[m]), 1:length(X))
end

# Rebuild a sequence object from one column, used when writing to fasta.
function _sequence_from_column(::SequenceSample{S}, col) where {S<:AASequence}
    return AASequence(collect(col))
end
function _sequence_from_column(::SequenceSample{S}, col) where {S<:RNASequence}
    return RNASequence(collect(col))
end
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

"""
    regroup_samples(data::AbstractVector{<:SequenceSample})

Transpose a vector of samples against their labels.

Take `M` samples of size `L x N`, all covering the same `N` labels, and return `N` samples of
size `L x M`: one per label, in the order of `first(data).labels`. In other words,
`data[m]["label"]` becomes `out[n][m]`, where `n` is the position of `"label"` in that order.
The labels of each output sample are the indices `1:M` into `data`.

The `N` labels must be the same *set* in every element of `data`, but not necessarily in the
same order: each sample is looked up by label, not by position.

Used to regroup the output of repeated runs. For a vector of chains, labels are sampling times,
so the output is one sample per time, with chains as columns; see
`parallel_chains`. For a vector of tree runs, labels are node names, and
`pernode_alignment` wraps this to key the result by node label instead.
"""
function regroup_samples(data::AbstractVector{S}) where {S<:SequenceSample}
    @argcheck !isempty(data) "Cannot regroup an empty vector of samples"
    @argcheck allequal(X -> sort(X.labels), data) """
        All samples must cover the same set of labels
        """
    @argcheck allequal(X -> size(X, 1), data) """
        All samples must hold sequences of the same length
        """

    labels = first(data).labels
    L = size(first(data), 1)
    M = length(data)
    T = eltype(first(data).data)

    # one label -> column lookup per sample, instead of a linear scan per (label, sample)
    columns = map(X -> Dict(l => i for (i, l) in enumerate(X.labels)), data)

    return map(labels) do label
        D = Matrix{T}(undef, L, M)
        for m in 1:M
            D[:, m] .= view(data[m].data, :, columns[m][label])
        end
        return S(D, string.(1:M))
    end
end

"""
    hamming(x, y; normalize=true, positions=nothing, exclude_state=nothing)

Hamming distance between vectors of integers `x` and `y`.
Only sites in `positions` are considered; if `exclude_state` is given, sites where either
sequence is in that state are skipped entirely (they count towards neither the distance nor
the normalization).
"""
function hamming(
    X::AbstractVector{<:Integer},
    Y::AbstractVector{<:Integer};
    normalize=true,
    positions=nothing,
    exclude_state=nothing,
)
    @argcheck length(X) == length(Y) """
        Expect vectors of same length. Instead $(length(X)) != $(length(Y))
        """

    positions = isnothing(positions) ? eachindex(X) : positions
    H = 0
    Z = 0
    if isnothing(exclude_state)
        for i in positions
            Z += 1
            X[i] != Y[i] && (H += 1)
        end
    else
        for i in positions
            if X[i] != exclude_state && Y[i] != exclude_state
                Z += 1
                X[i] != Y[i] && (H += 1)
            end
        end
    end
    return normalize ? H / Z : H
end

"""
    hamming(x::AbstractSequence, y::AbstractSequence)
    hamming(x::CodonSequence, y::CodonSequence; source=:codon, kwargs...)

For `CodonSequence`, `source` picks whether the distance is computed on codons (`:codon`)
or on the translated amino acids (`:aa`).
"""
function hamming(x::AbstractSequence, y::AbstractSequence; source=nothing, kwargs...)
    # source kwarg to allow blind use of hamming
    return hamming(x.seq, y.seq; kwargs...)
end
function hamming(x::CodonSequence, y::CodonSequence; source=:codon, kwargs...)
    return if source == :codon
        hamming(x.seq, y.seq; kwargs...)
    elseif source == :aa
        hamming(x.aaseq, y.aaseq; kwargs...)
    else
        error("Valid `source` values: `:codon` or `:aa`. Instead $source")
    end
end

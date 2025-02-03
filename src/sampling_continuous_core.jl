#============================================================#
##################### Energy differences #####################
#============================================================#


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

In the first form: the new reference `b` sequence is equal to `refseq` at positions different from `i`, and to `x` at position `i`.
For all mutations `(j,z)`, giving rise to sequence `c`, compute `E(c)-E(b)` and store the result in output matrix `ΔE` of dimensions `qxL`.
Keyword argument `ΔE` can be provided to avoid allocation.


In the second form: compute the energy difference between any sequence `b` at distance one from `a` and `a` itself.
The output is a matrix of dimensions `qxL`, always allocated in the function.

*Note*: For `CodonSequence`, the output is of size `65xL` where `L` is the length.
Stop lead to infinite energy. Providing a sequence that contains a stop codon is undefined.
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

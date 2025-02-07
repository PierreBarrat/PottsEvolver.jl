"""
    PottsGraph{T}

- Array `J` of dimensions `q x q x L x L` and eltype `T`
- Array `h` of dimensions `q x L` and eltype `T`
- Inverse temperature `β`
- `alphabet`
"""
@kwdef mutable struct PottsGraph{T<:AbstractFloat}
    J::Array{T,4}
    h::Array{T,2}
    β::T = 1.0
    alphabet::Union{Nothing,Alphabet{Char,<:Integer}} = aa_alphabet
    function PottsGraph(J::Array{T,4}, h::Array{T,2}, β, alphabet) where {T}
        @argcheck size(h, 1) == size(J, 1) == size(J, 2) """
            Inconsistent sizes for `J` and `h`: $(size(J)) - $(size(h))
            """
        @argcheck size(h, 2) == size(J, 3) == size(J, 4) """
            Inconsistent sizes for `J` and `h`: $(size(J)) - $(size(h))
            """
        @argcheck isnothing(alphabet) || size(h, 1) == length(alphabet) """
            Inconsistent alphabet size: $(length(alphabet)) - h: $(size(h))
            """
        return new{T}(J, h, β, alphabet)
    end
end

"""
    PottsGraph(L, q[, T]; init = :null, alphabet)

Return a `PottsGraph{T}` of the required size.
- `init == :null`: parameters are intialized to zero.
- `init == :rand`: parameters are randomly sampled using `Jrand` and `hrand` keywords.

`alphabet` is `aa_alphabet` if `q=21`, `nothing` otherwise.

## Random initialization

- `hrand` should be a function `q -> h`.
- `Jrand` should be a function `q -> J`.

`Jrand` does not have to return a symetric matrix.
The output matrix is made symetric with zeroes on the diagonal blocks.

"""
function PottsGraph(
    L,
    q,
    T=FloatType;
    init=:null,
    Jrand=N -> 1 / L * randn(N, N),
    hrand=N -> 1 / sqrt(L) * randn(N),
    alphabet=(q == 21 ? aa_alphabet : nothing),
)
    return if init == :null
        PottsGraph(; J=zeros(T, q, q, L, L), h=zeros(T, q, L), alphabet)
    elseif init == :rand
        J, h = _random_graph(L, q)
        PottsGraph(; J, h, alphabet)
    end
end

function Base.size(g::PottsGraph)
    return (L=size(g.h, 2), q=size(g.h, 1))
end

function _random_graph(L, q)
    h = reshape(randn(L * q), q, L) / sqrt(L)
    J = reshape(randn(L * L * q * q), q, q, L, L) / L
    for i in 1:L
        J[:, :, i, i] .= 0
        for j in (i + 1):L
            J[:, :, i, j] .= J[:, :, j, i]'
        end
    end
    return J, h
end

"""
    couplings_as_matrix(J::Array{<:T, 4})
    matrix_as_couplings(Jm::Matrix{T}, L, q)

First version: convert `(q, q, L, L)` to `(L*q, L*q)`
Second version: opposite
"""
function couplings_as_matrix(J::Array{<:T,4}) where {T<:AbstractFloat}
    q, L = size(J, 1), size(J, 3)
    return reshape(permutedims(J, (1, 3, 2, 4)), L * q, L * q)
end
"""
    couplings_as_matrix(J::Array{<:T, 4})
    matrix_as_couplings(Jm::Matrix{T}, L, q)

First version: convert `(q, q, L, L)` to `(L*q, L*q)`
Second version: opposite
"""
function matrix_as_couplings(Jm::Matrix{T}, L, q) where {T<:AbstractFloat}
    return permutedims(reshape(Jm, q, L, q, L), (1, 3, 2, 4))
end

#==========================================#
############### Gauge change ###############
#==========================================#

const _gauges = [:zero_sum, :lattice_gas]
"""
    set_gauge!(g::PottsGraph, gauge = :zero_sum)
"""
function set_gauge!(g::PottsGraph, gauge=:zero_sum)
    if !in(gauge, _gauges)
        error("Gauge must be in $_gauges. Instead $gauge")
    end

    if gauge == :zero_sum
        set_gauge_zero_sum!(g)
    elseif gauge == :latticegas
        set_gauge_lattice_gas!(g)
    else
        error("Unknown gauge $gauge")
    end

    return g
end

function set_gauge_zero_sum!(g::PottsGraph{T}) where {T}
    L, q = (size(g).L, size(g).q)
    h0 = zeros(T, q, L)
    J0 = zeros(T, q, q, L, L)
    t = zeros(T, q, q)

    for i in 1:L, j in (i + 1):L
        t .= g.J[:, :, i, j]
        J0[:, :, i, j] .=
            t - repeat(mean(t; dims=1), q, 1) .- repeat(mean(t; dims=2), 1, q) .+ mean(t)
        J0[:, :, j, i] .= J0[:, :, i, j]'
    end

    for i in 1:L
        h0[:, i] .= g.h[:, i] .- mean(g.h[:, i])
        for j in 1:L
            h0[:, i] .+= mean(g.J[:, :, i, j]; dims=2) .- mean(g.J[:, :, i, j])
        end
    end

    g.J .= J0
    g.h .= h0
    return g
end

function set_gauge_lattice_gas!(g)
    error("Not implemented yet")
    return g
end

#==================#
####### Misc #######
#==================#

"""
    energy(s, g::PottsGraph)

Energy of sequence `s` in `g`.
"""
function energy(s::AbstractVector{<:Integer}, g::PottsGraph)
    (; L, q) = size(g)
    E = 0.0
    for i in 1:L
        E += g.h[s[i], i]
        for j in (i + 1):L
            E += g.J[s[i], s[j], i, j]
        end
    end
    return -E
end
energy(s::AbstractSequence, g) = energy(sequence(s), g)
energy(s::CodonSequence, g) = energy(s.aaseq, g)

function Base.show(io::IO, g::PottsGraph{T}) where {T}
    (; L, q) = size(g)
    return print(io, "PottsGraph{$T} (L=$L, q=$q)")
end
function Base.show(io::IO, x::MIME"text/plain", g::PottsGraph{T}) where {T}
    (; L, q) = size(g)
    return print(io, "PottsGraph{$T}: dimensions (L=$L, q=$q) -- β=$(g.β) -- $(g.alphabet)")
end

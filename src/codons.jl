#==================================================================#
####################### Alphabets and Codons #######################
#==================================================================#

const nucleotides = ['A', 'C', 'G', 'T']
const _nucleotides_gap = ['A', 'C', 'G', 'T', '-']

#=
`aa_alphabet` and `codon_alphabet` are *functions*, not objects: they map a symbol to its
index and an index back to its symbol, dispatching on the type of their argument.
Both mappings are fixed, so there is nothing to configure and nothing to carry around.
Use `symbols(alphabet)` for the vector of symbols and `Q_AA`/`Q_CODON` for their number.
=#

const AA_SYMBOLS = collect("-ACDEFGHIKLMNPQRSTVWY")
const AA_INDEX = Dict{Char,IntType}(c => i for (i, c) in enumerate(AA_SYMBOLS))
"""
    Q_AA

Number of amino acid symbols: the 20 amino acids and the gap.
"""
const Q_AA = IntType(length(AA_SYMBOLS))

"""
    aa_alphabet(c::AbstractChar) -> Integer
    aa_alphabet(i::Integer) -> Char

Map an amino acid symbol to its index, or an index back to its symbol.
"""
function aa_alphabet(c::AbstractChar)
    i = get(AA_INDEX, c, nothing)
    isnothing(i) && throw(
        ArgumentError(
            "'$c' is not an amino acid symbol - expected one of \"$(prod(AA_SYMBOLS))\""
        ),
    )
    return i
end
aa_alphabet(i::Integer) = AA_SYMBOLS[i]

@kwdef struct Codon
    b1::Char
    b2::Char
    b3::Char
    function Codon(b1, b2, b3)
        @argcheck b1 in _nucleotides_gap && b2 in _nucleotides_gap && b3 in _nucleotides_gap
        return new(b1, b2, b3)
    end
end
function Codon(s::AbstractString)
    @argcheck length(s) == 3
    return Codon(s[1], s[2], s[3])
end

"""
    bases(codon)

Return iterator on the bases of `codon`.
"""
bases(codon::Codon) = Iterators.map(i -> getfield(codon, i), 1:3)

const CODON_SYMBOLS = let
    nt = nucleotides
    C = vec(map(x -> Codon(x...), Iterators.product(nt, nt, nt))) # AAA CAA GAAA etc... (first changes fastest)
    pushfirst!(C, Codon('-', '-', '-')) # the gap codon therefore has index 1
    C
end
const CODON_INDEX = Dict{Codon,IntType}(c => i for (i, c) in enumerate(CODON_SYMBOLS))
"""
    Q_CODON

Number of codon symbols: the 64 nucleotide triplets and the gap codon.
"""
const Q_CODON = IntType(length(CODON_SYMBOLS))

"""
    codon_alphabet(c::Codon) -> Integer
    codon_alphabet(i::Integer) -> Codon

Map a codon to its index, or an index back to its codon.
"""
function codon_alphabet(c::Codon)
    i = get(CODON_INDEX, c, nothing)
    isnothing(i) && throw(ArgumentError("$c is not in the codon alphabet"))
    return i
end
codon_alphabet(i::Integer) = CODON_SYMBOLS[i]

"""
    symbols(alphabet)

Return the vector of symbols used by `aa_alphabet` or `codon_alphabet`.
"""
symbols(::typeof(aa_alphabet)) = AA_SYMBOLS
symbols(::typeof(codon_alphabet)) = CODON_SYMBOLS

#==========================================#
############### Genetic code ###############
#==========================================#

#! format: off
const aa_order = [
'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
]
#! format: on

# Dictionary from Codon to Char
const _genetic_code_struct = let
    code = Dict{Codon,Char}()
    i = 1
    for a in nucleotides, b in nucleotides, c in nucleotides # last changes fastest, this is ok
        code[Codon(a, b, c)] = aa_order[i]
        i += 1
    end
    code[Codon('-', '-', '-')] = '-'
    code
end
# Dictionary from Codon to Int
const _genetic_code_integers = let
    code = Dict{IntType,Union{Nothing,IntType}}()
    for (codon, aa) in _genetic_code_struct
        code[codon_alphabet(codon)] = aa == '*' ? nothing : aa_alphabet(aa)
    end
    code
end
"""
    genetic_code(x::Integer)

Translate the `i`th codon and return the index of the corresponding amino acid, using
the default `aa_alphabet`
"""
function genetic_code(codon::T) where {T<:Integer}
    aa = _genetic_code_integers[codon]
    return isnothing(aa) ? aa : T(aa)
end
"""
    genetic_code(c::Codon)

Translate `c` and return the amino acid as a `Char`.
"""
genetic_code(codon::Codon) = _genetic_code_struct[codon]

const _reverse_code_integers = let
    rcode = Vector{Vector{IntType}}(undef, Q_AA)
    for aa in 1:Q_AA
        rcode[aa] = IntType[
            i for (i, c) in enumerate(CODON_SYMBOLS) if genetic_code(c) == aa_alphabet(aa)
        ]
    end
    rcode
end
const _reverse_code_struct = let
    rcode = Dict{Char,Vector{Codon}}()
    for aa in symbols(aa_alphabet)
        rcode[aa] = map(codon_alphabet, _reverse_code_integers[aa_alphabet(aa)])
    end
    rcode
end

"""
    reverse_code(aa)

Return the set of codons coding for `aa`.
"""
reverse_code(aa::T) where {T<:Integer} = _reverse_code_integers[aa]
reverse_code(aa::AbstractChar) = _reverse_code_struct[aa]
"""
    reverse_code_rand(aa; rng)

Return a random codon coding for `aa`
"""
function reverse_code_rand(aa::Integer; rng=Random.default_rng())
    codons = reverse_code(aa)
    # @info codons
    if isempty(codons)
        error("No codon corresponds to amino acid $aa ($aa_alphabet(aa))")
    end
    return rand(rng, codons)
end
function reverse_code_rand(aa::AbstractChar; kwargs...)
    return codon_alphabet(reverse_code_rand(aa_alphabet(aa)); kwargs...)
end

#========================================================================#
######################### Codon helper functions #########################
#========================================================================#

Base.show(io::IO, c::Codon) = print(io, "\"$(c.b1)$(c.b2)$(c.b3)\"")

function Base.show(io::IO, x::MIME"text/plain", c::Codon)
    # `Codon` accepts a gap in any position, so it can hold frameshifts that are absent
    # from the genetic code: check before looking anything up.
    isvalid(c) || return println(io, "Codon \"$(c.b1)$(c.b2)$(c.b3)\": invalid (frameshift)")

    aa, iaa = if isstop(c)
        '*', "STOP"
    else
        a = genetic_code(c)
        a, aa_alphabet(a)
    end
    return println(
        io, "Codon $(codon_alphabet(c)): \"$(c.b1)$(c.b2)$(c.b3)\" --> $aa($iaa)"
    )
end

# Only all gaps or all nt codons are valid
# Anything else means a frameshift and I do not deal with that here
isgap(c::Codon) = all(==('-'), bases(c))
function Base.isvalid(c::Codon)
    return if isgap(c)
        true
    elseif all(in(nucleotides), bases(c))
        true
    else
        false
    end
end
isstop(c::Codon) = genetic_code(c) == '*'
iscoding(c::Codon) = !isgap(c) && !isstop(c) && isvalid(c)

# pre-computed for faster calculations on integers

const stop_codon_indices = IntType.(findall(isstop, symbols(codon_alphabet)))
isstop(i::Integer) = in(i, stop_codon_indices)

# isgap should also work on amino acids : pass alphabet
const gap_codon_index = IntType(findfirst(isgap, symbols(codon_alphabet)))
isgap_codon(i::Integer) = (i == gap_codon_index)
isgap(c::AbstractChar) = (c == '-')

iscoding(c::Integer) = !isgap_codon(c) && !isstop(c)
const coding_codons = IntType.(findall(iscoding, 1:Q_CODON))

# Number of non-stop / non-gaps codons
const n_aa_codons = count(c -> !isgap(c) && !isstop(c), symbols(codon_alphabet))

#===========================================================================#
########################## Codon accessibility map ##########################
#===========================================================================#

#=
For a given codon `c` and a given position `i` inside the codon,
what other codons are accessible by one mutation?

`_codon_access_map`: let `c::Int` be a codon and `i ∈ [1,3]` a position. `codon_access_map[c,i]` returns a tuple with:
    - the list of codons accessible by mutating `c` at position `i`, as integers. Stop codons or invalid codons are filtered out.
    - the list of corresponding amino acids, again as integers
!!! Only nucleotide mutations are considered. For this reason the gap codon never appears in this dictionary.
!!! I consider here that a codon is always accessible from itself! *i.e.* `c` will appear in the list
=#

function _build_codon_access_map()
    M = Dict{Tuple{IntType,IntType},Tuple{ReadOnlyVector{IntType},ReadOnlyVector{IntType}}}()
    for c in 1:Q_CODON, i in 1:3
        codon = codon_alphabet(c)
        if !isgap(codon) && isvalid(codon)
            accessible_codons = map(nucleotides) do a
                nts = [codon.b1, codon.b2, codon.b3]
                nts[i] = a
                return codon_alphabet(Codon(nts...))
            end
            filter!(!isstop, accessible_codons)
            M[c, i] = (
                ReadOnlyArray(accessible_codons),
                ReadOnlyArray(map(genetic_code, accessible_codons)),
            )
        end
    end
    return M
end
const _codon_access_map = _build_codon_access_map()

#=
Similar to the above, with the following differences.
- Stores all codons accessible by any nucleotide or gap mutation. Consequently, keys are `IntType` (just the codon)
- Gap mutations are counted: the gap codon is accessible from all codons, and all codons are accessible from the gap.
- A codon is not accessible from itself (this is used for continuous time sampling).
=#
function _build_codon_access_map_2()
    M = Dict{IntType,Tuple{ReadOnlyVector{IntType},ReadOnlyVector{IntType}}}()

    for c in 1:Q_CODON
        codon = codon_alphabet(c)
        isstop(codon) && continue

        # Gap case first
        if isgap(codon)
            accessible_codons = collect(1:Q_CODON)
            filter!(c -> iscoding(codon_alphabet(c)), accessible_codons) # remove all non-coding (i.e. gap and stop)
            M[c] = (
                ReadOnlyArray(accessible_codons),
                ReadOnlyArray(map(genetic_code, accessible_codons)),
            )
            continue
        end

        # general case
        accessible_codons = IntType[]
        # check all mutations, filter for coding codons
        for b in 1:3, nt in nucleotides
            nts = collect(bases(codon))
            nts[b] = nt
            new_codon = Codon(nts...)
            if iscoding(new_codon) && new_codon != codon
                push!(accessible_codons, codon_alphabet(new_codon))
            end
        end
        # add gap codon
        push!(accessible_codons, codon_alphabet(Codon("---")))
        # store
        M[c] = (
            ReadOnlyArray(accessible_codons),
            ReadOnlyArray(map(genetic_code, accessible_codons)),
        )
    end

    return M
end

const _codon_access_map_2 = _build_codon_access_map_2()

"""
    accessible_codons(codon, b::Integer)

Return all codons/amino-acids accessible by mutating `codon` at base `b`.
Value returned is a `Tuple` whose first/second elements represent codons/amino-acids.
`codon` itself (and the corresponding amino-acid) is included in the result.

# Examples
```jldoctest
julia> seq = CodonSequence([1,2,3]; source=:codon); # three first codons: gap, AAA, CAA

julia> PottsEvolver.accessible_codons(seq[1], 1) # mutating the gap codon is undefined
(nothing, nothing)

julia> PottsEvolver.accessible_codons(seq[2], 1) # mutating codon 2 at base 1 gives access to 2 others (AAA, CAA, GAA); TAA is stop.
([2, 3, 4], [10, 15, 5])
```
"""
function accessible_codons(codon::Integer, b::Integer)
    return get(_codon_access_map, (codon, b), (nothing, nothing))
end
function accessible_codons(codon::Codon, b::Integer)
    Cs = get(_codon_access_map, (codon_alphabet(codon), b), nothing)
    return if isnothing(Cs)
        Codon[], Char[]
    else
        map(codon_alphabet, Cs[1]), map(aa_alphabet, Cs[2])
    end
end

"""
    accessible_codons(codon)

Return all codons/amino-acids accessible by mutating `codon` at any base.
Value returned is a `Tuple` whose first/second elements represent codons/amino-acids.
`codon` itself (and the corresponding amino-acid) is **not** included in the result.
This differs from the two argument version.

This function differs from the two argument version:
- the gap codon is accessible from all other codons;
- a codon is not accessible from itself;
- any codon is accessible from the gap codon.

All non-gap codons are considered accessible from the gap codon.
"""
function accessible_codons(codon::Integer)
    return get(_codon_access_map_2, codon, (nothing, nothing))
end
function accessible_codons(codon::Codon)
    Cs = get(_codon_access_map_2, codon_alphabet(codon), nothing)
    return if isnothing(Cs)
        Codon[], Char[]
    else
        map(codon_alphabet, Cs[1]), map(aa_alphabet, Cs[2])
    end
end

#===========================================================================#
########################## Amino acid degeneracies ##########################
#===========================================================================#

const _aa_degeneracy = Dict{IntType,FloatType}(
    a => log(length(reverse_code(a))) for a in 1:Q_AA
)
aa_degeneracy(a::Integer) = get(_aa_degeneracy, a, -Inf)

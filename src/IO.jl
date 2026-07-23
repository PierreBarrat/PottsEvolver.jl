#==========================================================#
##################### Writing sequences #####################
#==========================================================#

"""
    write_fasta(file, sequences; labels=nothing, kwargs...)

Write `sequences` to `file` in fasta format.
`labels` defaults to the index of each sequence; extra keyword arguments are forwarded to
`string`, so `as_aa=true` writes a `CodonSequence` vector as amino acids.

`PottsEvolver` does not read fasta files: use a dedicated package and the string
constructors to build sequences, *e.g.*
```julia
using FASTX
seqs = [AASequence(sequence(rec)) for rec in FASTAReader(open(file))]
```
"""
function write_fasta(
    file::AbstractString,
    sequences::AbstractVector{<:AbstractSequence};
    labels=nothing,
    kwargs...,
)
    labels = isnothing(labels) ? (1:length(sequences)) : labels
    @argcheck length(labels) == length(sequences) """
        Got $(length(labels)) labels for $(length(sequences)) sequences.
        """
    open(file, "w") do io
        for (label, seq) in zip(labels, sequences)
            write(io, ">", string(label), "\n", string(seq; kwargs...), "\n")
        end
    end
    return nothing
end

"""
    write_fasta(file, sample::SequenceSample; labels=sample.labels, kwargs...)

Write a sample obtained from `mcmc_sample` to `file` in fasta format.
Labels default to those carried by the sample: sampling times for a chain, node labels for
a tree.
"""
function write_fasta(
    file::AbstractString, sample::SequenceSample; labels=sample.labels, kwargs...
)
    sequences = map(m -> _sequence_from_column(sample, sample[m]), 1:length(sample))
    return write_fasta(file, sequences; labels, kwargs...)
end

#============================================================#
##################### Reading PottsGraph #####################
#============================================================#

"""
    read_graph
    read_potts_graph

Read `PottsGraph` object from a file.
"""
function read_graph(file::AbstractString, T=FloatType)
    first_line = readline(file)
    format = infer_format_from_line(first_line)
    if format == :numerical
        @debug "Numerical format detected"
        return read_graph_numerical(file, T)
    elseif format == :symbolic
        @debug "Symbolic format detected"
        return read_graph_symbolic(file, T)
    else
        throw(ArgumentError("Unknown format for first line of $file\n $(first_line)"))
    end
end
"""
    read_potts_graph

Alias for `read_graph`.
"""
read_potts_graph = read_graph

function read_graph_symbolic(file, T=FloatType)
    ## Go through file twice: first to get L and the alphabet, the second to store parameters
    # The alphabet (and hence q) is inferred from the state letters: the file lists every
    # symbol at each position, so the `h` lines cover exactly the alphabet.
    L = 0
    min_idx = Inf
    chars = Set{Char}()
    for (n, line) in enumerate(eachline(file))
        if !is_valid_line(line, :symbolic)
            throw(ArgumentError("""
                Format problem with line $n in $file
                Expected format `J i j a b` or `h i a`, with symbolic states.\
                Instead $line"""))
        end
        if !isempty(line) && line[1] == 'h'
            s = split(line, " ")
            i = parse(Int, s[2])
            L = max(L, i)
            min_idx = min(min_idx, i)
            push!(chars, s[3][1])
        end
    end
    if min_idx != 1 && min_idx != 0
        throw(
            ArgumentError("Issue with indexing in $file: smallest index found is $min_idx")
        )
    end
    index_style = (min_idx == 0 ? 0 : 1)
    if index_style == 0
        L += 1
    end
    alphabet = _symbolic_alphabet_from_chars(chars)
    q = length(symbols(alphabet))
    @debug "Index style: $index_style"
    @debug L, q, alphabet

    g = PottsGraph(L, q, T)
    for line in eachline(file)
        if line[1] == 'J'
            i, j, a, b, val = parse_coupling_line_symbolic(line, T, alphabet)
            index_style == 0 && (i += 1; j += 1)
            g.J[a, b, i, j] = val
            g.J[b, a, j, i] = val
        elseif line[1] == 'h'
            i, a, val = parse_field_line_symbolic(line, T, alphabet)
            index_style == 0 && (i += 1)
            g.h[a, i] = val
        end
    end

    return g
end

function read_graph_numerical(file, T=FloatType)
    ## Go through file twice: first to get L and q, the second to store parameters
    q = 0
    L = 0
    min_idx = Inf
    index_style = 1
    for (n, line) in enumerate(eachline(file))
        if !is_valid_line(line, :numerical)
            throw(ArgumentError("""
                Format problem with line $n in $file
                Expected format `J i j a b` or `h i a` (with numerical symbols).
                Instead $line"""))
        end
        if !isempty(line) && line[1] == 'h'
            i, a, val = parse_field_line_numerical(line, T)
            if i > L
                L = i
            end
            if a > q
                q = a
            end
            if i < min_idx || a < min_idx
                min_idx = min(i, a)
            end
        end
    end
    if min_idx != 1 && min_idx != 0
        throw(
            ArgumentError("Issue with indexing in $file: smallest index found is $min_idx")
        )
    end
    index_style = (min_idx == 0 ? 0 : 1)
    if index_style == 0
        L += 1
        q += 1
    end

    g = PottsGraph(L, q, T)
    for line in eachline(file)
        if line[1] == 'J'
            i, j, a, b, val = parse_coupling_line_numerical(line, T)
            index_style == 0 && (i += 1; j += 1; a += 1; b += 1)
            g.J[a, b, i, j] = val
            g.J[b, a, j, i] = val
        elseif line[1] == 'h'
            i, a, val = parse_field_line_numerical(line, T)
            index_style == 0 && (i += 1; a += 1)
            g.h[a, i] = val
        end
    end

    return g
end

# Alphabets whose states can appear as single letters in a *symbolic* PottsGraph file.
# Codons are excluded: they are three characters, not one.
const _SYMBOLIC_ALPHABETS = (aa_alphabet, rna_alphabet)

# A graph file lists every symbol at each position, so the set of state letters seen in the
# `h` lines is exactly the alphabet. Match it to decide which one it is.
function _symbolic_alphabet_from_chars(chars)
    cset = Set(chars)
    for alphabet in _SYMBOLIC_ALPHABETS
        Set(symbols(alphabet)) == cset && return alphabet
    end
    throw(
        ArgumentError("""
        State symbols $(sort(collect(cset))) match no known symbolic alphabet.
        Known: $(join(map(a -> "\"$(prod(symbols(a)))\"", _SYMBOLIC_ALPHABETS), ", ")).
        """),
    )
end

# The reverse, for writing: the alphabet is fixed by the number of states.
function _symbolic_alphabet_from_q(q)
    for alphabet in _SYMBOLIC_ALPHABETS
        length(symbols(alphabet)) == q && return alphabet
    end
    throw(
        ArgumentError("""
        No symbolic alphabet has q=$q states, cannot write this graph in symbolic format.
        Known: aa (21), rna (5). Use `format=:numerical` instead.
        """),
    )
end

let
    # union of the letters of all symbolic alphabets; `-` stays leading so it is literal
    letters = prod(unique(reduce(vcat, map(symbols, _SYMBOLIC_ALPHABETS))))
    patterns = Dict(
        :numerical => Regex.(["J [0-9]+ [0-9]+ [0-9]+ [0-9]+", "h [0-9]+ [0-9]+"]),
        :symbolic =>
            Regex.(["J [0-9]+ [0-9]+ [$(letters)]+ [$(letters)]+", "h [0-9]+ [$(letters)]"]),
    )
    global get_line_patterns() = patterns
end

function infer_format_from_line(line)
    pattern_dict = get_line_patterns()
    for (name, patterns) in pattern_dict
        if !isnothing(match(patterns[1], line)) || !isnothing(match(patterns[2], line))
            return name
        end
    end
    return nothing
end
function is_valid_line(line, format::Symbol)
    patterns = get_line_patterns()[format]
    return !isnothing(match(patterns[1], line)) || !isnothing(match(patterns[2], line))
end

# function is_valid_line(line)
#     if isnothing(match(r"J [0-9]+ [0-9]+ [0-9]+ [0-9]+", line)) &&
#         isnothing(match(r"h [0-9]+ [0-9]+", line))
#         return false
#     else
#         return true
#     end
# end
function parse_field_line_symbolic(line, T, alphabet)
    s = split(line, " ")
    i = parse(Int, s[2])
    a = Int(alphabet(s[3][1]))
    val = parse(T, s[4])
    return i, a, val
end
function parse_field_line_numerical(line, T)
    s = split(line, " ")
    i = parse(Int, s[2])
    a = parse(Int, s[3])
    val = parse(T, s[4])
    return i, a, val
end
function parse_coupling_line_symbolic(line, T, alphabet)
    s = split(line, " ")
    i = parse(Int, s[2])
    j = parse(Int, s[3])
    a = Int(alphabet(s[4][1]))
    b = Int(alphabet(s[5][1]))
    val = parse(T, s[6])
    return i, j, a, b, val
end
function parse_coupling_line_numerical(line, T)
    s = split(line, " ")
    i = parse(Int, s[2])
    j = parse(Int, s[3])
    a = parse(Int, s[4])
    b = parse(Int, s[5])
    val = parse(T, s[6])
    return i, j, a, b, val
end

#============================================================#
##################### Writing PottsGraph #####################
#============================================================#

"""
    write(file::AbstractString, g::PottsGraph; sigdigits)

Write parameters of `g` to `file` using the format `J i j a b value`.
"""
function write(
    file::AbstractString, g::PottsGraph; sigdigits=5, index_style=1, format=:numerical
)
    return write_graph_extended(file, g, sigdigits, index_style, format)
end

function write_graph_extended(
    file::AbstractString, g::PottsGraph, sigdigits, index_style, format
)
    @argcheck index_style == 0 || index_style == 1 "Got `index_style==`$(index_style)"
    @argcheck format in (:numerical, :symbolic) "Got format=$format"
    L, q = size(g)
    # symbolic format writes states as letters: the alphabet is fixed by q (aa or rna)
    alphabet = format == :symbolic ? _symbolic_alphabet_from_q(q) : nothing
    open(file, "w") do f
        for i in 1:L, j in (i + 1):L, a in 1:q, b in 1:q
            val = round(g.J[a, b, i, j]; sigdigits)
            if format == :symbolic
                a, b = alphabet.([a, b])
            end
            if index_style == 0 && format == :numerical
                write(f, "J $(i-1) $(j-1) $(a-1) $(b-1) $val\n")
            elseif index_style == 0 && format == :symbolic
                write(f, "J $(i-1) $(j-1) $(a) $(b) $val\n")
            elseif index_style == 1
                write(f, "J $i $j $a $b $val\n")
            else
                throw(
                    ArgumentError(
                        "Invalid arguments index_style=$(index_style) & format=$(format)"
                    ),
                )
            end
        end
        for i in 1:L, a in 1:q
            val = round(g.h[a, i]; sigdigits)
            if format == :symbolic
                a = alphabet(a)
            end
            if index_style == 0 && format == :numerical
                write(f, "h $(i-1) $(a-1) $val\n")
            elseif index_style == 0 && format == :symbolic
                write(f, "h $(i-1) $(a) $val\n")
            elseif index_style == 1
                write(f, "h $i $a $val\n")
            end
        end
    end
end

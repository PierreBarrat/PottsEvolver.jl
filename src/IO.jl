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
    ## Go through file twice: first to get L, the second to store parameters
    # Note that this function only works for amino acids --> q = 21
    L = 0
    q = 21
    min_idx = Inf
    index_style = 1
    for (n, line) in enumerate(eachline(file))
        if !is_valid_line(line, :symbolic)
            throw(ArgumentError("""
                Format problem with line $n in $file
                Expected format `J i j a b` or `h i a`, with amino acid symbols.\
                Instead $line"""))
        end
        if !isempty(line) && line[1] == 'h'
            i, a, val = parse_field_line_symbolic(line, T)
            L = max(L, i)
            min_idx = min(min_idx, i)
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
    @debug "Index style: $index_style"
    @debug L, q

    g = PottsGraph(L, q, T)
    for line in eachline(file)
        if line[1] == 'J'
            i, j, a, b, val = parse_coupling_line_symbolic(line, T)
            index_style == 0 && (i += 1; j += 1)
            g.J[a, b, i, j] = val
            g.J[b, a, j, i] = val
        elseif line[1] == 'h'
            i, a, val = parse_field_line_symbolic(line, T)
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

let
    aas = prod(symbols(aa_alphabet))
    patterns = Dict(
        :numerical => Regex.(["J [0-9]+ [0-9]+ [0-9]+ [0-9]+", "h [0-9]+ [0-9]+"]),
        :symbolic => Regex.(["J [0-9]+ [0-9]+ [$(aas)]+ [$(aas)]+", "h [0-9]+ [$(aas)]"]),
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
function parse_field_line_symbolic(line, T)
    s = split(line, " ")
    i = parse(Int, s[2])
    a = Int(aa_alphabet(s[3][1]))
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
function parse_coupling_line_symbolic(line, T)
    s = split(line, " ")
    i = parse(Int, s[2])
    j = parse(Int, s[3])
    a = Int(aa_alphabet(s[4][1]))
    b = Int(aa_alphabet(s[5][1]))
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
    open(file, "w") do f
        for i in 1:L, j in (i + 1):L, a in 1:q, b in 1:q
            val = round(g.J[a, b, i, j]; sigdigits)
            if format == :symbolic
                a, b = aa_alphabet.([a, b])
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
                a = aa_alphabet(a)
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

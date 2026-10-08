```@meta
CurrentModule = PottsEvolver
```

# CLI option reference

Every option of every `pottsevolver` subcommand. For what the options *mean*, and for the
schedule and parameter-file formats, see [Command line interface](@ref).

The tables below are generated from the CLI's own argument definitions when these docs are
built, so they cannot drift from `pottsevolver --help`.

Each table is complete for its subcommand: the options shared by both appear in both. Every
command also accepts `-h`/`--help`, which prints the same information and exits.

```@eval
using PottsEvolver, Markdown

settings = PottsEvolver.CLI.cli_settings()

# `--outdir`'s default is `pwd()`, captured when the settings are built: showing the literal
# path would print the directory the docs were built in
const HERE = pwd()

function option_names(f)
    if isempty(f.long_opt_name) && isempty(f.short_opt_name)
        return "`$(f.dest_name)`" # positional
    end
    names = vcat(["-$o" for o in f.short_opt_name], ["--$o" for o in f.long_opt_name])
    return join(("`$n`" for n in names), ", ")
end

function argument_of(f)
    f.action == :store_true && return "none (flag)"
    f.nargs == '*' && return "numbers"
    f.arg_type == Any && return "string"
    return "`$(f.arg_type)`"
end

function default_of(f)
    f.required && return "**required**"
    f.action == :store_true && return "off"
    d = f.default
    d == HERE && return "current directory"
    (isnothing(d) || (d isa AbstractVector && isempty(d))) && return "–"
    return "`$(d)`"
end

# Help strings are prose: their underscores and asterisks would read as Markdown emphasis, and
# a literal `|` would end the table cell. Escape only outside backtick spans, since a backslash
# inside one would be shown verbatim.
function cell(s)
    parts = split(string(s), '`')
    for i in eachindex(parts)
        isodd(i) || continue # even segments are inside a code span
        seg = replace(parts[i], "_" => "\\_", "*" => "\\*", "|" => "\\|")
        # a bare `--flag` in prose: render it as code, which also keeps Markdown from
        # turning the double hyphen into an en dash
        parts[i] = replace(seg, r"(?<![\w\\])--([A-Za-z][\w-]*)" => s"`--\1`")
    end
    return join(parts, '`')
end

function table(command)
    rows = ["| Option | Argument | Default | Description |", "|---|---|---|---|"]
    for f in settings[command].args_table.fields
        f.fake && continue
        push!(
            rows,
            "| $(option_names(f)) | $(argument_of(f)) | $(default_of(f)) | $(cell(f.help)) |",
        )
    end
    return join(rows, "\n")
end

# subcommands are the `%COMMAND%` entries of the top-level table: `constant` holds the name
# and `help` its one-line description, in declaration order
commands = [
    (f.constant, f.help) for f in settings.args_table.fields if f.action == :command_arg
]

sections = map(commands) do (command, description)
    join(["## `$command`", "", description, "", table(command)], "\n")
end

Markdown.parse(join(sections, "\n\n"))
```

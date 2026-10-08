```@meta
CurrentModule = PottsEvolver
```

# Command line interface

`pottsevolver` is a command line program that runs the sampling described in
[Sampling from Potts Models](@ref) without writing any Julia code: you give it a Potts model, a
parameter file, and it writes sequences to disk.
It covers the two main use cases, sampling along a Markov chain and sampling along the branches
of a phylogenetic tree.

## Installation

The CLI is distributed as a [Julia app](https://pkgdocs.julialang.org/v1/apps/): installing it
puts a `pottsevolver` executable on your system.

```julia
using Pkg
Pkg.Apps.add("PottsEvolver")
```

!!! warning
    The CLI requires PottsEvolver v0.2.0 or later. Earlier versions ship no app, and installing
    them gives you no executable. Until v0.2.0 is released, install from the repository:
    ```julia
    Pkg.Apps.add(; url="https://github.com/PierreBarrat/PottsEvolver.jl", rev="v0.2")
    ```

The executable is written to `~/.julia/bin`, which has to be on your `PATH`.
`Pkg` warns you at install time if it is not. To add it, put this in your `~/.bashrc`,
`~/.zshrc` or equivalent:

```sh
export PATH="$HOME/.julia/bin:$PATH"
```

Check that it worked:

```sh
pottsevolver --help
pottsevolver sample-chain --help
```

The two `--help` screens are the authoritative list of options: this page explains the ones
whose behaviour needs more than a line.

### Re-installing after a Julia upgrade

The launcher records the full path of the Julia binary that installed it, so it stops working
if that Julia is moved or removed — which is what happens when `juliaup` upgrades you to a new
version:

```
pottsevolver: exec: /home/you/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia: not found
```

To fix it, re-run the install command with the Julia you now want to use; the launcher is
regenerated against it.

```julia
using Pkg
Pkg.Apps.add("PottsEvolver")
```

`Pkg.Apps.status()` lists the apps you have installed, and `Pkg.Apps.rm("PottsEvolver")`
removes this one. `Pkg.Apps.update()` also regenerates the launcher, and in addition updates
the app itself.

## Quick start

Write the sampling parameters to a TOML file:

```toml
# params.toml
sampling_type = "discrete"
step_type = "gibbs"

[schedule]
spacing = 50     # time between samples
burnin = 100     # time before the first sample
```

Then sample 100 sequences from a Potts model stored in `model.dat` (any file that
[`read_graph`](@ref) accepts):

```sh
pottsevolver sample-chain model.dat -p params.toml -o results -M 100
```

This writes three files in `results/`:

```
results/model.fasta         the 100 sequences
results/model.params.toml   the parameters the run actually used
results/model.log           what happened
```

The output files are named after the model file; use `--prefix` to name them differently, and
`-o/--outdir` to put them elsewhere (it defaults to the current directory, and is created if
needed).

## Options

The tables below cover the options you are likely to need. `pottsevolver sample-chain --help`
and `pottsevolver sample-tree --help` remain the complete and authoritative list.

Both commands take these:

| Option | Default | Description |
|--------|---------|-------------|
| `model` | required | Path to the Potts model, as read by [`read_graph`](@ref) |
| `-p`, `--params` | required | Path to the TOML parameter file, see [The parameter file](@ref) |
| `-o`, `--outdir` | current directory | Where to write the output; created if it does not exist |
| `--prefix` | model file name, without extension | Base name of the output files, see [Output files](@ref) |
| `--init` | `random_aa` | Initial or root sequence, see [The initial sequence](@ref) |
| `--init-type` | inferred | Force the sequence type of an `--init` fasta: `aa`, `codon` or `rna` |
| `--burnin` | `5 * spacing`, or `0` on a tree | Time before the first sample; on a tree, time spent equilibrating the root |
| `--seed` | drawn, and recorded | Integer rng seed, see [Reproducibility](@ref) |
| `-v`, `--verbose` | `0` | Terminal verbosity, see [Logging and verbosity](@ref) |
| `--log-verbose` | `1` | Verbosity of the log file, same scale |
| `--compute-omega` | off | Compute `substitution_rate` when the parameter file omits it |
| `--translate` | off | Write amino acids instead of codons, for a `CodonSequence` run |

`sample-chain` additionally takes the sampling schedule, in one of three mutually exclusive
ways, see [`sample-chain`](@ref):

| Option | Default | Description |
|--------|---------|-------------|
| `-M`, `--n-samples` | `1` | Number of samples |
| `--spacing` | `0` | Time between samples |
| `--tvals` | — | Explicit schedule as space-separated numbers, the first being the burn-in |
| `--tvals-file` | — | The same, read from a file with one value per line |

`sample-tree` takes the tree and one flag, see [`sample-tree`](@ref):

| Option | Default | Description |
|--------|---------|-------------|
| `tree` | required | Path to the Newick tree file, given after `model` |
| `--internals` | off | Also write the internal node sequences and the sampled tree |

## The parameter file

`--params/-p` is required and points at a TOML file whose keys are the fields of
[`SamplingParameters`](@ref) — see its docstring for what each one means and what the defaults
are. Any field may be omitted, in which case `SamplingParameters`'s own default applies.

Fields that are Julia symbols are written as strings:

```toml
sampling_type = "continuous"   # or "discrete"
step_type = "glauber"          # :metropolis, :glauber, :sqrt, :gibbs, :mutsel
step_meaning = "accepted"
```

`branchlength_meaning` is a subtable:

```toml
[branchlength_meaning]
type = "step"      # or "sweep"
length = "exact"   # or "round", "poisson"
```

and `mutation_matrix` is given as an array **of rows**, so that `mutation_matrix[a][b]` is the
relative rate of mutating from state `a` to state `b`:

```toml
mutation_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
```

Invalid combinations are reported by `SamplingParameters` itself, so the error you get is the
same one the library would give.

## The initial sequence

`--init` chooses the sequence the chain starts from, or the sequence placed at the root of the
tree. It is either one of four keywords

```sh
pottsevolver sample-chain model.dat -p params.toml --init random_aa
```

(`random_aa`, `random_codon`, `random_rna`, `random_num` — a random sequence of the
corresponding type, with the length taken from the model), or a fasta file:

```sh
pottsevolver sample-chain model.dat -p params.toml --init root.fasta
pottsevolver sample-chain model.dat -p params.toml --init root.fasta,my_sequence
```

The first form requires the file to hold exactly one sequence. The second picks the record
whose name is exactly `my_sequence`, mirroring IQ-TREE `alisim`'s `--root-seq` syntax.

The type of sequence (amino acids, codons or RNA) is inferred from the characters in the file,
falling back to its length when the characters alone are ambiguous — `-ACGT` could be any of
the three alphabets. When inference cannot decide, it says so and you pick explicitly:

```sh
pottsevolver sample-chain model.dat -p params.toml --init root.fasta --init-type codon
```

`--init-type` takes `aa`, `codon` or `rna`. There is no fasta form for `NumSequence`, which has
no symbols to write: use `--init random_num`.

## `sample-chain`

Samples along a single Markov chain. The schedule — when samples are taken — is set in one of
three mutually exclusive ways.

### Equispaced samples

`-M/--n-samples`, `--spacing` and `--burnin` take `M` samples `spacing` apart, after `burnin`.
This is the common case:

```sh
pottsevolver sample-chain model.dat -p params.toml -M 100 --spacing 50 --burnin 100
```

The same schedule can live in the parameter file, which is handy when you reuse it:

```toml
[schedule]
spacing = 50
burnin = 100
n_samples = 100
```

Flags win over the file, entry by entry, so the line above with `-M 200` reuses the file's
`spacing` and `burnin` and takes 200 samples instead. `--burnin` defaults to `5 * spacing`, and
`--spacing` to 0 — which puts every sample at the same time, and warns, since no mcmc step
happens between them.

The file also accepts the same schedule said the other way round — first and last sample
instead of the step between them:

```toml
[schedule]
kind = "linear"   # the default, so it may be omitted
start = 0
stop = 5000
n_samples = 100
```

### Log-spaced schedules

Equispaced samples are a poor fit for watching a chain equilibrate, where you want the early
times in detail. `kind = "log"` takes the same three entries but spaces the times
geometrically:

```toml
[schedule]
kind = "log"
start = 1
stop = 100000
n_samples = 100
```

`start` must be positive, since a geometric schedule has no ratio from zero. This has no flag
equivalent: `--params` is required anyway, so an unusual schedule can simply live there.

Note that a schedule is stored as this handful of numbers rather than as its values, so asking
for 100 000 samples costs four lines in the file, not 100 000.

### Arbitrary schedules

`--tvals` gives the times explicitly instead. The first value is the burn-in, and later
values are *absolute* times, not gaps between samples:

```sh
pottsevolver sample-chain model.dat -p params.toml --tvals 100 150 200 400
```

`--tvals-file` reads the same thing from a file, one value per line, skipping blank lines and
lines starting with `#`:

```sh
pottsevolver sample-chain model.dat -p params.toml --tvals-file schedule.txt
```

Giving the times directly *and* describing them with `-M`/`--spacing`/`--burnin` is an error
rather than one silently winning: they are two ways of saying the same thing. Note that for
`sampling_type = "discrete"`, time is a number of MCMC steps, so an explicit schedule has to be
whole numbers.

An explicit schedule is the one case that cannot be stored as a few numbers. `--tvals` is
written into `<prefix>.params.toml` as it stands, while `--tvals-file` is copied next to the
other outputs as `<prefix>.tvals.txt` and referenced from the parameter file - so the file
stays readable however long the schedule is, and the output directory stays self-contained.

## `sample-tree`

Samples along the branches of a tree, which is given as a Newick file after the model:

```sh
pottsevolver sample-tree model.dat tree.nwk -p params.toml -o results
```

The root sequence comes from `--init`, and `--burnin` equilibrates it first — useful with a
random initial sequence, and left at 0 to place the given sequence at the root exactly as it
is. There is no schedule to give otherwise: sampling times come from the branch lengths of the
tree, see [Interpreting branch lengths](@ref).
See [Interpreting branch lengths](@ref) for how branch lengths map onto MCMC time.

By default only the leaf sequences are written. `--internals` additionally writes the internal
node sequences and the tree itself:

```sh
pottsevolver sample-tree model.dat tree.nwk -p params.toml --internals
```

The tree is written out because internal nodes that are unnamed in your input file are
automatically given names, and those are the names used to label the internal sequences. The
written tree is therefore the one the labels refer to.

## Output files

With `--prefix` unset, `<prefix>` below is the name of the model file without its extension.

| File | Written | Contents |
|------|---------|----------|
| `<prefix>.fasta` | always | the sample; for `sample-tree`, the leaf sequences |
| `<prefix>.params.toml` | always | the parameters of the run, plus a `[run]` table |
| `<prefix>.log` | always | the whole run, see [Logging and verbosity](@ref) |
| `<prefix>.substitutions.csv` | `track_substitutions = true` | one row per substitution: `sample,time,position,old,new` |
| `<prefix>.internals.fasta` | `sample-tree --internals` | internal node sequences |
| `<prefix>.nwk` | `sample-tree --internals` | the sampled tree, carrying the labels used above |

Sequences in the fasta files are labelled by *sampling time* for `sample-chain` — the values of
the schedule above — and by *node name* for `sample-tree`.

Existing files are overwritten, with a warning naming each one.

`track_substitutions` is only meaningful for continuous time sampling along a chain; the times
in `<prefix>.substitutions.csv` are counted from the start of the run, burn-in included.

## Reproducibility

Every run can be reproduced from its own output directory.

`--seed` sets the rng seed. If you do not give one, a seed is drawn for you — and either way it
is recorded in the `[run]` table of `<prefix>.params.toml`, together with the version of
PottsEvolver, the model file and the `--init` value:

```toml
[run]
init = "random_aa"
model = "model.dat"
seed = 1785483533
sequence_type = "AASequence{Int64}"
version = "0.2.0"
```

`load_params` ignores that table, so the file is a valid parameter file: you can hand a run's
own output straight back to `-p` and sample again with identical settings.

```sh
pottsevolver sample-chain model.dat -p results/model.params.toml -o rerun -M 100 --seed 1785483533
```

This matters most for continuous time sampling, where `substitution_rate` has to be set.
`--compute-omega` computes it for you when the parameter file leaves it out:

```sh
pottsevolver sample-chain model.dat -p continuous.toml --compute-omega -M 100
```

The computed value is written to `<prefix>.params.toml`, so a rerun from that file uses the
same rate instead of computing a new one. If the parameter file already sets
`substitution_rate`, the file wins and `--compute-omega` warns that it had no effect.

## Logging and verbosity

Every run writes `<prefix>.log`, covering the whole thing: reading the inputs, the sampling
itself, and the files written at the end. `--log-verbose` sets how much detail goes in it
(default `1`).

`-v/--verbose` does the same for the terminal, and defaults to `0`, which prints warnings and
errors only. Use `-v 1` to follow a run, including progress meters:

```sh
pottsevolver sample-chain model.dat -p params.toml -M 1000 -v 1
```

Both take the same scale: `-1` errors only, `0` warnings, `1` information, `2` and above debug.

## Exit codes

`pottsevolver` returns `0` on success and `1` for an error you can fix — an unreadable file,
invalid parameters, conflicting options. The message says what went wrong, and the log file
holds the details of a run that got far enough to start one. Anything else is a bug, and exits
with a Julia stacktrace.

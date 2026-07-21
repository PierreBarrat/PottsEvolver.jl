```@meta
CurrentModule = PottsEvolver
```

# PottsEvolver

Documentation for [PottsEvolver](https://github.com/PierreBarrat/PottsEvolver.jl).
In construction. 

```julia
using Pkg
Pkg.add("PottsEvolver")
using PottsEvolver
```

Simulate evolution of protein sequences (or other) using a Potts model. 
- Discrete or continuous time sampling.
- Sample single MCMC chains, or along branches of a tree.
- Sampling can take the genetic code into account. 

_Note_: the package relies on [`TreeTools.jl`](https://github.com/PierreBarrat/TreeTools.jl) to handle phylogenetic trees.

Sequences are read from and written to strings with `AASequence(::AbstractString)` and
`CodonSequence(::AbstractString)`, and samples are written to fasta with `write_fasta`.
`PottsEvolver` does not read fasta itself: use a dedicated package such as
[`FASTX.jl`](https://github.com/BioJulia/FASTX.jl) and convert.

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

_Note_: the package relies on two personal packages
- [`TreeTools.jl`](https://github.com/PierreBarrat/TreeTools.jl) to handle phylogenetic trees
- [`BioSequenceMappings.jl`](https://github.com/PierreBarrat/BioSequenceMappings.jl) to handle sequence alignments and to convert amino acids to digits.

_Remaining issues_: 
- sampling not reproducible when picking a random seed.
- use of `Alphabet` not clear: it's part of the `PottsGraph`, but only the default alphabets can really be used. 
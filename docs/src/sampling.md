```@setup sampling_1
using PottsEvolver
```

```@meta
DocTestSetup = quote
    using PottsEvolver
    L, q, M = 20, 21, 5; # lenth of the sequence, 21 amino acids, 5 samples
    g = PottsGraph(L, q, Float64; init=:rand); # random graph
    params = SamplingParameters(; sampling_type=:discrete, Teq=50, burnin=100)
    initial_sequence = AASequence(L)
    result = mcmc_sample(g, M, params; init=initial_sequence);
end
```
# Sampling from Potts Models

The [`mcmc_sample`](@ref) function allows you to sample sequences from a Potts model using Markov Chain Monte Carlo (MCMC) methods. This page provides an overview of the sampling process, including how to use the `mcmc_sample` function, its inputs, and outputs. 
The docstring of [`mcmc_sample`](@ref) provides extra details. 

## Overview

The `mcmc_sample` function supports different modes of sampling.
1. **Chain**: sample sequences along a Markov chain for a specified number of steps. 
2. **Tree**: sample sequences along the branches of a phylogenetic tree.
3. **Discrete**: the evolutionary time corresponds to a discrete number of mcmc steps. 
4. **Continuous**: the evolutionary time corresponds to a continuous time markov chain, as in many phylogenetics methods. 

It is of course possible to mix modes (1,2) with (3,4). 
The sampling process is controlled by the  structure [`SamplingParameters`](@ref), which allows you to customize the behavior of the MCMC algorithm.

---

## Chain Sampling

Start an evolutionary trajectory, _i.e._ a montecarlo chain, which is sampled at specified times. 
The code below samples `M=5` sequences, each separated by 50 discrete Gibbs steps. 
```@repl sampling_1
L, q, M = 20, 21, 5; # sequence of lenth 20, 21 amino acids, 5 samples
g = PottsGraph(L, q, Float64; init=:rand); # random graph
params = SamplingParameters(; sampling_type=:discrete, Teq=50, burnin=100)
initial_sequence = AASequence(L)
result = mcmc_sample(g, M, params; init=initial_sequence);
```

The output is a named tuple with the following fields:
- `sequences`: an alignment or vector of sampled sequences.
- `tvals`: a vector with the number of steps at each sample.
- `info`: information about the run.
- `params`: parameters used for the run as a `Dict`. 

```@repl sampling_1
result.sequences # a BioSequenceMappings.Alignment
write("/tmp/example.fasta", result.sequences) # can be written to fasta format
run(`cat /tmp/example.fasta`) # headers correspond to the time of sampling
```

The `burnin` parameter means that 100 steps are performed before taking the first sample in the chain. 
The first sequence of the chain is therefore different from `initial_sequence`: 
```jldoctest sampling_1
julia> AASequence(result.sequences[1]) != initial_sequence
true
```
If we want to start the chain from a specific sequence, we just set `burnin` to 0.
```jldoctest sampling_1
julia> params = SamplingParameters(; sampling_type=:discrete, Teq=50, burnin=0);

julia> result = mcmc_sample(g, M, params; init=initial_sequence);

julia> AASequence(result.sequences[1]) == initial_sequence
true
```

It is also possible to sample at non regularly spaced times by providing a vector of time values instead of the number `M` of sequences. 
```jldoctest sampling_1
julia> time_values = round.(Int, logrange(10, 1000, length=M)) # times at which the chain is sampled
5-element Vector{Int64}:
   10
   32
  100
  316
 1000

julia> result = mcmc_sample(g, time_values, params; init=initial_sequence);

julia> length(result.sequences)
5

julia> result.tvals == time_values
true

julia> result.sequences.names # sequence labels in the alignment correspond to sampling time
5-element Vector{String}:
 "10"
 "32"
 "100"
 "316"
 "1000"
```

---

## Tree Sampling

Sampling along the branches of the tree is also done using `mcmc_sample`. 
The input tree is provided either from a newick file, either as a `TreeTools.Tree` objet, see [TreeTools](https://github.com/PierreBarrat/TreeTools.jl). 
```@repl sampling_1
params = SamplingParameters(; sampling_type=:discrete, Teq=50, burnin=0);
result = mcmc_sample(g, "../../example/small_tree_integers.nwk", params; init=initial_sequence);
result.tree # the tree object, with sequences stored at each node
result.tree.root.data.seq == initial_sequence # the initial sequence is placed at the root
```
The return value contains alignments of leaf and internal node sequences. 
Sequences in the alignments are labeled according to the labels of the nodes in the input tree.
If nodes are not labeled in the input tree, they are automatically assigned labels which can be read in the output tree. 
```@repl sampling_1
result.internal_sequences
result.leaf_sequences
result.leaf_sequences.names
```
Alternatively, sequences can be stored in dictionaries indexed by node labels: 
```@repl sampling_1
result = mcmc_sample(
    g, "../../example/small_tree_integers.nwk", params; 
    init=initial_sequence, alignment_output=false
);
result.leaf_sequences
```

It is sometimes useful to perform sampling several times in a row. 
This is done by calling `mcmc_sample(graph, tree, M, params)`, where `M` is the number of repetitions. 
In this case, it is can be practical to have all sequences sampled at a given node grouped in a single alignment, instead of having a list of alignments that each correspond to one tree. 
This is achieved by calling `PottsEvolver.pernode_alignment` on the output of `mcmc_sample`:
```@repl sampling_1
M = 5 # five repeats
result = mcmc_sample(
    g, "../../example/small_tree_integers.nwk", M, params; init=initial_sequence
);
result = PottsEvolver.pernode_alignment(result);
result.sequences["1"] # all sequences sampled at node "1"
```

### Interpreting branch lengths

When sampling from a continuous time model, any real valued branch length is valid.
Things get more complicated when sampling from a discrete model, where evolutionary time is measured in number of MCMC steps. 
The [`BranchLengthMeaning`](@ref) type is used to facilitate this. 
It is provided as a field of `SamplingParameters`. For instance
```julia
SamplingParameters(; branchlength_meaning=BranchLengthMeaning(; type=:sweep, length=:exact))
```
will first multiply each branch in the tree by a factor `L` (length of the sequence), and expects the result to be exactly an integer. 

## Continuous time sampling

Discrete sampling is quite straightforward: the time corresponds to a discrete number of mcmc steps. 
In continuous time however, it is customary to scale time so that the expected number of mutations per site per unit of time is 1. 
This is done by setting the `average_substitution_rate` field of `SamplingParameters`, see [Generative continuous time model reveals epistatic signatures in protein evolution](https://doi.org/10.1101/2025.09.17.676821).
This average substituion rate can be computed either by sampling the potts model or by using an existing alignment of sequences. 

```@repl sampling_1
n_samples = 100 # compute the average transition rate from 100 samples
Teq = 10*L; # discrete steps between samples
step_type = :glauber;

Ω_1 = let 
    # first generate sample with discrete time, then use it to compute the average rate
    s0 = AASequence(L)
    params = SamplingParameters(; sampling_type=:discrete, Teq) 
    aln = mcmc_sample(g, M, params; init=s0).sequences
    write("/tmp/tmp.fasta", aln) # save the alignment for use below
    aln = map(AASequence, aln)
    # compute the average transition rate from the sampled alignment
    PottsEvolver.average_transition_rate(g, step_type, aln)
end

Ω_2 = let
    # Read the previously sampled alignment, stored in a fasta    
    # Ω_1 == Ω_2
    PottsEvolver.average_transition_rate(g, step_type, "/tmp/tmp.fasta")
end

Ω_3 = let
    # call the average_transition_rate function directly, which will sample the potts model
    s0 = AASequence(L)
    PottsEvolver.average_transition_rate(g, step_type, s0; n_samples, Teq)
end
```

Now, define `SamplingParameters(; sampling_type=:continuous, average_substitution_rate=Ω_1)` to have a time-normalized MCMC. 

---

## Additional Utilities

### Logging

You can control the verbosity of the sampling process using the `verbose` and `logfile` parameters:

```@repl sampling_1
result = mcmc_sample(g, 1000, params; init=:random_aa, verbose=1);
```

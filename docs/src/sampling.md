```@setup sampling_1
using PottsEvolver
```

```@meta
DocTestSetup = quote
    using PottsEvolver
    L, q, M = 20, 21, 5; # lenth of the sequence, 21 amino acids, 5 samples
    g = PottsGraph(L, q, Float64; init=:rand); # random graph
    params = SamplingParameters(; sampling_type=:discrete)
    initial_sequence = AASequence(L)
    result = mcmc_sample(g, M, params; spacing=50, burnin=100, init=initial_sequence);
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
`SamplingParameters` describes *how* a step is taken; *when* samples are taken - the spacing
between them, and the burn-in before the first - is given to `mcmc_sample` itself.
You are expected to provide a schedule: `spacing` defaults to 1, meaning a single mcmc step
between samples, which is seldom enough for successive sequences to be independent.
```@repl sampling_1
L, q, M = 20, 21, 5; # sequence of lenth 20, 21 amino acids, 5 samples
g = PottsGraph(L, q, Float64; init=:rand); # random graph
params = SamplingParameters(; sampling_type=:discrete) # how a step is taken
initial_sequence = AASequence(L)
# `spacing` and `burnin` say when samples are taken
result = mcmc_sample(g, M, params; spacing=50, burnin=100, init=initial_sequence);
```

The output is a named tuple with the following fields:
- `sequences`: a `SequenceSample`, or a vector of sampled sequences if `pack_output=false`.
- `tvals`: a vector with the number of steps at each sample.
- `info`: information about the run.
- `params`: parameters used for the run as a `Dict`. 

```@repl sampling_1
result.sequences # a SequenceSample: a matrix of integers plus one label per sequence
write_fasta("/tmp/example.fasta", result.sequences) # can be written to fasta format
run(`cat /tmp/example.fasta`) # headers correspond to the time of sampling
```

The `burnin` argument means that 100 steps are performed before taking the first sample in the chain. 
The first sequence of the chain is therefore different from `initial_sequence`: 
```jldoctest sampling_1
julia> AASequence(result.sequences[1]) != initial_sequence
true
```
If we want to start the chain from a specific sequence, we just set `burnin` to 0.
```jldoctest sampling_1
julia> result = mcmc_sample(g, M, params; spacing=50, burnin=0, init=initial_sequence);

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

julia> result.sequences.labels # sequence labels correspond to sampling time
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
params = SamplingParameters(; sampling_type=:discrete);
result = mcmc_sample(g, "../../example/small_tree_integers.nwk", params; init=initial_sequence);
result.tree # the tree object, with sequences stored at each node
result.tree.root.data.seq == initial_sequence # the initial sequence is placed at the root
```
The return value contains samples of leaf and internal node sequences. 
Sequences in the samples are labeled according to the labels of the nodes in the input tree.
If nodes are not labeled in the input tree, they are automatically assigned labels which can be read in the output tree. 
```@repl sampling_1
result.internal_sequences
result.leaf_sequences
result.leaf_sequences.labels
```
Alternatively, sequences can be stored in dictionaries indexed by node labels: 
```@repl sampling_1
result = mcmc_sample(
    g, "../../example/small_tree_integers.nwk", params; 
    init=initial_sequence, pack_output=false
);
result.leaf_sequences
```

It is sometimes useful to perform sampling several times in a row. 
This is done by calling `mcmc_sample(graph, tree, M, params)`, where `M` is the number of repetitions. 
In this case, it is can be practical to have all sequences sampled at a given node grouped in a single sample, instead of having a list of samples that each correspond to one tree. 
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
This is done by setting the `substitution_rate` field of `SamplingParameters`, see [Generative continuous time model reveals epistatic signatures in protein evolution](https://doi.org/10.1101/2025.09.17.676821).
That rate is usually written ``\Omega``, and has to be computed for the model at hand before
sampling from it: see [The substitution rate Ω](@ref).

## The substitution rate Ω

### What it is

```math
\Omega = \Big\langle \sum_{i=1}^{L} \sum_{a=1}^{q} Q(\mathbf{s} \to \mathbf{s}^{(i,a)}) \Big\rangle_{\mathbf{s}}
```

the total rate at which a sequence leaves its current state, averaged over sequences.
Passing it as `substitution_rate` rescales time so that one unit of time is one expected
substitution per site, which is what makes times — and the branch lengths of a tree —
comparable to those of standard phylogenetic models. Continuous sampling refuses to run
without it, since otherwise time has no defined unit.

The average has to be taken over the sequences the process actually visits, *i.e.* over its
equilibrium distribution. Two things can stand in for that: a sample of the model itself, or
an alignment of natural sequences, if you are willing to assume they are at equilibrium under
the model.

### Computing it

[`average_transition_rate`](@ref) does the averaging. Either hand it the sequences, or let it
draw them for you:

```@repl sampling_1
n_samples = 100; # average over 100 sequences
step_type = :glauber;

Ω_1 = let # sample the model first, then average over the result
    s0 = AASequence(L)
    params = SamplingParameters(; sampling_type=:discrete)
    S = mcmc_sample(g, n_samples, params; spacing=10*L, init=s0, pack_output=false).sequences
    average_transition_rate(g, step_type, S)
end

Ω_2 = let # same thing, with the sampling done for you
    average_transition_rate(g, step_type, AASequence(L); n_samples)
end
```

For natural sequences, read them yourself and pass the vector:

```julia
using FASTX
S = [AASequence(sequence(rec)) for rec in FASTAReader(open(fastafile))]
Ω = average_transition_rate(g, :glauber, S)
```

Then sample with `SamplingParameters(; sampling_type=:continuous, substitution_rate=Ω_1)`.

### With a mutation matrix

An asymmetric `mutation_matrix` ``\mu`` moves the equilibrium away from the Potts model and
towards mutation-selection balance. This breaks the second form above: it draws its sequences
with a *discrete* chain, and discrete sampling ignores `mutation_matrix` entirely, so it would
average ``\mu``-weighted rates over Potts-equilibrium sequences. Use the first form, and pass
the matrix to both halves — `average_transition_rate` defaults to `mutation_matrix=nothing`
and will otherwise return the wrong number without complaining.

That leaves an apparent chicken-and-egg problem: sampling the right equilibrium needs the
``\Omega`` you are trying to compute. It is milder than it looks, because ``\Omega`` only
rescales time. The substitution made at each step is drawn from the *unscaled* rates, so the
sequence of states visited — and hence the equilibrium being sampled — does not depend on the
value supplied. Any positive one gives a correctly distributed sample.

What it does change is the clock, and so how much sampling a requested time buys:

```math
n_{\rm subs} \approx T \, L \, \Omega_{\rm true} / \Omega_{\rm provisional}
```

An overestimated ``\Omega`` therefore does proportionally less work than asked, and the sample
may not be equilibrated at all. In practice:

1. Sample with a provisional ``\Omega`` — the Potts value above is a reasonable starting point,
   and deliberately scaling it *down* is the safe direction, since underestimating oversamples.
2. Judge equilibration in substitutions rather than in time, with
   `sum(x -> x.number_substitutions, result.info)` against however many substitutions per site
   you expect to need.
3. Compute ``\Omega`` from that sample, then sample again with it. Only now can `spacing` be
   chosen larger than the equilibration time in meaningful units, so the second sample is the
   one to use. Recomputing ``\Omega`` from it is a free consistency check: it should agree.

```julia
params = SamplingParameters(;
    sampling_type=:continuous, step_type=:glauber,
    mutation_matrix=μ, substitution_rate=Ω_potts/10, # deliberately low: oversample
    track_substitutions=true,
)
result = mcmc_sample(g, 500, params; spacing=10, burnin=1000, init=s0, pack_output=false)
sum(x -> x.number_substitutions, result.info) # enough substitutions to equilibrate?
Ω = average_transition_rate(g, :glauber, result.sequences; mutation_matrix=μ)
```

---

## Additional Utilities

### Logging

You can control the verbosity of the sampling process using the `verbose` and `logfile` parameters:

```@repl sampling_1
result = mcmc_sample(g, 1000, params; init=:random_aa, verbose=1);
```

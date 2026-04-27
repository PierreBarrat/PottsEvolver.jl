# Using the genetic code during sampling

Codon-based sampling is activated by providing a `CodonSequence` as the initial sequence (or using `init=:random_codon`).
It is based on the algorithm of [Emergent time scales of epistasis in protein evolution](https://www.pnas.org/doi/10.1073/pnas.2406807121).

The Potts model is still defined on amino acids, but mutations happen at the nucleotide level.
The energy of a codon sequence is

$$E(\mathbf{c}) = E_{\rm Potts}(\mathbf{a}) - \sum_i \log d(a_i)$$

where $\mathbf{a}$ is the amino acid sequence encoded by codons $\mathbf{c}$, and $d(a_i)$ is the number of synonymous codons for amino acid $a_i$.
The degeneracy terms ensures that the correct equilibrium distribution is recovered by downweighting amino acids that are degenerate in the genetic code.

## Discrete sampling

Uses the algorithm of [di Bari *et. al.* 2024](https://doi.org/10.1073/pnas.2406807121).
Each MCMC step is either an amino acid step or a gap step, controlled by `fraction_gap_step` in `SamplingParameters` (default: 0.9, meaning 90% of steps are gap steps). 

- **Amino acid step** (Gibbs): pick a sequence position and a base within its codon, then sample from the Boltzmann distribution over all coding codons accessible by a single nucleotide mutation at that base.
- **Gap step** (Metropolis): propose inserting or removing a gap codon at a random position.

## Continuous time sampling

Works for `step_type ∈ {:glauber, :metropolis, :sqrt}`.
The accessible neighbours of a codon are all coding codons reachable by a single nucleotide mutation, plus the gap codon.
`step_type=:gibbs` is **not** supported for codon sequences in continuous time.

## Example

```julia
params = SamplingParameters(; sampling_type=:discrete, Teq=100)
result = mcmc_sample(potts, 10, params; init=:random_codon, translate_output=true)
result.sequences  # amino acid alignment
```

Set `translate_output=false` to keep the output as codon sequences.

!!! note
    Full documentation for this feature is under construction.

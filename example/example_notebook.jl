### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ a9825e4c-92e6-11ef-168d-0f2ba863e26d
begin
    using Pkg
    Pkg.activate(".")
    using UnPack
    using PlutoUI
    using PottsEvolver
end

# ╔═╡ 992f54a2-c37a-47ce-9395-88dda801ecbc
begin
    using TreeTools # at this point it's just useful to have this
    using BioSequenceMappings
end

# ╔═╡ ab2d6cdb-b5be-438c-98c0-27bcb81f2f95
PlutoUI.TableOfContents()

# ╔═╡ 1c8c3060-3b4e-4407-a041-fc4f1c024332
md"# `PottsGraph`"

# ╔═╡ aed8bb56-d5e3-40ae-8133-11225486b62c
md"""
Potts models are contained in the `PottsGraph` type.
You can read them from a file with the following format
```
J i j a b value # coupling lines
...
h i a value # field lines
```
where `i,j` represent positions in a sequence and `a, b` amino acids (or other kind of states). 
"""

# ╔═╡ ccfa1d4e-0bc7-4b10-9315-66d404f57574
# Read from file
potts = read_graph("parameters_PF00076.dat") # alignment with L=69 sites 

# ╔═╡ aa766af7-c050-4928-8c90-b1af5cbd3527
md"""
The `PottsGraph` object essentially two fields: 
- `h`: matrix of size $(q \times L)$
- `J`: dimension 4 array of size $(q \times q \times L \times L)$
where $q$ is the number of possible states for a position (21 for amino acids and gaps) and $L$ is the length of the sequence.


"""

# ╔═╡ 907c9922-b45c-4973-8bf3-7dcc8219a61f
[
    size(potts.h), # q*L
    size(potts.J), # q*q*L*L
    size(potts), # named tuple (; L, q)
]

# ╔═╡ 95d11aca-45ef-4720-8bc8-0812319c6854
(; L, q) = size(potts)

# ╔═╡ 5d878adb-7813-4105-8be1-9e919467ac1a
md"""
Another field `alphabet` contains the mapping from integers (*e.g.* used to index `h` and `J`) and biological symbols (amino acids, codons...). It is an object of type `Alphabet`, coming from the `BioSequenceMappings.jl` package. 
Ideally, I would like this to allow the user to seamlessly change the meaning of the indices `a, b, ...`. However, it is not functional yet. I encourage the user to 
- not rely on the `alphabet` field for now
- if working with amino acids, ensure that the Potts model is consistent with the default alphabet defined by the string `"-ACDEFGHIKLMNPQRSTVWY"`. This is just assumed in some parts of the code. 
"""

# ╔═╡ 5045ef14-8547-43a2-a105-d4ec0962b04f
potts.alphabet 

# ╔═╡ 58197d9e-3b67-40db-9852-6451cf4b254e
md"""
One can use the alphabet as a function to map from symbols to integers or the reverse. 
"""

# ╔═╡ b8431fc3-d619-4ce3-bbf8-2e33e4963288
potts.alphabet(4) # convert integers to chars

# ╔═╡ 25e2f982-9a3c-4499-b2b9-dee74dc36951
potts.alphabet('D') # or the reverse

# ╔═╡ be5ed71f-c3ac-464c-839f-4829f36058cf
md"""
It is also possible to create graphs with arbitrary length and number of states. It can be useful for toy-models. 
The following cell creates a graph with `L=3` and `q=2` states, with random parameters. 
"""

# ╔═╡ 4d85b9a9-4c3c-4bcb-a956-9a8e31f5a525
toy_potts = PottsGraph(3, 2; init=:rand)

# ╔═╡ 6c962a9f-0324-4368-8c9d-2f6b033e1875
md"""
Writing the model to a file is easy
"""

# ╔═╡ 70164a3b-d708-4c02-b462-1187d3de3079
write("toy_potts.dat", toy_potts)

# ╔═╡ ed8559ec-2c67-4338-bcbb-20803700c24c
md"# Sampling"

# ╔═╡ b1408a8a-f62d-46f2-ad05-a611324793ad
md"""
The main sampling function is `mcmc_sample`. See `?mcmc_sample` for details. 
There are different ways to sample, namely: 
- continuous or discrete MCMC; 
- MCMC chain, or sample along branches of a tree;
- based on the genetic code or not. 

The first option is controlled by the `SamplingParameters` type. 
Whether to use the genetic code is decided based on the sequence used to initiate the sampling. If it is an `AASequence` or `NumSequence`, standard MCMC is used. If it is a `CodonSequence`, the genetic code will be used. In this latter case, the algorithm used is the one of *Emergent time scales of epistasis in protein evolution*, di Bari et. al., PNAS 2024.


"""

# ╔═╡ de0d517c-ca5c-4ef5-970d-b9ab71748c48
md"""
## Chains
"""

# ╔═╡ 5414be7f-8631-4885-bdd5-7fc35c0a33ca
md"""
As an example, the following will sample a chain of length `M=10` sequences, starting from a random amino acid sequence, and taking `Teq=100` *discrete* steps between each sample. Finally, the output is translated back to amino acids
"""

# ╔═╡ 259ed51a-75cc-4d02-a42c-b35054c4d1d0
result = let
    M = 10
    parameters = SamplingParameters(; Teq=100, burnin=0)
    mcmc_sample(potts, M, parameters; init=:random_aa)
end;

# ╔═╡ 7b174b8a-13bd-4534-bc0e-595149f9411f
md"""
The output contains the sample, the "time" of each sequence within the chain (*i.e.* how many steps away from the initial sequence), and other information.
"""

# ╔═╡ 7b455214-13e5-4646-a9c5-11e538475bca
begin
    @unpack sequences, tvals = result
    println("Time values: $tvals")
    sequences
end

# ╔═╡ f871ce6f-c80a-433e-a954-d34002d7b5c1
md"""
The sequences are in an `Alignment` object (from `BioSequenceMappings`). This essentially wraps the integer vectors in a `data` matrix, and the alphabet to convert them to biological symbols.

If the output was translated from codons to amino acids (as is the case here), the alphabet is automatically the same as the one from the input Potts model. Otherwise, it is `codon_alphabet`. 

We can use it to iterate over sequences, write the result to a fasta file, etc...
"""

# ╔═╡ 63fd2b58-4a71-4258-9eb4-96fb491bdd37
sequences[1] # the first sample (here also the initial sequence)

# ╔═╡ 42e624ed-b38d-4908-a4fa-7bc635c78846
size(sequences) # length L x 10 sequences

# ╔═╡ 97aa28b3-012b-4614-9e2a-7706d9d48402
md"""
It's also easy to write / read alignments. 
"""

# ╔═╡ 0017a217-7cd6-4974-b936-30f93631396e
let
    write("example_alignment.fasta", sequences)
    X = read_fasta("example_alignment.fasta")
    all(x -> x[1] == x[2], zip(X, sequences))
end

# ╔═╡ 72d67f50-de92-4795-bf09-c7a6e91114e5
md"""
Sometimes, it is easier to just have an array of sequences instead of an alignment. This is easily done:
"""

# ╔═╡ 0aab35e9-99d8-41ff-83d6-94cbc420468c
let
    M = 3
    parameters = SamplingParameters(; Teq=100, burnin=0)
    chain = mcmc_sample(
		potts, M, parameters; init=:random_aa, alignment_output=false
	)
	chain.sequences
end

# ╔═╡ 558f08ed-dfa0-4b18-8dc1-5c7c816b3534
md"""
### Parameters of the sampling
"""

# ╔═╡ 4a07a258-7a6f-4307-a759-ee9115020fea
md"""
Many properties of the sampling are controlled by a `SamplingParameters` object. The docstring `?SamplingParameters` should help. 
"""

# ╔═╡ edefb624-3dce-4891-9b39-0653ba39570d
md"""
### Continuous time sampling
"""

# ╔═╡ 8a4eaf2f-9ad9-45c3-9dd3-341ba4a7a688
md"""
To sample in continuous time we need to define a timescale based on the average substitution rate of the model, which must be computed. This is done by the `average_transition_rate` function. It samples from the input potts model and uses the resulting sequences to compute the average. This has to be done only once per model.
"""

# ╔═╡ 7a148263-9275-43ea-a1c1-08fa6a06c542
Ω_glauber = PottsEvolver.average_transition_rate(
	potts, 
	:glauber, # av. sub. rate depends on the step type
	AASequence(L) # random AA sequence to initialize the sampling 
)

# ╔═╡ c9fda537-4d44-4093-8878-783bd8b74d6e
md"""
Alternatively, the average transition rate can be computed using natural sequences.  
"""

# ╔═╡ 3bb57e16-361e-4825-86d5-87665cce7671
let
	parameters = SamplingParameters(
		sampling_type=:continuous,
		step_type=:glauber,
		substitution_rate=Ω_glauber,
		Teq=0.25, # time between samples
	)

	mcmc_sample(potts, 5, parameters; init=:random_aa).sequences
end

# ╔═╡ 9f7a28b5-6fad-45d3-a294-fe6217996206
md"""
Alternatively, a set of time values can be provided
"""

# ╔═╡ d1752c83-496b-415b-9eae-93dea204fdd7
let
	parameters = SamplingParameters(
		sampling_type=:continuous,
		step_type=:glauber,
		substitution_rate=Ω_glauber,
	)
	time_values = [0, 0.05, 0.1, 42.] # four samples
	mcmc_sample(potts, time_values, parameters; init=:random_aa).sequences
end

# ╔═╡ 8c2d4679-5c9a-46da-b053-2885fd4cf4a5
md"### Choosing the initial sequence"

# ╔═╡ 04fd2366-b9b5-4c45-ab0a-aee0d821fe1e
md"""
The initial sequence is provided to `mcmc_sample` either
- as a third positional argument
- as a kwarg `init`

In the first case, something of the type `AbstractSequence` should be used (see below). In the second case, the sequence is generated by calling `PottsEvolver.get_init_sequence(init, potts)`
the docstring of `PottsEvolver.get_init_sequence` will be helpful.

**Important**: the *type* of the initial sequence will determine the algorithm used for sampling:
- if it is a `CodonSequence`, then the MCMC is based on codons
- otherwise, it is just a standard MCMC
"""

# ╔═╡ aabe38f7-8f6e-4b66-949f-5da952193237
let
	# filter out integers that stand for stop codons in `codon_alphabet`. 
	q = 21
    not_stop = filter(x -> x != 13 && x != 5, 1:q)

	# Set of L states in [1:21]. Not a valid input because it's a Vector{Int}
    s_not_valid = rand(not_stop, L) 

    # Valid initial sequences
    aa_seq = AASequence(s_not_valid) # assumes the sequence is amino acids
    codon_seq = CodonSequence(s_not_valid) # assumes the sequence is codons
    num_seq = PottsEvolver.NumSequence(s_not_valid, q) # assumes the sequence has no biological meaning

    # Using this for sampling
    parameters = SamplingParameters(; Teq=100, burnin=0)
    mcmc_sample(potts, 10, aa_seq, parameters) # use "normal" amino acid based mcmc 
    mcmc_sample(potts, 10, codon_seq, parameters) # use codon based mcmc

    mcmc_sample(potts, 10, parameters; init=aa_seq) # this also works
    mcmc_sample(potts, 10, parameters; init=s_not_valid) # should work, but it's less clear
    mcmc_sample(potts, 10, parameters; init=:random_aa) # operates on aa, but starting from a random sequence
end

# ╔═╡ 65c82b80-cdfb-47f4-b398-989ec42ebbb3
let
    # Creating a random sequence is also easy
    rand_codon_seq = AASequence(L)
    # etc...
end

# ╔═╡ c5f720ac-6efb-4092-860e-4284cd510736
md"""
### Using the genetic code

If `init=:random_codon` or `init=sequence` where `sequence :: CodonSequence`, sampling will  use the genetic code.  
Note that this is **not** yet compatible with continuous time sampling (will implement in the future).
"""

# ╔═╡ 6576b8ba-762c-4e74-90b5-76f3725557d7
result_codons = let
    M = 10
    parameters = SamplingParameters(; Teq=100, burnin=0)
    mcmc_sample(potts, M, parameters; init=:random_codon, translate_output=true)
end;

# ╔═╡ 863b9af8-7d93-4633-929c-bdfe7e934512
md"""
The `translate_output` argument means that the sequences in the final alignment are translated back to amino acids. If you need the result as a sequence of codons, set `translate_output=false`. The default alphabet for codons is `codon_alphabet`.
"""

# ╔═╡ a1031609-9ee0-42da-9ab8-ea36fe2da040
codon_alphabet

# ╔═╡ 58c13678-8cc7-45f2-9720-6b9b47596ccd
md"### Controlling the output"

# ╔═╡ 0fd48be5-b07f-4376-a915-c224990e2c97
md"""
The following arguments control the output of `run_mcmc`: 
- `alignment_output` (default: `true`): return an `Alignment` containing the sequences. If `false`, will return a vector of sequences. 
  The advantage of the alignment format is mainly that it is easy to convert to fasta, and that some convenient functions of `BioSequenceMappings` operate on alignment. 
- `translate_output` (default: `true`): if a `CodonSequence` was provided as input, the output would normally also consist of codon sequences (or a codon based alignment).  This flag causes the output to be translated back to amino acids. 
"""

# ╔═╡ 578b2fea-f418-4f3c-8423-5be10412c884
md"""
If you do not translate the output, here is what happens: 
"""

# ╔═╡ d9e5a4a1-7fab-4140-94dd-a4943bb16a0e
sequences_codon = let
    M = 10
    parameters = SamplingParameters(; Teq=100, burnin=0)
    mcmc_sample(potts, M, parameters; init=:random_codon).sequences
end;

# ╔═╡ 2ddc7737-8ae2-4725-91dd-5076c3bbd350
sequences_codon[1] # integers larger than 21 --> codons

# ╔═╡ 9bb03d89-8156-468d-a8e1-b7c8f753f9af
# Here are three ways to map this to amino acids
let
    alphabet = sequences_codon.alphabet # same as potts.alphabet
    s0 = sequences_codon[1]
    @info "Integer codons" s0 # Vector{Integer} standing for codons
    @info "Codons" alphabet(s0) # vector of PottsEvolver.Codon
    @info "Integer AA" map(x -> genetic_code(x), s0) # Vector{Integer} standing for amino acids
    @info "Char AA" map(x -> genetic_code(alphabet(x)), s0) # vector of Char, for amino acids
end

# ╔═╡ 434db5aa-92d4-45dc-abda-30757305aa5e
md"## Sampling on a tree"

# ╔═╡ 25215879-eb19-46a6-bbdf-15bd96e18000
md"""
*Note*: The package uses `TreeTools.jl` to manipulate trees, please have a look at its documentation if you want to handle trees. The `read_tree` function is re-exported from `TreeTools`. 
"""

# ╔═╡ fe24fb58-bfc3-44e0-94c5-964bf0f71dc8
md"""
First, we read a small tree. It is a balanced binary tree with `8` leaves and randomly chosen branch lengths.
"""

# ╔═╡ 4a522c90-0f74-4b35-a2c8-f64726f4875d
tree = read_tree("small_tree.nwk")

# ╔═╡ 8ce3440d-353d-47f6-a669-dc7f65fd57d4
md"""
### Continuous sampling
"""

# ╔═╡ c463db7d-cb3c-407b-b8ec-03ce04ade53c
md"""
Sampling on this tree is easy: just call `mcmc_sample` with the tree object instead of the number of sequences `M`. 
If the continuous time version is used, sampling time on each branch is equal to the length of the branch. 
"""

# ╔═╡ 0ea5a273-b01f-4266-af55-212cf0b5870f
map(branch_length, traversal(tree; root=false)) # branch lengths

# ╔═╡ b1aea1ba-c13d-4205-8761-3c6bf557da05
result_tree = let
    parameters = SamplingParameters(
		sampling_type=:continuous,
		step_type=:glauber,
		substitution_rate=Ω_glauber,
	)
    mcmc_sample(potts, tree, parameters; init=:random_aa)
end

# ╔═╡ 91103fb3-896b-420b-ba28-51871f725a11
md"""
The output contains
- the tree,
- an alignment of leaf sequences,
- an alignment of internal nodes sequences.

The sequences are labeled using the label of tree nodes. 
"""

# ╔═╡ ac91aca7-1ea4-481f-9be0-be7cad82a892
result_tree.leaf_sequences

# ╔═╡ 6f661d71-923b-43c3-838c-ce7add792842
interesting_leaf = label(first(leaves(tree)))

# ╔═╡ 77ed2514-772c-45be-8b04-f4593d80a462
let
    # Find a sequence given a leaf label
    @info "Let's look at leaf $(interesting_leaf)"
    seq = find_sequence(interesting_leaf, result_tree.leaf_sequences)[2]
	@info "The numerical sequence: $seq"
    @info "Corresponding AA sequence: $(result_tree.leaf_sequences.alphabet(seq))"
end

# ╔═╡ e7e26419-29d1-4444-a248-5979c4ba755f
md"### Discrete sampling"

# ╔═╡ 1d30ae13-bfc2-4da1-9563-8fbe770b4734
md"""
If sampling is discrete, a convertion from branch lengths to integer number of steps needs to be done. 
"""

# ╔═╡ 53ad75d0-bb98-4195-b6ba-fb69e4b53b55


# ╔═╡ Cell order:
# ╠═a9825e4c-92e6-11ef-168d-0f2ba863e26d
# ╠═ab2d6cdb-b5be-438c-98c0-27bcb81f2f95
# ╟─1c8c3060-3b4e-4407-a041-fc4f1c024332
# ╟─aed8bb56-d5e3-40ae-8133-11225486b62c
# ╠═ccfa1d4e-0bc7-4b10-9315-66d404f57574
# ╟─aa766af7-c050-4928-8c90-b1af5cbd3527
# ╠═907c9922-b45c-4973-8bf3-7dcc8219a61f
# ╠═95d11aca-45ef-4720-8bc8-0812319c6854
# ╟─5d878adb-7813-4105-8be1-9e919467ac1a
# ╠═5045ef14-8547-43a2-a105-d4ec0962b04f
# ╟─58197d9e-3b67-40db-9852-6451cf4b254e
# ╠═b8431fc3-d619-4ce3-bbf8-2e33e4963288
# ╠═25e2f982-9a3c-4499-b2b9-dee74dc36951
# ╟─be5ed71f-c3ac-464c-839f-4829f36058cf
# ╠═4d85b9a9-4c3c-4bcb-a956-9a8e31f5a525
# ╟─6c962a9f-0324-4368-8c9d-2f6b033e1875
# ╠═70164a3b-d708-4c02-b462-1187d3de3079
# ╟─ed8559ec-2c67-4338-bcbb-20803700c24c
# ╟─b1408a8a-f62d-46f2-ad05-a611324793ad
# ╟─de0d517c-ca5c-4ef5-970d-b9ab71748c48
# ╟─5414be7f-8631-4885-bdd5-7fc35c0a33ca
# ╠═259ed51a-75cc-4d02-a42c-b35054c4d1d0
# ╟─7b174b8a-13bd-4534-bc0e-595149f9411f
# ╠═7b455214-13e5-4646-a9c5-11e538475bca
# ╟─f871ce6f-c80a-433e-a954-d34002d7b5c1
# ╠═63fd2b58-4a71-4258-9eb4-96fb491bdd37
# ╠═42e624ed-b38d-4908-a4fa-7bc635c78846
# ╟─97aa28b3-012b-4614-9e2a-7706d9d48402
# ╠═0017a217-7cd6-4974-b936-30f93631396e
# ╟─72d67f50-de92-4795-bf09-c7a6e91114e5
# ╠═0aab35e9-99d8-41ff-83d6-94cbc420468c
# ╟─558f08ed-dfa0-4b18-8dc1-5c7c816b3534
# ╟─4a07a258-7a6f-4307-a759-ee9115020fea
# ╟─edefb624-3dce-4891-9b39-0653ba39570d
# ╟─8a4eaf2f-9ad9-45c3-9dd3-341ba4a7a688
# ╠═7a148263-9275-43ea-a1c1-08fa6a06c542
# ╟─c9fda537-4d44-4093-8878-783bd8b74d6e
# ╠═3bb57e16-361e-4825-86d5-87665cce7671
# ╟─9f7a28b5-6fad-45d3-a294-fe6217996206
# ╠═d1752c83-496b-415b-9eae-93dea204fdd7
# ╟─8c2d4679-5c9a-46da-b053-2885fd4cf4a5
# ╟─04fd2366-b9b5-4c45-ab0a-aee0d821fe1e
# ╠═aabe38f7-8f6e-4b66-949f-5da952193237
# ╠═65c82b80-cdfb-47f4-b398-989ec42ebbb3
# ╟─c5f720ac-6efb-4092-860e-4284cd510736
# ╠═6576b8ba-762c-4e74-90b5-76f3725557d7
# ╟─863b9af8-7d93-4633-929c-bdfe7e934512
# ╠═a1031609-9ee0-42da-9ab8-ea36fe2da040
# ╟─58c13678-8cc7-45f2-9720-6b9b47596ccd
# ╟─0fd48be5-b07f-4376-a915-c224990e2c97
# ╟─578b2fea-f418-4f3c-8423-5be10412c884
# ╠═d9e5a4a1-7fab-4140-94dd-a4943bb16a0e
# ╠═2ddc7737-8ae2-4725-91dd-5076c3bbd350
# ╠═9bb03d89-8156-468d-a8e1-b7c8f753f9af
# ╟─434db5aa-92d4-45dc-abda-30757305aa5e
# ╟─25215879-eb19-46a6-bbdf-15bd96e18000
# ╟─fe24fb58-bfc3-44e0-94c5-964bf0f71dc8
# ╠═4a522c90-0f74-4b35-a2c8-f64726f4875d
# ╟─8ce3440d-353d-47f6-a669-dc7f65fd57d4
# ╟─c463db7d-cb3c-407b-b8ec-03ce04ade53c
# ╠═0ea5a273-b01f-4266-af55-212cf0b5870f
# ╠═b1aea1ba-c13d-4205-8761-3c6bf557da05
# ╟─91103fb3-896b-420b-ba28-51871f725a11
# ╠═ac91aca7-1ea4-481f-9be0-be7cad82a892
# ╠═992f54a2-c37a-47ce-9395-88dda801ecbc
# ╠═6f661d71-923b-43c3-838c-ce7add792842
# ╠═77ed2514-772c-45be-8b04-f4593d80a462
# ╟─e7e26419-29d1-4444-a248-5979c4ba755f
# ╟─1d30ae13-bfc2-4da1-9563-8fbe770b4734
# ╠═53ad75d0-bb98-4195-b6ba-fb69e4b53b55

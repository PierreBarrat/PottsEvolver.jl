module PottsEvolver

using ArgCheck
using Distributions
using Logging
using LoggingExtras
using PoissonRandom
using ProgressMeter
using Random
using StatsBase
using TreeTools
using UnPack

export read_tree # from TreeTools

import Base: ==, hash, isvalid
import Base: convert, copy, copy!, show, write
import Base: getindex, setindex!
import Base: iterate, length, eltype, size

export hamming

# Default types for numerical quantities
const IntType = Int64
const FloatType = Float64

include("codons.jl")
export codon_alphabet, aa_alphabet, rna_alphabet, symbols
export genetic_code
#! format: off
# public Q_AA, Q_RNA, Q_CODON
#! format: on

include("sequences.jl")
export AbstractSequence, AASequence, RNASequence, CodonSequence, NumSequence

include("misc.jl")

include("sample_output.jl")
export SequenceSample
#! format: off
# public sequence_type
#! format: on
#! format: off
# public translate
#! format: on
include("pottsgraph.jl")
export PottsGraph
export energy
#! format: off
# public set_gauge!
#! format: on

include("sampling_parameters.jl")
export BranchLengthMeaning, SamplingParameters

include("sampling_core.jl")
#! format: off
# public mcmc_steps!, steps_from_branchlength
#! format: on

include("sampling_continuous_core.jl")
export average_transition_rate

include("sampling_chain.jl")
#! format: off
# public mcmc_sample_chain
#! format: on

include("sample_tree.jl")
#! format: off
# public mcmc_sample_tree, pernode_alignment
#! format: on

include("sampling.jl")
export mcmc_sample
#! format: off
# public get_init_sequence
#! format: on

include("IO.jl")
export read_graph, read_potts_graph
export write_fasta

#=
- codons.jl: alphabets and genetic code
- sequences.jl: contain only a vector of Int (or two for CodonSequence).
  Conversion is done through alphabets
- IO.jl: for reading/writing Potts models, and for writing sequences to fasta.
  Reading fasta is left to the user: build sequences from strings with
  `AASequence(::AbstractString)` / `CodonSequence(::AbstractString)`.
=#

include("Parallel/Parallel.jl")

end

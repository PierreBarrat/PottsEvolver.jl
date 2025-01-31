module PottsEvolver

using ArgCheck
using BioSequenceMappings
using Logging
using LoggingExtras
using PoissonRandom
using ProgressMeter
using Random
using StatsBase
using TreeTools
using UnPack

export read_fasta, symbols # from BioSequenceMappings
export read_tree # from TreeTools

import Base: ==, hash
import Base: copy, show, write
import Base: getindex, setindex!
import Base: iterate, length, eltype, size

import BioSequenceMappings: Alignment, to_string, hamming
export Alignment, Alphabet

# Default types for numerical quantities
const IntType = UInt8
const FloatType = Float64

include("codons.jl")
export codon_alphabet, aa_alphabet, nt_alphabet
export bases, genetic_code

include("sequences.jl")
export AbstractSequence, AASequence, CodonSequence
#! format: off
# public translate
#! format: on
include("pottsgraph.jl")
export PottsGraph
export energy
#! format: off
# public set_gauge!
#! format: on

include("sampling_core.jl")
export BranchLengthMeaning, SamplingParameters
#! format: off
# public mcmc_steps!, steps_from_branchlength
#! format: on

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

#=
- codons.jl: alphabets and genetic code
- sequences.jl: contain only a vector of Int (or two for CodonSequence).
  Conversion is done through alphabets
- IO.jl: for reading Potts models -- alignment/sequence IO is done through BioSequenceMappings
=#

end

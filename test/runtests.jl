using Test
using PottsEvolver

using Accessors
using BioSequenceMappings
using Logging
using Random
using StatsBase
using TreeTools

@testset "PottsEvolver.jl" begin
    @testset "Codons" begin
        include("codons/test.jl")
    end

    @testset "PottsGraph" begin
        include("pottsgraph/test.jl")
    end

    @testset "Sequence" begin
        include("sequences/test.jl")
    end

    @testset "Sampling" begin
        include("sampling/test.jl")
    end

    @testset "Sampling on tree" begin
        include("sampling_tree/test.jl")
    end
end

@testset "Aqua.jl" begin
    include("Aqua.jl")
end

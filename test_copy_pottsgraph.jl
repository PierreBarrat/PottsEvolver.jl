# Unit tests for PottsGraph copy functionality
using Test
using Random
using PottsEvolver

@testset "PottsGraph copy tests" begin
    # Test 1: Basic copy functionality
    @testset "Basic copy" begin
        rng = MersenneTwister(123)
        L, q = 5, 21
        g_original = PottsGraph(L, q; init=:rand, rng=rng)
        
        g_copy = copy(g_original)
        
        # Test that they are equal
        @test g_copy.J == g_original.J
        @test g_copy.h == g_original.h
        @test g_copy.β == g_original.β
        @test g_copy.alphabet == g_original.alphabet
        
        # Test that they are different objects
        @test g_copy !== g_original
        @test g_copy.J !== g_original.J
        @test g_copy.h !== g_original.h
        @test g_copy.alphabet !== g_original.alphabet
    end
    
    # Test 2: Modification independence
    @testset "Modification independence" begin
        rng = MersenneTwister(456)
        L, q = 3, 4
        g_original = PottsGraph(L, q; init=:rand, rng=rng)
        g_copy = copy(g_original)
        
        # Store original values for comparison
        original_J_val = g_original.J[1, 2, 1, 2]
        original_h_val = g_original.h[1, 1]
        original_β_val = g_original.β
        
        # Modify the copy
        g_copy.J[1, 2, 1, 2] = 999.0
        g_copy.h[1, 1] = -555.0
        g_copy.β = 2.5
        
        # Original should be unchanged
        @test g_original.J[1, 2, 1, 2] == original_J_val
        @test g_original.h[1, 1] == original_h_val
        @test g_original.β == original_β_val
        
        # Copy should have the new values
        @test g_copy.J[1, 2, 1, 2] == 999.0
        @test g_copy.h[1, 1] == -555.0
        @test g_copy.β == 2.5
        
        # Test that they are different objects
        @test g_copy.J !== g_original.J
        @test g_copy.h !== g_original.h
    end
    
    # Test 3: Different initialization methods
    @testset "Different initialization methods" begin
        L, q = 4, 3
        
        # Null initialization
        g_null = PottsGraph(L, q; init=:null)
        g_null_copy = copy(g_null)
        @test g_null_copy.J == g_null.J
        @test g_null_copy.h == g_null.h
        
        # Random initialization
        rng = MersenneTwister(789)
        g_rand = PottsGraph(L, q; init=:rand, rng=rng)
        g_rand_copy = copy(g_rand)
        @test g_rand_copy.J == g_rand.J
        @test g_rand_copy.h == g_rand.h
    end
    
    # Test 4: Alphabet handling
    @testset "Alphabet handling" begin
        # With alphabet (q=21)
        L = 6
        g_with_alphabet = PottsGraph(L, 21; init=:null)
        g_with_alphabet_copy = copy(g_with_alphabet)
        @test g_with_alphabet_copy.alphabet == g_with_alphabet.alphabet
        @test g_with_alphabet_copy.alphabet !== g_with_alphabet.alphabet
        
        # Without alphabet (q≠21)
        g_no_alphabet = PottsGraph(L, 5; init=:null)
        g_no_alphabet_copy = copy(g_no_alphabet)
        @test g_no_alphabet_copy.alphabet === g_no_alphabet.alphabet === nothing
    end
    
    # Test 5: Type preservation
    @testset "Type preservation" begin
        rng = MersenneTwister(101112)
        L, q = 3, 3
        
        # Float64 (default)
        g_float64 = PottsGraph(L, q; init=:rand, rng=rng)
        g_float64_copy = copy(g_float64)
        @test eltype(g_float64_copy.J) === Float64
        @test eltype(g_float64_copy.h) === Float64
        @test typeof(g_float64_copy.β) === Float64
        
        # Float32
        g_float32 = PottsGraph(L, q, Float32; init=:rand, rng=rng)
        g_float32_copy = copy(g_float32)
        @test eltype(g_float32_copy.J) === Float32
        @test eltype(g_float32_copy.h) === Float32
        @test typeof(g_float32_copy.β) === Float32
    end
    
    # Test 6: Energy consistency
    @testset "Energy consistency" begin
        rng = MersenneTwister(131415)
        L, q = 4, 3
        g_original = PottsGraph(L, q; init=:rand, rng=rng)
        g_copy = copy(g_original)
        
        # Create a random sequence
        seq = rand(1:q, L)
        
        # Energies should be identical
        energy_original = PottsEvolver.energy(seq, g_original)
        energy_copy = PottsEvolver.energy(seq, g_copy)
        @test energy_original ≈ energy_copy
    end
end

println("All PottsGraph copy tests passed!")
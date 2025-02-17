@testset "Energy differences" begin
    L = 5
    @testset "AA sequences" begin
        q = 21
        g = PottsGraph(L, q; init=:rand)

        a = AASequence(L)
        a.seq[1] = 1
        ref_ΔE = PottsEvolver.compute_energy_differences(a, g)

        b = copy(a)
        i = 1
        x = rand(2:q)
        b.seq[i] = x
        ΔE_true = PottsEvolver.compute_energy_differences(b, g)

        ΔE_buffer = similar(ref_ΔE)
        ΔE_test = PottsEvolver.compute_energy_differences!(ΔE_buffer, ref_ΔE, a, i, x, g)

        @test all(ΔE_test .≈ ΔE_true)
    end

    @testset "Codon sequences" begin
        q = length(codon_alphabet)
        g = PottsGraph(L, q; init=:rand)

        a = CodonSequence(L)
        a.seq[1] = 1
        a.aaseq[1] = 1
        ref_ΔE = PottsEvolver.compute_energy_differences(a, g)

        b = copy(a)
        i = 1
        x = rand(2:q)
        b.seq[i] = x
        b.aaseq[i] = genetic_code(x)
        ΔE_true = PottsEvolver.compute_energy_differences(b, g)

        ΔE_buffer = copy(ref_ΔE)
        ΔE_test = PottsEvolver.compute_energy_differences!(ΔE_buffer, ref_ΔE, a, i, x, g)

        @test all(ΔE_test .≈ ΔE_true)
    end

    @testset "_delta_energy" begin
        L = 3
        q = 21
        g = PottsGraph(L, q; init=:rand)

        @testset "AASequence" begin
            # Create two sequences differing at position 1
            seq = AASequence([1, 2, 3])
            refseq = AASequence([4, 2, 3])
            
            # Calculate energy difference at position 1
            dE = PottsEvolver._delta_energy(seq, refseq, 1, g)
            
            # Verify by calculating full energies
            E1 = energy(seq, g)
            E2 = energy(refseq, g)
            @test dE ≈ E1 - E2
            
            # Test that sequences differing at multiple positions throw an error
            seq_invalid = AASequence([1, 5, 3])
            @test_throws ArgumentError PottsEvolver._delta_energy(seq_invalid, refseq, 1, g)
        end

        @testset "CodonSequence" begin
            # Create two sequences differing at position 1
            seq = CodonSequence([1, 2, 3])
            refseq = copy(seq)
            refseq.seq[1] = 4  # Different codon at position 1
            refseq.aaseq[1] = genetic_code(4)
            
            # Calculate energy difference at position 1
            dE = PottsEvolver._delta_energy(seq, refseq, 1, g)
            
            # Verify by calculating full energies
            E1 = energy(seq, g)
            E2 = energy(refseq, g)
            @test dE ≈ E1 - E2
            
            # Test that sequences differing at multiple positions throw an error
            seq_invalid = copy(seq)
            seq_invalid.seq[2] = 4
            seq_invalid.aaseq[2] = genetic_code(4)
            @test_throws ArgumentError PottsEvolver._delta_energy(seq_invalid, refseq, 1, g)
        end
    end
end

@testset "transition rates" begin
    L = 3
    q = 21
    g = PottsGraph(L, q; init=:rand)

    @testset "CodonSequence" begin
        # Testing all three ways to compute transition rates
        # - working on a pre-allocated CTMCState
        # - from scratch, but still using the buffer from the CTMCState
        # - from an uninitialized CTMCState 
        refseq = CodonSequence([1, 2, 3])
        state = PottsEvolver.CTMCState(refseq)

        Q1, R1 = PottsEvolver.transition_rates!(state, g, :glauber)
        Q2, R2 = PottsEvolver.transition_rates!(state, g, :glauber; from_scratch=true)
        state.previous_seq = nothing
        Q3, R3 = PottsEvolver.transition_rates!(state, g, :glauber)

        @test Q1 ≈ Q2
        @test Q2 ≈ Q3
        @test R1 ≈ R2
        @test R2 ≈ R3
    end
end
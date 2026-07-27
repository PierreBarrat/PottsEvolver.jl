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
        q = PottsEvolver.Q_CODON
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

@testset "mutation_matrix" begin
    L, q = 5, 5
    g = PottsGraph(L, q; init=:rand)
    seq = NumSequence(rand(1:q, L), q)

    @testset "rate mechanic" begin
        # μ multiplies the rate (i, a) -> (i, b) by μ[a, b], with `a` the current state at i
        μ = rand(q, q)
        Q0 = PottsEvolver.transition_rates(seq, g, :glauber)
        Qμ = PottsEvolver.transition_rates(seq, g, :glauber; mutation_matrix=μ)
        cur = PottsEvolver.sequence(seq)
        @test all(Qμ[b, i] ≈ Q0[b, i] * μ[cur[i], b] for b in 1:q, i in 1:L)
    end

    @testset "μ of ones is a no-op" begin
        μ = ones(q, q)
        Q0 = PottsEvolver.transition_rates(seq, g, :glauber)
        Qμ = PottsEvolver.transition_rates(seq, g, :glauber; mutation_matrix=μ)
        @test Qμ ≈ Q0
    end

    @testset "wrong size rejected at sampling time" begin
        params = SamplingParameters(;
            sampling_type=:continuous, step_type=:glauber, Teq=1.0, substitution_rate=1.0,
            mutation_matrix=ones(q + 1, q + 1),
        )
        @test_throws ArgumentError mcmc_sample(g, 3, params; init=:random_num)
    end

    @testset "codon sequences rejected" begin
        gc = PottsGraph(L, 21; init=:null)
        params = SamplingParameters(;
            sampling_type=:continuous, step_type=:glauber, Teq=1.0, substitution_rate=1.0,
            mutation_matrix=ones(PottsEvolver.Q_CODON, PottsEvolver.Q_CODON),
        )
        @test_throws ArgumentError mcmc_sample(gc, 3, params; init=:random_codon)
    end

    @testset "sampling with μ is reproducible" begin
        μ = rand(q, q)
        params = SamplingParameters(;
            sampling_type=:continuous, step_type=:glauber, Teq=1.0, substitution_rate=1.0,
            mutation_matrix=μ,
        )
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, 10, params; rng, init=:random_num).sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, 10, params; rng, init=:random_num).sequences
        for i in 1:length(aln_1)
            @test aln_1[i] == aln_2[i]
        end
    end
end

@testset "Reproducibility" begin
    @testset "Continuous - Glauber" begin
        L, q = 10, 21
        g = PottsGraph(L, q; init=:null)
        params = SamplingParameters(;
            sampling_type=:continuous, step_type=:glauber, Teq=1.0
        )
        init = :random_aa
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, 10, params; rng, init).sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, 10, params; rng, init).sequences
        for i in 1:length(aln_1)
            @test aln_1[i] == aln_2[i]
        end
    end
    @testset "Continuous - Gibbs" begin
        L, q = 10, 21
        g = PottsGraph(L, q; init=:null)
        params = SamplingParameters(; sampling_type=:continuous, step_type=:gibbs, Teq=1.0)
        init = :random_aa
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, 10, params; rng, init).sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, 10, params; rng, init).sequences
        for i in 1:length(aln_1)
            @test aln_1[i] == aln_2[i]
        end
    end
    @testset "Continuous - Sqrt - Codon" begin
        L, q = 10, 21
        g = PottsGraph(L, q; init=:null)
        params = SamplingParameters(; sampling_type=:continuous, step_type=:sqrt, Teq=1.0)
        init = :random_codon
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, 10, params; rng, init).sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, 10, params; rng, init).sequences
        for i in 1:length(aln_1)
            @test aln_1[i] == aln_2[i]
        end
    end
    @testset "Continuous - metropolis - Num" begin
        L, q = 10, 21
        g = PottsGraph(L, q; init=:null)
        params = SamplingParameters(;
            sampling_type=:continuous, step_type=:metropolis, Teq=1.0
        )
        init = :random_num
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, 10, params; rng, init).sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, 10, params; rng, init).sequences
        for i in 1:length(aln_1)
            @test aln_1[i] == aln_2[i]
        end
    end

    @testset "Continuous - glauber - Codon - Tree" begin
        L, q = 10, 21
        g = PottsGraph(L, q; init=:null)
        tree = TreeTools.Generate.balanced_binary_tree(8, 1.0)
        params = SamplingParameters(;
            sampling_type=:continuous, step_type=:glauber, Teq=1.0
        )
        init = :random_codon
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, tree, params; rng, init).leaf_sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, tree, params; rng, init).leaf_sequences

        for i in 1:length(aln_1)
            @test aln_1[i] == aln_2[i]
        end
    end
end

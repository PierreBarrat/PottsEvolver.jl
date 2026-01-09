@testset "SamplingParameters" begin
    @test SamplingParameters(; Teq=5) isa Any
    @test_throws ArgumentError SamplingParameters(; Teq=5, step_type=:dubstep)
    @test_throws ArgumentError SamplingParameters(; Teq=5, step_meaning=:olive_tree)

    @testset "Construction" begin
        sampling_type = :discrete
        @test_throws ArgumentError SamplingParameters(; sampling_type, Teq=5.2) # discrete needs integer times
        @test_throws ArgumentError SamplingParameters(; sampling_type, burnin=1.5)  # discrete needs integer times
        @test_throws ArgumentError SamplingParameters(; sampling_type, Teq=5.2, burnin=1.5)  # type of Teq determines type of burnin
        @test SamplingParameters(; sampling_type, Teq=5.0, burnin=1.0) isa
            SamplingParameters{Int}

        sampling_type = :continuous # everything should work and give SamplingParameters{Float64}
        @test SamplingParameters(; sampling_type, Teq=5.2) isa SamplingParameters{Float64}
        @test SamplingParameters(; sampling_type, burnin=1.5) isa
            SamplingParameters{Float64}
        @test SamplingParameters(; sampling_type, Teq=5.2, burnin=1.5) isa
            SamplingParameters{Float64}
        @test SamplingParameters(; sampling_type, Teq=50) isa SamplingParameters{Float64}
        @test SamplingParameters(; sampling_type, Teq=5, burnin=1.5) isa
            SamplingParameters{Float64}
    end
end

@testset "Accepted steps" begin
    L, q, M = (5, 3, 100)
    g = PottsGraph(L, q; init=:rand)

    # the difference between two samples should always be exactly one.
    params = SamplingParameters(; step_meaning=:changed, Teq=1, burnin=3)
    S, tvals, _ = mcmc_sample(g, M, params)
    @test unique(map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)) == [1]
    @test tvals == collect(3 .+ (0:((M - 1))))

    # difference from init should be 1, then no changes
    params = SamplingParameters(; step_meaning=:changed, Teq=0, burnin=1)
    s0 = NumSequence(L, q)
    S, tvals, _ = mcmc_sample(g, M, s0, params)
    @test tvals == 1 .+ zeros(Int, M)
    @test allequal(S)
    @test hamming(S[1], s0.seq; normalize=false) == 1

    # the difference between two samples should be strictly less than one
    # (!! this is a statistically true test, but could be false!)
    params = SamplingParameters(; step_meaning=:changed, Teq=1, burnin=0)
    params = @set params.step_meaning = :accepted
    S, tvals, _ = mcmc_sample(g, M, params)
    H = map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)
    acc_ratio_1 = sum(H) / length(H)
    @test acc_ratio_1 < 1
    @test tvals == collect(0 .+ (0:(M - 1)))

    # The difference between two samples should be 0
    params = SamplingParameters() # Teq=0, burnin=0
    S, _ = mcmc_sample(g, M, params)
    @test all(==(0), map(i -> hamming(S[i - 1], S[i]), 2:M))
end

@testset "Output values" begin
    L, q, M = (4, 21, 2)
    g = PottsGraph(L, q; init=:rand)
    @test g.alphabet == aa_alphabet

    params = SamplingParameters(; Teq=1)

    ## Num sequence
    # because asked explicitely
    S, _ = mcmc_sample(g, M, params; init=:random_num, alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.NumSequence}

    # because g has no alphabet
    g_noalphabet = PottsGraph(L, q; init=:rand, alphabet=nothing)
    S, _ = mcmc_sample(g_noalphabet, M, params; init=[1, 2, 3, 4], alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.NumSequence}

    # because g has an alphabet that is not the default aa
    alphabet = Alphabet("ACDEFGHIKLMNPQRSTVWY-")
    g_strangealphabet = PottsGraph(L, q; init=:rand, alphabet)
    S, _ = mcmc_sample(
        g_strangealphabet, M, params; init=[1, 2, 3, 4], alignment_output=false
    )
    @test S isa AbstractVector{<:PottsEvolver.NumSequence}

    ## AA sequence
    # because asked explicitely
    S, _ = mcmc_sample(g, M, params; init=:random_aa, alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.AASequence}

    # because g has aa_alphabet and init vector has elements <= 21
    S, _ = mcmc_sample(g, M, params; init=[1, 2, 3, 21], alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.AASequence}

    ## Codon sequence
    # because asked explicitely
    S, _ = mcmc_sample(g, M, params; init=:random_codon, alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.CodonSequence}

    # because g has aa_alphabet and init vector has elements > 21
    S, _ = mcmc_sample(g, M, params; init=[1, 2, 3, 22], alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.CodonSequence}
end

@testset "Continuous/discrete dispatch" begin
    L, q, M = (4, 21, 2)
    g = PottsGraph(L, q; init=:rand)

    # Test discrete case - Integer Teq and :discrete sampling_type
    params_discrete = SamplingParameters(; Teq=Int(5), sampling_type=:discrete)
    S_discrete = mcmc_sample(g, M, params_discrete; alignment_output=false).sequences
    @test S_discrete isa AbstractVector{<:PottsEvolver.NumSequence}

    # Test continuous case - Float Teq and :continuous sampling_type
    params_continuous = SamplingParameters(;
        Teq=5.0, sampling_type=:continuous, step_type=:glauber
    )
    S_continuous =
        mcmc_sample(
            g, M, params_continuous; alignment_output=false, init=:random_aa
        ).sequences
    @test S_continuous isa AbstractVector{<:PottsEvolver.AASequence}
    #= NEED TO CHECK TESTS BELOW  =#
    # Test mixed case 1 - Integer Teq but :continuous sampling_type
    # Should convert to continuous
    params_mixed1 = SamplingParameters(;
        Teq=Int(5), sampling_type=:continuous, step_type=:metropolis
    )
    S_mixed1 =
        mcmc_sample(
            g, M, params_mixed1; alignment_output=false, init=:random_codon
        ).sequences
    @test S_mixed1 isa AbstractVector{<:PottsEvolver.CodonSequence}

    # Test mixed case 2 - Float Teq (convertible to Int) and :discrete sampling_type
    # Should work by converting to Int
    params_mixed2 = SamplingParameters(; Teq=5.0, sampling_type=:discrete)
    S_mixed2, _ = mcmc_sample(g, M, params_mixed2; alignment_output=false, init=:random_num)
    @test S_mixed2 isa AbstractVector{<:PottsEvolver.NumSequence}
end

@testset "Sampling with time values" begin
    L, q = 10, 5
    g = PottsGraph(L, q)

    @testset "Errors" begin
        params_discrete = SamplingParameters(; sampling_type=:discrete)
        params_continuous = SamplingParameters(; sampling_type=:continuous)

        tvals = [0, 2, 1]
        @test_throws ArgumentError mcmc_sample(g, tvals, params_discrete)
        @test_throws MethodError mcmc_sample(g, tvals, params_continuous) # tvals is Ints
        @test_throws ArgumentError mcmc_sample(g, Float64.(tvals), params_continuous)

        tvals = [-1, 0, 1]
        @test_throws ArgumentError mcmc_sample(g, tvals, params_discrete)
        @test_throws MethodError mcmc_sample(g, tvals, params_continuous)
        @test_throws ArgumentError mcmc_sample(g, Float64.(tvals), params_continuous)

        tvals = [0.0, 1.0, 2.0]
        @test_throws MethodError mcmc_sample(g, tvals, params_discrete)
    end

    # test in discrete case
    @testset "Discrete" begin
        sampling_type = :discrete

        time_steps = 0:20
        M = length(time_steps)

        # the difference between two samples should always be exactly one.
        params = SamplingParameters(; sampling_type, step_meaning=:changed, Teq=4, burnin=0)
        S, tvals, _ = mcmc_sample(g, time_steps, params)
        @test length(S) == M
        @test unique(map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)) == [1]
        @test tvals == time_steps

        # the difference between two samples should be strictly less than one
        # (!! this is a statistically true test, but could be false!)
        params = SamplingParameters(; step_meaning=:accepted, sampling_type)
        S, tvals, _ = mcmc_sample(g, time_steps, params)
        @test length(S) == M
        @test tvals == time_steps
        H = map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)
        acc_ratio_1 = sum(H) / length(H)
        @test acc_ratio_1 < 1

        # The difference between two samples should be 0
        time_steps = [0, 0, 0]
        M = 3
        params = SamplingParameters(; sampling_type)
        S, _ = mcmc_sample(g, time_steps, params)
        @test length(S) == M
        @test all(==(0), map(i -> hamming(S[i - 1], S[i]), 2:M))

        # burnin and changed
        params = SamplingParameters(; sampling_type, step_meaning=:changed)
        time_steps = [1]
        s0 = NumSequence(L, q)
        S, tvals, _ = mcmc_sample(g, time_steps, s0, params)
        @test hamming(S[1], s0.seq; normalize=false) == 1
    end

    @testset "Continuous" begin
        sampling_type = :continuous
        params = SamplingParameters(; sampling_type, step_meaning=:changed, Teq=4, burnin=0)

        time_steps = range(0, 3; length=5)
        S, tvals, _ = mcmc_sample(g, time_steps, params)
        @test length(S) == length(time_steps)
        @test tvals == time_steps

        time_steps = 3.0 .+ zeros(Float64, 5)
        s0 = NumSequence(L, q)
        S, tvals, _ = mcmc_sample(g, time_steps, s0, params)
        @test length(S) == length(time_steps)
        @test tvals == time_steps
        @test allequal(S)
        @test S[1] != s0.seq
    end
end

@testset "Reproducibility" begin
    L, q = 10, 5
    g = PottsGraph(L, q; init=:null)
    s0 = NumSequence(L, q)
    params = SamplingParameters(; sampling_type=:discrete, Teq=10)

    rng = Random.seed!(Xoshiro(123), 42)
    aln_1 = mcmc_sample(g, 2, params; init=s0, rng).sequences

    rng = Random.seed!(Xoshiro(123), 42)
    aln_1 = mcmc_sample(g, 2, params; init=s0, rng).sequences

    @test_broken aln_1[1] == aln_2[1] && aln_1[2] == aln_2[2]
end

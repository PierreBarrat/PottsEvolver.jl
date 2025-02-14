@testset "SamplingParameters" begin
    @test SamplingParameters(; Teq=5) isa Any
    @test_throws ArgumentError SamplingParameters(; Teq=5, step_type=:dubstep)
    @test_throws ArgumentError SamplingParameters(; Teq=5, step_meaning=:olive_tree)

    @testset "Conversion" begin
        p_int = SamplingParameters(; Teq=5)
        p_float = SamplingParameters(; Teq=5., sampling_type=:continuous, step_type=:sqrt)
        p_float_2 = SamplingParameters(; Teq=5.5, sampling_type=:continuous, step_type=:sqrt)

        # Test no-op cases
        p_int_converted = convert(SamplingParameters{Int}, p_int)
        p_float_converted = convert(SamplingParameters{Float64}, p_float)
        for p in propertynames(p_int)
            @test getproperty(p_int_converted, p) == getproperty(p_int, p)
            @test getproperty(p_float_converted, p) == getproperty(p_float, p)
        end

        # Test actual conversion
        p_int_converted = convert(SamplingParameters{Float64}, p_int)
        p_float_converted = convert(SamplingParameters{Int}, p_float)
        @test p_int_converted isa SamplingParameters{Float64}
        @test p_float_converted isa SamplingParameters{Int}
        for p in propertynames(p_int)
            @test getproperty(p_int_converted, p) == getproperty(p_int, p)
            @test getproperty(p_float_converted, p) == getproperty(p_float, p)
        end        

        # Test fail
        @test_throws InexactError convert(SamplingParameters{Int}, p_float_2)
    end
end

@testset "Accepted steps" begin
    L, q, M = (5, 3, 500)
    g = PottsGraph(L, q; init=:rand)

    # the difference between two samples should always be exactly one.
    params = SamplingParameters(; step_meaning=:changed, Teq=1)
    S, _ = mcmc_sample(g, M, params)
    @test unique(map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)) == [1]

    # the difference between two samples should be less than one (!! this is a statistically true test, but could be false!)
    params = @set params.step_meaning = :accepted
    S, _ = mcmc_sample(g, M, params)
    H = map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)
    acc_ratio_1 = sum(H) / length(H)
    @test acc_ratio_1 < 1

    # The difference between two samples should be 0
    params = SamplingParameters() # Teq=0, burnin=0
    S, _ = mcmc_sample(g, M, params)
    @test all(==(0), map(i -> hamming(S[i-1], S[i]),2:M))
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
    params_continuous = SamplingParameters(; Teq=5.0, sampling_type=:continuous, step_type=:glauber)
    S_continuous = mcmc_sample(
        g, M, params_continuous; alignment_output=false, init=:random_aa
    ).sequences
    @test S_continuous isa AbstractVector{<:PottsEvolver.AASequence}
    #= NEED TO CHECK TESTS BELOW  =#
    # Test mixed case 1 - Integer Teq but :continuous sampling_type
    # Should convert to continuous
    params_mixed1 = SamplingParameters(; Teq=Int(5), sampling_type=:continuous, step_type=:metropolis)
    S_mixed1 = mcmc_sample(
        g, M, params_mixed1; alignment_output=false, init=:random_codon
    ).sequences
    @test S_mixed1 isa AbstractVector{<:PottsEvolver.CodonSequence}
    
    # Test mixed case 2 - Float Teq (convertible to Int) and :discrete sampling_type
    # Should work by converting to Int
    params_mixed2 = SamplingParameters(; Teq=5.0, sampling_type=:discrete)
    S_mixed2, _ = mcmc_sample(g, M, params_mixed2; alignment_output=false, init=:random_num)
    @test S_mixed2 isa AbstractVector{<:PottsEvolver.NumSequence}
end

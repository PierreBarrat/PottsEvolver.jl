using BioSequenceMappings

@testset "SamplingParameters" begin
    @test SamplingParameters(; Teq=5) isa Any
    @test_throws AssertionError SamplingParameters(; Teq=5, step_type=:dubstep)
    @test_throws AssertionError SamplingParameters(; Teq=5, step_meaning=:olive_tree)
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
end

@testset "Output values" begin
    L, q, M = (5, 21, 2)
    g = PottsGraph(L, q; init=:rand)
    @test g.alphabet == aa_alphabet # used q = 21

    params = SamplingParameters(; Teq=1)

    ## Num sequence
    # because asked explicitely
    S, _ = mcmc_sample(g, M, params; init=:random_num, alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.NumSequence}

    # because g has no alphabet
    g_noalphabet = let
        x = deepcopy(g)
        x.alphabet = nothing
        x
    end
    S, _ = mcmc_sample(g_noalphabet, M, params; init=[1, 2, 3], alignment_output=false)
    @test S isa AbstractVector{<:PottsEvolver.NumSequence}

    # because g has an alphabet that is not the default aa
    g_strangealphabet = let
        x = deepcopy(g)
        x.alphabet = Alphabet("ACDEFGHIKLMNPQRSTVWY-")
        x
    end
    S, _ = mcmc_sample(g_strangealphabet, M, params; init=[1, 2, 3], alignment_output=false)
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

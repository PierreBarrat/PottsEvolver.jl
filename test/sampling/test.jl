@testset "SamplingParameters" begin
    @test SamplingParameters() isa Any
    @test_throws ArgumentError SamplingParameters(; step_type=:dubstep)
    @test_throws ArgumentError SamplingParameters(; step_meaning=:olive_tree)

    @testset "the schedule is not a sampling parameter" begin
        # `Teq` and `burnin` moved to `mcmc_sample`, as `spacing` and `burnin`
        p = SamplingParameters()
        @test !hasproperty(p, :Teq)
        @test !hasproperty(p, :burnin)
        @test_throws MethodError SamplingParameters(; Teq=5)
        @test_throws MethodError SamplingParameters(; burnin=5)
    end

    @testset "mutation_matrix" begin
        cont = (; sampling_type=:continuous, step_type=:glauber)

        # default is nothing
        @test SamplingParameters(; cont...).mutation_matrix === nothing

        # symmetric matrix: stored, converted to Float64
        μ = [0 1 1; 1 0 1; 1 1 0]
        p = SamplingParameters(; cont..., mutation_matrix=μ)
        @test p.mutation_matrix isa Matrix{Float64}
        @test p.mutation_matrix == μ

        # a genuine generator matrix (negative diagonal) is fine; diagonal is ignored
        @test SamplingParameters(;
            cont..., mutation_matrix=[-2.0 1 1; 1 -2 1; 1 1 -2]
        ) isa Any

        # invalid: non-square, or negative off-diagonal
        @test_throws ArgumentError SamplingParameters(;
            cont..., mutation_matrix=[1.0 2 3; 4 5 6]
        )
        @test_throws ArgumentError SamplingParameters(;
            cont..., mutation_matrix=[0.0 -1; 1 0]
        )

        # warnings: asymmetric mu (equilibrium shift), and mu on a discrete run (ignored)
        μasym = [0 2 1; 1 0 1; 1 1 0]
        @test_logs (:warn,) match_mode = :any SamplingParameters(;
            cont..., mutation_matrix=μasym
        )
        @test_logs (:warn,) match_mode = :any SamplingParameters(;
            sampling_type=:discrete, mutation_matrix=μ
        )
    end

    @testset "Construction" begin
        # `SamplingParameters` no longer carries a numeric type parameter: it held only
        # `Teq` and `burnin`, and nothing ever dispatched on it
        @test isconcretetype(SamplingParameters)
    end
end

@testset "Sampling schedule" begin
    L, q, M = 4, 21, 3
    g = PottsGraph(L, q; init=:rand)
    s0 = AASequence(L)
    quiet = (; verbose=-1, progress_meter=false)
    discrete = SamplingParameters(; sampling_type=:discrete)
    continuous = SamplingParameters(;
        sampling_type=:continuous, step_type=:glauber, substitution_rate=1.0
    )

    # `spacing` and `burnin` build the schedule
    @test mcmc_sample(g, 3, s0, discrete; spacing=5, burnin=2, quiet...).tvals == [2, 7, 12]
    # `burnin` defaults to 5 * spacing, as the old `SamplingParameters` default did
    @test mcmc_sample(g, 2, s0, discrete; spacing=4, quiet...).tvals == [20, 24]

    # discrete time counts steps: whole numbers only, wherever they come from
    @test_throws ArgumentError mcmc_sample(g, M, s0, discrete; spacing=5.2, quiet...)
    @test_throws ArgumentError mcmc_sample(
        g, M, s0, discrete; spacing=5, burnin=1.5, quiet...
    )
    @test_throws ArgumentError PottsEvolver.discrete_times([0.0, 1.5])
    # integer-valued floats are fine, and are converted
    @test PottsEvolver.discrete_times([0.0, 2.0]) == [0, 2]
    @test PottsEvolver.discrete_times([0.0, 2.0]) isa Vector{Int}
    @test mcmc_sample(g, M, s0, discrete; spacing=5.0, burnin=1.0, quiet...) isa Any

    # continuous time takes any non-negative real
    @test mcmc_sample(g, M, s0, continuous; spacing=0.5, quiet...) isa Any
    @test_throws ArgumentError mcmc_sample(g, M, s0, continuous; spacing=-1, quiet...)
end

@testset "Accepted steps" begin
    L, q, M = (5, 3, 100)
    g = PottsGraph(L, q; init=:rand)

    # the difference between two samples should always be exactly one.
    params = SamplingParameters(; step_meaning=:changed)
    S, tvals, _ = mcmc_sample(g, M, params; spacing=1, burnin=3)
    @test unique(map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)) == [1]
    @test tvals == collect(3 .+ (0:((M - 1))))

    # difference from init should be 1, then no changes
    params = SamplingParameters(; step_meaning=:changed)
    s0 = NumSequence(L, q)
    S, tvals, _ = mcmc_sample(g, M, s0, params; spacing=0, burnin=1, verbose=-1)
    @test tvals == 1 .+ zeros(Int, M)
    @test allequal(S)
    @test hamming(S[1], s0.seq; normalize=false) == 1

    # the difference between two samples should be strictly less than one
    # (!! this is a statistically true test, but could be false!)
    params = SamplingParameters(; step_meaning=:changed)
    params = @set params.step_meaning = :accepted
    S, tvals, _ = mcmc_sample(g, M, params; spacing=1, burnin=0)
    H = map(i -> hamming(S[i - 1], S[i]; normalize=false), 2:M)
    acc_ratio_1 = sum(H) / length(H)
    @test acc_ratio_1 < 1
    @test tvals == collect(0 .+ (0:(M - 1)))

    # The difference between two samples should be 0
    params = SamplingParameters() # no spacing: every sample at the same time
    S, _ = mcmc_sample(g, M, params; verbose=-1) # warns about the repeated times
    @test all(==(0), map(i -> hamming(S[i - 1], S[i]), 2:M))
end

@testset "Output values" begin
    L, q, M = (4, 21, 2)
    g = PottsGraph(L, q; init=:rand)

    params = SamplingParameters()
    sched = (; spacing=1)

    ## Num sequence
    # because asked explicitely
    S, _ = mcmc_sample(g, M, params; sched..., init=:random_num, pack_output=false)
    @test S isa AbstractVector{<:PottsEvolver.NumSequence}

    # because g is not an amino acid model (q != 21)
    g_num = PottsGraph(L, 5; init=:rand)
    S, _ = mcmc_sample(g_num, M, params; sched..., init=[1, 2, 3, 4], pack_output=false)
    @test S isa AbstractVector{<:PottsEvolver.NumSequence}

    ## AA sequence
    # because asked explicitely
    S, _ = mcmc_sample(g, M, params; sched..., init=:random_aa, pack_output=false)
    @test S isa AbstractVector{<:PottsEvolver.AASequence}

    # because g has q=21 and init vector has elements <= 21
    S, _ = mcmc_sample(g, M, params; sched..., init=[1, 2, 3, 21], pack_output=false)
    @test S isa AbstractVector{<:PottsEvolver.AASequence}

    ## Codon sequence
    # because asked explicitely
    S, _ = mcmc_sample(g, M, params; sched..., init=:random_codon, pack_output=false)
    @test S isa AbstractVector{<:PottsEvolver.CodonSequence}

    # because g has q=21 and init vector has elements > 21
    S, _ = mcmc_sample(g, M, params; sched..., init=[1, 2, 3, 22], pack_output=false)
    @test S isa AbstractVector{<:PottsEvolver.CodonSequence}
end

@testset "Continuous/discrete dispatch" begin
    L, q, M = (4, 21, 2)
    g = PottsGraph(L, q; init=:rand)

    # Test discrete case - integer spacing and :discrete sampling_type
    params_discrete = SamplingParameters(; sampling_type=:discrete)
    S_discrete =
        mcmc_sample(g, M, params_discrete; spacing=Int(5), pack_output=false).sequences
    @test S_discrete isa AbstractVector{<:PottsEvolver.NumSequence}

    # Test continuous case - float spacing and :continuous sampling_type
    params_continuous = SamplingParameters(;
        sampling_type=:continuous, step_type=:glauber, substitution_rate=1.0
    )
    S_continuous =
        mcmc_sample(
            g, M, params_continuous; spacing=5.0, pack_output=false, init=:random_aa
        ).sequences
    @test S_continuous isa AbstractVector{<:PottsEvolver.AASequence}
    #= NEED TO CHECK TESTS BELOW  =#
    # Test mixed case 1 - integer spacing but :continuous sampling_type
    # Should convert to continuous
    params_mixed1 = SamplingParameters(;
        sampling_type=:continuous, step_type=:metropolis, substitution_rate=1.0
    )
    S_mixed1 =
        mcmc_sample(
            g, M, params_mixed1; spacing=Int(5), pack_output=false, init=:random_codon
        ).sequences
    @test S_mixed1 isa AbstractVector{<:PottsEvolver.CodonSequence}

    # Test mixed case 2 - float spacing (convertible to Int) and :discrete sampling_type
    # Should work by converting to Int
    params_mixed2 = SamplingParameters(; sampling_type=:discrete)
    S_mixed2, _ = mcmc_sample(
        g, M, params_mixed2; spacing=5.0, pack_output=false, init=:random_num
    )
    @test S_mixed2 isa AbstractVector{<:PottsEvolver.NumSequence}
end

@testset "Sampling with time values" begin
    L, q = 10, 5
    g = PottsGraph(L, q)

    @testset "Errors" begin
        params_discrete = SamplingParameters(; sampling_type=:discrete)
        params_continuous = SamplingParameters(; sampling_type=:continuous)

        # unsorted, and negative
        for tvals in ([0, 2, 1], [-1, 0, 1])
            @test_throws ArgumentError mcmc_sample(g, tvals, params_discrete)
            @test_throws ArgumentError mcmc_sample(g, Float64.(tvals), params_discrete)
            # continuous: the schedule is fine to coerce, so the error is the missing
            # substitution_rate or the schedule itself - either way an ArgumentError
            @test_throws ArgumentError mcmc_sample(g, tvals, params_continuous)
            @test_throws ArgumentError mcmc_sample(g, Float64.(tvals), params_continuous)
        end

        # integer-valued floats are accepted for discrete sampling, and converted
        @test mcmc_sample(g, [0.0, 1.0, 2.0], params_discrete; verbose=-1) isa Any
        # fractional ones are not
        @test_throws ArgumentError mcmc_sample(g, [0.0, 1.5], params_discrete)

        # continuous sampling requires substitution_rate to be set
        @test_throws ArgumentError mcmc_sample(g, [0.0, 1.0, 2.0], params_continuous)
    end

    # test in discrete case
    @testset "Discrete" begin
        sampling_type = :discrete

        time_steps = 0:20
        M = length(time_steps)

        # the difference between two samples should always be exactly one.
        params = SamplingParameters(; sampling_type, step_meaning=:changed)
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
        params = SamplingParameters(;
            sampling_type, step_meaning=:changed, substitution_rate=1.0
        )

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
    @testset "Discrete - Gibbs" begin
        L, q = 10, 5
        g = PottsGraph(L, q; init=:null)
        params = SamplingParameters(; sampling_type=:discrete)
        sched = (; spacing=10)
        @testset "Input given" begin
            s0 = NumSequence(L, q)
            rng = Random.seed!(Xoshiro(123), 42)
            aln_1 = mcmc_sample(g, 2, s0, params; sched..., rng).sequences
            rng = Random.seed!(Xoshiro(123), 42)
            aln_2 = mcmc_sample(g, 2, s0, params; sched..., rng).sequences
            @test aln_1[1] == aln_2[1] && aln_1[2] == aln_2[2]
        end
        @testset "Input not given" begin
            init = :random_num
            rng = Random.seed!(Xoshiro(123), 42)
            aln_1 = mcmc_sample(g, 2, params; sched..., rng, init).sequences
            rng = Random.seed!(Xoshiro(123), 42)
            aln_2 = mcmc_sample(g, 2, params; sched..., rng, init).sequences
            @test aln_1[1] == aln_2[1] && aln_1[2] == aln_2[2]
        end
    end
    @testset "Discrete - Gibbs - Codon" begin
        L, q = 10, 21 # need 21 to sample from CodonSequence
        g = PottsGraph(L, q; init=:null)

        params = SamplingParameters(; sampling_type=:discrete)
        sched = (; spacing=10, burnin=0)
        init = :random_codon
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, 10, params; sched..., rng, init).sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, 10, params; sched..., rng, init).sequences
        for i in 1:10
            @test aln_1[i] == aln_2[i]
        end
    end

    # Add a test with a tree as input, both discrete and continuous, with different type of branch length interpretation
    @testset "Discrete - Gibbs - Codon - Tree exact" begin
        L, q = 10, 21 # need 21 to sample from CodonSequence
        g = PottsGraph(L, q; init=:null)
        tree = TreeTools.Generate.balanced_binary_tree(8, 1.0)

        params = SamplingParameters(; sampling_type=:discrete)
        init = :random_codon
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, tree, params; rng, init).leaf_sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, tree, params; rng, init).leaf_sequences
        for i in 1:8
            @test aln_1[i] == aln_2[i]
        end
    end

    @testset "Discrete - Gibbs - Codon - Tree poisson" begin
        L, q = 10, 21 # need 21 to sample from CodonSequence
        g = PottsGraph(L, q; init=:null)
        tree = TreeTools.Generate.balanced_binary_tree(8, 0.25)
        params = SamplingParameters(;
            sampling_type=:discrete,
            branchlength_meaning=BranchLengthMeaning(:sweep, :poisson),
        )
        init = :random_codon
        rng = Random.seed!(Xoshiro(123), 42)
        aln_1 = mcmc_sample(g, tree, params; rng, init).leaf_sequences
        rng = Random.seed!(Xoshiro(123), 42)
        aln_2 = mcmc_sample(g, tree, params; rng, init).leaf_sequences
        for i in 1:8
            @test aln_1[i] == aln_2[i]
        end
    end

    @testset "parallel_chains basic functionality" begin
        # Setup minimal test data
        L, q, n_chains = 5, 3, 3  # Small values for fast testing
        g = PottsGraph(L, q; init=:rand)
        params = SamplingParameters()
        time_steps = [0, 5, 10]  # Few time points
        inits = [NumSequence(L, q) for _ in 1:n_chains]

        # Execute the function
        result = PottsEvolver.Parallel.parallel_chains(
            g, time_steps, inits, params; verbose=0
        )

        # Verify output structure
        @test result isa Vector{<:Vector{<:NumSequence}}
        @test length(result) == n_chains  # One chain per input

        # Verify each chain structure
        @test all(length(chain) == length(time_steps) for chain in result)
        @test all(length(chain[1].seq) == L for chain in result)
        @test all(eltype(chain) <: NumSequence for chain in result)
    end
end

@testset "warn_repeated_times" begin
    # no mcmc step happens between samples at the same time: they come out identical
    @test_logs (:warn,) match_mode = :any PottsEvolver.warn_repeated_times([5, 5, 5])
    @test_logs (:warn,) match_mode = :any PottsEvolver.warn_repeated_times([0, 2, 2, 4])
    @test_logs min_level = Logging.Warn PottsEvolver.warn_repeated_times([0, 2, 4])
    @test_logs min_level = Logging.Warn PottsEvolver.warn_repeated_times([3])
    @test_logs min_level = Logging.Warn PottsEvolver.warn_repeated_times(Float64[])

    # which is what a `spacing` of 0 produces, the default
    g = PottsGraph(5, 21; init=:rand)
    params = SamplingParameters(; sampling_type=:discrete)
    result = mcmc_sample(
        g, 3, AASequence(5), params; spacing=0, burnin=1, verbose=-1, progress_meter=false
    )
    @test allequal(eachcol(result.sequences.data))
end

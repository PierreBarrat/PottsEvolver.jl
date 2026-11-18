@testset "regroup_samples" begin
    # M samples of L x N  ->  N samples of L x M
    L, M, N = 4, 3, 2
    data = map(1:M) do m
        SequenceSample{NumSequence{Int,4},Int}(
            reshape([10m + 10n + i for i in 1:L, n in 1:N], L, N), ["a", "b"]
        )
    end
    out = PottsEvolver.regroup_samples(data)
    @test out isa Vector{<:SequenceSample}
    @test length(out) == N
    @test all(X -> size(X) == (L, M), out)
    # out[n][m] == data[m][n], and output labels are the index into `data`
    for n in 1:N, m in 1:M
        @test out[n][m] == data[m][n]
    end
    @test out[1].labels == string.(1:M)

    # the label *order* may differ between samples: lookup is by label, not by position
    shuffled = SequenceSample{NumSequence{Int,4},Int}(data[2].data[:, [2, 1]], ["b", "a"])
    @test PottsEvolver.regroup_samples([data[1], shuffled, data[3]]) == out

    # but the label *set* and the sequence length must agree
    @test_throws ArgumentError PottsEvolver.regroup_samples([
        data[1], SequenceSample{NumSequence{Int,4},Int}(data[2].data, ["a", "c"])
    ])
    @test_throws ArgumentError PottsEvolver.regroup_samples([
        data[1], SequenceSample{NumSequence{Int,4},Int}(data[2].data[1:2, :], ["a", "b"])
    ])
    @test_throws ArgumentError PottsEvolver.regroup_samples(SequenceSample[])
end

@testset "parallel_chains" begin
    L, q, n_chains = 5, 4, 3
    g = PottsGraph(L, q; init=:rand)
    params = SamplingParameters(; sampling_type=:discrete)
    tvals = [0, 5, 10]
    inits = [NumSequence(L, q) for _ in 1:n_chains]

    chains = PottsEvolver.parallel_chains(g, tvals, inits, params; seed=1, verbose=-1)
    @test chains isa Vector{<:SequenceSample}
    @test length(chains) == n_chains            # one entry per chain
    @test all(X -> size(X) == (L, length(tvals)), chains)
    @test all(X -> X.labels == string.(tvals), chains)   # labelled by sampling time

    # `seed` alone reproduces the run, and chain `i` does not depend on the others
    @test PottsEvolver.parallel_chains(g, tvals, inits, params; seed=1, verbose=-1) ==
        chains
    @test PottsEvolver.parallel_chains(g, tvals, inits, params; seed=2, verbose=-1) !=
        chains
    @test first(
        PottsEvolver.parallel_chains(g, tvals, inits[1:1], params; seed=1, verbose=-1)
    ) == first(chains)

    # tvals=[0, ...]: the first sample is the initial sequence, so chains differ by their init
    @test all(i -> chains[i][1] == PottsEvolver.sequence(inits[i]), 1:n_chains)

    @testset "M/spacing/burnin schedule" begin
        chains = PottsEvolver.parallel_chains(
            g, 4, inits, params; spacing=3, burnin=2, seed=1, verbose=-1
        )
        @test all(X -> X.labels == ["2", "5", "8", "11"], chains)
    end

    @testset "drawing the initial sequences" begin
        chains = PottsEvolver.parallel_chains(
            g, tvals, params; n_chains, init=:random_num, seed=1, verbose=-1
        )
        @test length(chains) == n_chains
        # independent inits: the chains are not all the same trajectory
        @test !allequal(chains)
    end

    @testset "regrouped by time" begin
        pertime = PottsEvolver.regroup_samples(chains)
        @test length(pertime) == length(tvals)
        @test all(X -> size(X) == (L, n_chains), pertime)
        for (n, t) in enumerate(tvals), i in 1:n_chains
            @test pertime[n][i] == chains[i][n]
        end
    end

    @test_throws ArgumentError PottsEvolver.parallel_chains(
        g, tvals, params; n_chains=0, verbose=-1
    )
end

@testset "equilibrium_sample" begin
    L, q, n_chains = 5, 4, 10
    g = PottsGraph(L, q; init=:rand)
    params = SamplingParameters(; sampling_type=:discrete)

    X = PottsEvolver.equilibrium_sample(g, n_chains, params; burnin=20, seed=1, verbose=-1)
    @test X isa SequenceSample
    @test size(X) == (L, n_chains)        # one sequence per chain
    @test X.labels == string.(1:n_chains) # labelled by chain index
    @test PottsEvolver.equilibrium_sample(
        g, n_chains, params; burnin=20, seed=1, verbose=-1
    ) == X

    # burnin=0 returns the initial sequences untouched
    X0 = PottsEvolver.equilibrium_sample(
        g, n_chains, params; burnin=0, init=:random_num, seed=1, verbose=-1
    )
    @test size(X0) == (L, n_chains)

    @test_throws ArgumentError PottsEvolver.equilibrium_sample(
        g, n_chains, params; burnin=-1, verbose=-1
    )
end

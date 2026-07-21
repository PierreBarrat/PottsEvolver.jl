@testset "IO" begin
    numerical_graph_file_1 = joinpath(dirname(@__FILE__), "test_numerical_1.pottsgraph")
    g1 = read_graph(numerical_graph_file_1)
    @test isa(g1, PottsGraph)
    L, q = size(g1)
    @test L == 2
    @test q == 21
    for i in 1:L, j in (i + 1):L, a in 1:q, b in 1:q
        if g1.J[a, b, i, j] == 0 # random graph, so technically it is a rand based test
            @test false
            break
        end
        g1.h[a, i] == 0 && (@test false; break)
    end

    numerical_graph_file_0 = joinpath(dirname(@__FILE__), "test_numerical_0.pottsgraph")
    g0 = read_graph(numerical_graph_file_0)
    @test g0.J == g1.J
    @test g0.h == g1.h

    symbolic_graph_file_1 = joinpath(dirname(@__FILE__), "test_symbolic_1.pottsgraph")
    g2 = read_graph(symbolic_graph_file_1)
    @test isa(g2, PottsGraph)
    L, q = size(g2)
    @test L == 2
    @test q == 21

    symbolic_graph_file_0 = joinpath(dirname(@__FILE__), "test_symbolic_0.pottsgraph")
    g3 = read_graph(symbolic_graph_file_0)
    @test g3.J == g2.J
    @test g3.h == g2.h

    # I could (should?) add a test where lines of the file are shuffled?
end

@testset "Gauge change" begin
    g = PottsGraph(5, 3; init=:rand)
    PottsEvolver.set_gauge!(g, :zero_sum)
    for i in 1:5
        @test abs(sum(g.h[:, i])) < 1e-5
        for j in 1:5, b in 1:3
            @test abs(sum(g.J[:, b, i, j])) < 1e-5
        end
    end

    @test_throws ErrorException PottsEvolver.set_gauge!(g, :some_random_gauge)
end

@testset "copy" begin
    rng = MersenneTwister(123)
    g = PottsGraph(5, 21; init=:rand, rng)
    gc = copy(g)

    @test gc.J == g.J && gc.h == g.h && gc.β == g.β
    @test gc !== g && gc.J !== g.J && gc.h !== g.h

    # the copy is independent of the original
    gc.J[1, 2, 1, 2] = 999.0
    gc.h[1, 1] = -555.0
    gc.β = 2.5
    @test g.J[1, 2, 1, 2] != 999.0
    @test g.h[1, 1] != -555.0
    @test g.β == 1.0

    seq = rand(1:21, 5)
    @test PottsEvolver.energy(seq, g) ≈ PottsEvolver.energy(seq, copy(g))
end

@testset "Numerical type" begin
    # `T` must be respected for both `init` values: `:rand` used to ignore it, since
    # `_random_graph` was called without it and `randn` defaults to Float64
    for T in (Float64, Float32), init in (:null, :rand)
        g = PottsGraph(3, 3, T; init)
        @test g isa PottsGraph{T}
        @test eltype(g.J) === eltype(g.h) === typeof(g.β) === T
        @test copy(g) isa PottsGraph{T}
    end

    @test_throws ArgumentError PottsGraph(3, 3; init=:not_a_thing)
end

function params_fields_equal(p, q)
    return all(f -> getproperty(p, f) == getproperty(q, f), propertynames(p))
end

@testset "load_params" begin
    @testset "defaults for omitted fields" begin
        mktemp() do path, io
            write(io, "Teq = 5\n")
            close(io)
            p = PottsEvolver.CLI.load_params(path)
            @test params_fields_equal(p, SamplingParameters(; Teq=5))
        end
    end

    @testset "discrete, all fields set" begin
        mktemp() do path, io
            write(
                io,
                """
                sampling_type = "discrete"
                step_type = "gibbs"
                Teq = 4
                burnin = 10
                step_meaning = "changed"
                fraction_gap_step = 0.8

                [branchlength_meaning]
                type = "sweep"
                length = "poisson"
                """,
            )
            close(io)
            p = PottsEvolver.CLI.load_params(path)
            expected = SamplingParameters(;
                sampling_type=:discrete,
                step_type=:gibbs,
                Teq=4,
                burnin=10,
                step_meaning=:changed,
                fraction_gap_step=0.8,
                branchlength_meaning=BranchLengthMeaning(:sweep, :poisson),
            )
            @test params_fields_equal(p, expected)
        end
    end

    @testset "continuous with mutation_matrix: row convention" begin
        mktemp() do path, io
            write(
                io,
                """
                sampling_type = "continuous"
                step_type = "glauber"
                Teq = 1.0
                substitution_rate = 12.3
                mutation_matrix = [[0, 2, 3], [1, 0, 1], [1, 5, 0]]
                """,
            )
            close(io)
            p = PottsEvolver.CLI.load_params(path)
            # row i of the TOML array must land as row i of the matrix, unpermuted
            @test p.mutation_matrix == [0 2 3; 1 0 1; 1 5 0]
            @test p.mutation_matrix isa Matrix{Float64}
            @test p.substitution_rate == 12.3
        end
    end

    @testset "invalid combination surfaces SamplingParameters's own error" begin
        mktemp() do path, io
            write(
                io,
                """
                sampling_type = "discrete"
                step_type = "glauber"
                """,
            )
            close(io)
            @test_throws ArgumentError PottsEvolver.CLI.load_params(path)
        end
    end
end

include("load_init_test.jl")
include("parse_cli_test.jl")
include("write_output_test.jl")

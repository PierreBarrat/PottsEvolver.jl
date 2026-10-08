@testset "load_tvals" begin
    mktemp() do path, io
        write(
            io,
            """
            # a comment
            0.0
              1.5

            3.0
            """,
        )
        close(io)
        @test CLI.load_tvals(path) == [0.0, 1.5, 3.0]
    end

    mktemp() do path, io
        write(io, "0.0\nnot_a_number\n")
        close(io)
        @test_throws ArgumentError CLI.load_tvals(path)
    end
end

@testset "sampling_schedule" begin
    discrete = SamplingParameters(; Teq=2)
    continuous = SamplingParameters(;
        sampling_type=:continuous, step_type=:glauber, substitution_rate=1.0
    )

    @testset "one of the three, or none" begin
        @test CLI.sampling_schedule(nothing, Float64[], nothing, discrete) == 1 # default
        @test CLI.sampling_schedule(7, Float64[], nothing, discrete) == 7
        @test CLI.sampling_schedule(nothing, [0.0, 1.5], nothing, continuous) == [0.0, 1.5]
        mktemp() do path, io
            write(io, "0.0\n2.0\n")
            close(io)
            @test CLI.sampling_schedule(nothing, Float64[], path, continuous) == [0.0, 2.0]
        end
    end

    @testset "element type follows the sampling type" begin
        # discrete sampling counts steps, and its methods dispatch on `Integer` time steps
        @test CLI.sampling_schedule(nothing, [0.0, 2.0], nothing, discrete) isa Vector{Int}
        @test CLI.sampling_schedule(nothing, [0.0, 2.0], nothing, continuous) isa
            Vector{Float64}
        @test_throws ArgumentError CLI.sampling_schedule(
            nothing, [0.0, 1.5], nothing, discrete
        )
    end

    @testset "combinations are refused" begin
        @test_throws ArgumentError CLI.sampling_schedule(3, [0.0, 1.0], nothing, discrete)
        @test_throws ArgumentError CLI.sampling_schedule(
            3, Float64[], "tvals.txt", discrete
        )
        @test_throws ArgumentError CLI.sampling_schedule(
            nothing, [0.0], "tvals.txt", discrete
        )
        # an explicit `-M 1` still conflicts: this is why the default is `nothing`
        @test_throws ArgumentError CLI.sampling_schedule(1, [0.0, 1.0], nothing, discrete)
    end
end

@testset "resolve_seed" begin
    seed, rng = CLI.resolve_seed(42)
    @test seed == 42
    @test rand(rng, 1:1000) == rand(Random.Xoshiro(42), 1:1000)

    # a seed is drawn when none is given, so that [run] can always record one
    drawn, rng = CLI.resolve_seed(nothing)
    @test drawn isa Integer
    @test rand(rng, 1:1000) == rand(Random.Xoshiro(drawn), 1:1000)
end

@testset "compute_omega!" begin
    g = PottsGraph(5, 21; init=:rand)
    s0 = AASequence(5)

    @testset "sets the rate when the file does not" begin
        params = SamplingParameters(; sampling_type=:continuous, step_type=:glauber)
        @test isnothing(params.substitution_rate)
        CLI.compute_omega!(params, g, s0; progress_meter=false)
        @test params.substitution_rate > 0
    end

    @testset "the parameter file wins" begin
        params = SamplingParameters(;
            sampling_type=:continuous, step_type=:glauber, substitution_rate=3.5
        )
        @test_logs (:warn,) match_mode = :any CLI.compute_omega!(
            params, g, s0; progress_meter=false
        )
        @test params.substitution_rate == 3.5
    end
end

@testset "cli_settings: usage string names the program" begin
    # without `prog`, ArgParse writes `<PROGRAM>` in the usage line and in --help
    usage = ArgParse.usage_string(CLI.cli_settings())
    @test occursin("pottsevolver", usage)
    @test occursin("sample-chain", usage)
end

@testset "parse errors report the message, not the exception" begin
    # `showerror` on an ArgParseError prints its constructor; `_error_message` unwraps it
    err = ArgParse.ArgParseError("no command given")
    @test CLI._error_message(err) == "no command given"
    # other errors keep `showerror`'s output, where the type is informative
    @test CLI._error_message(ArgumentError("boom")) == "ArgumentError: boom"

    # and they are user errors: exit code 1, no exception escaping
    @test CLI.main(String[]) == 1
    @test CLI.main(["sample-chain", "--no-such-flag"]) == 1
end

@testset "main: end to end" begin
    # a model and a parameter file on disk, as a user would have them
    function setup(dir; params_content="Teq = 2\nburnin = 1\n")
        model = joinpath(dir, "model.dat")
        PottsEvolver.write(model, PottsGraph(4, 21; init=:rand))
        paramfile = joinpath(dir, "params.toml")
        write(paramfile, params_content)
        return model, paramfile
    end

    @testset "sample-chain" begin
        mktempdir() do dir
            model, paramfile = setup(dir)
            out = joinpath(dir, "out")
            code = CLI.main([
                "sample-chain",
                model,
                "-p",
                paramfile,
                "-o",
                out,
                "-M",
                "4",
                "--seed",
                "1234",
                "--verbose=-1",
                "--log-verbose",
                "1",
            ])
            @test code == 0
            @test Set(readdir(out)) ==
                Set(["model.fasta", "model.log", "model.params.toml"])

            # the [run] table records what is needed to reproduce the run
            written = TOML.parsefile(joinpath(out, "model.params.toml"))
            @test written["run"]["seed"] == 1234
            @test written["run"]["model"] == model
            @test written["run"]["init"] == "random_aa"
            @test written["run"]["sequence_type"] == "AASequence{Int64}"
            # and the parameters round-trip back into the same run
            @test params_fields_equal(
                CLI.load_params(joinpath(out, "model.params.toml")),
                CLI.load_params(paramfile),
            )

            # 4 sequences, and the log covers the whole run
            @test length(fasta_labels(joinpath(out, "model.fasta"))) == 4
            log = read(joinpath(out, "model.log"), String)
            @test occursin("Reading Potts model", log)  # CLI preamble
            @test occursin("Sampling done", log)        # mcmc_sample
            @test occursin("Wrote 4 sequences", log)    # CLI output writing
        end
    end

    @testset "sample-chain with --prefix and --tvals" begin
        mktempdir() do dir
            model, paramfile = setup(dir)
            out = joinpath(dir, "out")
            code = CLI.main([
                "sample-chain",
                model,
                "-p",
                paramfile,
                "-o",
                out,
                "--prefix",
                "mychain",
                "--tvals",
                "0",
                "2",
                "4",
                "--verbose=-1",
            ])
            @test code == 0
            @test Set(readdir(out)) ==
                Set(["mychain.fasta", "mychain.log", "mychain.params.toml"])
            # discrete sampling: the schedule became Int, so the labels are step counts
            @test fasta_labels(joinpath(out, "mychain.fasta")) == ["0", "2", "4"]
        end
    end

    @testset "sample-tree --internals" begin
        mktempdir() do dir
            model, paramfile = setup(dir)
            treefile = joinpath(dir, "tree.nwk")
            TreeTools.write_newick(treefile, balanced_binary_tree(4, 1.0))
            out = joinpath(dir, "out")
            code = CLI.main([
                "sample-tree",
                model,
                treefile,
                "-p",
                paramfile,
                "-o",
                out,
                "--internals",
                "--verbose=-1",
            ])
            @test code == 0
            @test Set(readdir(out)) == Set([
                "model.fasta",
                "model.internals.fasta",
                "model.log",
                "model.nwk",
                "model.params.toml",
            ])
            @test TOML.parsefile(joinpath(out, "model.params.toml"))["run"]["tree"] ==
                treefile
            @test length(fasta_labels(joinpath(out, "model.fasta"))) == 4
        end
    end

    @testset "a seed is drawn and recorded when --seed is omitted" begin
        mktempdir() do dir
            model, paramfile = setup(dir)
            out = joinpath(dir, "out")
            @test CLI.main([
                "sample-chain", model, "-p", paramfile, "-o", out, "--verbose=-1"
            ]) == 0
            @test TOML.parsefile(joinpath(out, "model.params.toml"))["run"]["seed"] isa
                Integer
        end
    end

    @testset "user errors give exit code 1, not a stacktrace" begin
        mktempdir() do dir
            model, paramfile = setup(dir)
            out = joinpath(dir, "out")

            # conflicting schedule flags
            @test CLI.main([
                "sample-chain",
                model,
                "-p",
                paramfile,
                "-o",
                out,
                "-M",
                "3",
                "--tvals",
                "0",
                "1",
                "--verbose=-1",
            ]) == 1

            # missing model file
            @test CLI.main([
                "sample-chain",
                joinpath(dir, "nope.dat"),
                "-p",
                paramfile,
                "-o",
                out,
                "--verbose=-1",
            ]) == 1

            # continuous sampling without substitution_rate: the library's own ArgumentError
            continuous = joinpath(dir, "continuous.toml")
            write(continuous, "sampling_type = \"continuous\"\nTeq = 0.5\n")
            @test CLI.main([
                "sample-chain", model, "-p", continuous, "-o", out, "--verbose=-1"
            ]) == 1
        end
    end

    @testset "--compute-omega fills in the missing rate" begin
        mktempdir() do dir
            model, _ = setup(dir)
            paramfile = joinpath(dir, "continuous.toml")
            write(
                paramfile,
                """
                sampling_type = "continuous"
                step_type = "glauber"
                Teq = 0.5
                """,
            )
            out = joinpath(dir, "out")
            @test CLI.main([
                "sample-chain",
                model,
                "-p",
                paramfile,
                "-o",
                out,
                "-M",
                "2",
                "--compute-omega",
                "--verbose=-1",
            ]) == 0
            written = CLI.load_params(joinpath(out, "model.params.toml"))
            @test written.substitution_rate > 0
        end
    end
end

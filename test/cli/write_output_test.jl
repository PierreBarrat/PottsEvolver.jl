using TOML
import TreeTools.Generate: balanced_binary_tree

const CLI = PottsEvolver.CLI

function fasta_labels(file)
    return map(l -> chop(l; head=1, tail=0), filter(startswith('>'), readlines(file)))
end

@testset "default_prefix" begin
    @test CLI.default_prefix(nothing, joinpath("some", "dir", "graph.dat")) == "graph"
    @test CLI.default_prefix(nothing, "graph") == "graph"
    @test CLI.default_prefix("given", joinpath("some", "dir", "graph.dat")) == "given"
end

@testset "write_params" begin
    @testset "round-trip through load_params" begin
        params = SamplingParameters(;
            sampling_type=:continuous,
            step_type=:glauber,
            Teq=1.5,
            burnin=3.0,
            substitution_rate=2.5,
            track_substitutions=true,
            mutation_matrix=[0 2 3; 1 0 1; 1 5 0],
            branchlength_meaning=BranchLengthMeaning(:sweep, :poisson),
        )
        mktempdir() do dir
            path = CLI.write_params(
                joinpath(dir, "run.params.toml"),
                params;
                run=Dict(:model => "graph.dat", :seed => 42),
            )
            @test params_fields_equal(CLI.load_params(path), params)

            raw = TOML.parsefile(path)
            @test raw["run"]["model"] == "graph.dat"
            @test raw["run"]["seed"] == 42
            # the matrix is written as rows, the layout `load_params` reads
            @test raw["mutation_matrix"] ==
                [[0.0, 2.0, 3.0], [1.0, 0.0, 1.0], [1.0, 5.0, 0.0]]
            @test raw["branchlength_meaning"] ==
                Dict("type" => "sweep", "length" => "poisson")
        end
    end

    @testset "fields that are nothing are omitted" begin
        params = SamplingParameters(; Teq=2) # substitution_rate, mutation_matrix unset
        mktempdir() do dir
            path = CLI.write_params(joinpath(dir, "run.params.toml"), params)
            raw = TOML.parsefile(path)
            @test !haskey(raw, "substitution_rate")
            @test !haskey(raw, "mutation_matrix")
            @test !haskey(raw, "run") # nothing to put there
            @test params_fields_equal(CLI.load_params(path), params)
        end
    end

    @testset "sequence_type goes to [run], which load_params ignores" begin
        params = SamplingParameters(; Teq=2)
        s0 = AASequence(5)
        mktempdir() do dir
            path = CLI.write_params(
                joinpath(dir, "run.params.toml"), PottsEvolver.return_params(params, s0)
            )
            @test TOML.parsefile(path)["run"]["sequence_type"] == string(typeof(s0))
            @test params_fields_equal(CLI.load_params(path), params)
        end
    end
end

@testset "write_substitutions" begin
    info = [
        (;
            number_substitutions=2,
            substitutions=[
                (PottsEvolver.Mutation(3, 1, 5), 0.5), (PottsEvolver.Mutation(1, 2, 4), 0.9)
            ],
        ),
        (; number_substitutions=0, substitutions=Tuple{PottsEvolver.Mutation,Float64}[]),
    ]
    mktempdir() do dir
        path = CLI.write_substitutions(joinpath(dir, "run.csv"), info, [0.0, 1.0, 2.0])
        # `info[m]` is the interval ending at `tvals[m+1]`
        @test readlines(path) ==
            ["sample,time,position,old,new", "1.0,0.5,3,1,5", "1.0,0.9,1,2,4"]
    end
end

@testset "write_chain_output" begin
    L, q = 5, 21
    g = PottsGraph(L, q; init=:rand)

    @testset "discrete" begin
        params = SamplingParameters(; Teq=2, burnin=1)
        result = mcmc_sample(
            g, 4, params; init=:random_aa, verbose=-1, progress_meter=false
        )
        mktempdir() do dir
            files = CLI.write_chain_output(dir, "run", result)
            @test Set(basename.(files)) == Set(["run.fasta", "run.params.toml"])
            @test Set(readdir(dir)) == Set(["run.fasta", "run.params.toml"])
            # sequences are labelled by sampling time
            @test fasta_labels(joinpath(dir, "run.fasta")) == string.(result.tvals)
            # written parameters describe the run that produced them
            @test params_fields_equal(
                CLI.load_params(joinpath(dir, "run.params.toml")), params
            )
        end
    end

    @testset "continuous with tracked substitutions" begin
        params = SamplingParameters(;
            sampling_type=:continuous,
            step_type=:glauber,
            Teq=0.5,
            burnin=0.1,
            substitution_rate=1.0,
            track_substitutions=true,
        )
        result = mcmc_sample(
            g, 3, params; init=:random_aa, verbose=-1, progress_meter=false
        )
        mktempdir() do dir
            files = CLI.write_chain_output(dir, "run", result)
            @test Set(readdir(dir)) ==
                Set(["run.fasta", "run.params.toml", "run.substitutions.csv"])
            @test first(readlines(joinpath(dir, "run.substitutions.csv"))) ==
                "sample,time,position,old,new"
        end
    end

    @testset "existing files are overwritten with a warning" begin
        params = SamplingParameters(; Teq=2)
        result = mcmc_sample(
            g, 3, params; init=:random_aa, verbose=-1, progress_meter=false
        )
        mktempdir() do dir
            CLI.write_chain_output(dir, "run", result)
            @test_logs (:warn,) match_mode = :any CLI.write_chain_output(dir, "run", result)
            @test Set(readdir(dir)) == Set(["run.fasta", "run.params.toml"])
        end
    end
end

@testset "write_tree_output" begin
    g = PottsGraph(5, 21; init=:rand)
    tree = balanced_binary_tree(8, 1.0)
    params = SamplingParameters(; Teq=1, burnin=1)
    result = mcmc_sample(g, tree, params; init=:random_aa, verbose=-1)

    @testset "leaves only" begin
        mktempdir() do dir
            files = CLI.write_tree_output(dir, "run", result)
            @test Set(basename.(files)) == Set(["run.fasta", "run.params.toml"])
            @test Set(readdir(dir)) == Set(["run.fasta", "run.params.toml"])
        end
    end

    @testset "--internals writes the sampled tree too" begin
        mktempdir() do dir
            CLI.write_tree_output(dir, "run", result; internals=true)
            @test Set(readdir(dir)) ==
                Set(["run.fasta", "run.params.toml", "run.internals.fasta", "run.nwk"])

            # the written tree is the one the sequence labels refer to
            written = read_tree(joinpath(dir, "run.nwk"))
            @test Set(fasta_labels(joinpath(dir, "run.fasta"))) ==
                Set(map(label, leaves(written)))
            @test Set(fasta_labels(joinpath(dir, "run.internals.fasta"))) ==
                Set(map(label, internals(written)))
        end
    end

    @testset "track_substitutions warns: no substitutions on a tree" begin
        tracked = SamplingParameters(;
            sampling_type=:continuous,
            step_type=:glauber,
            substitution_rate=1.0,
            track_substitutions=true,
        )
        tracked_result = mcmc_sample(g, tree, tracked; init=:random_aa, verbose=-1)
        mktempdir() do dir
            @test_logs (:warn,) match_mode = :any CLI.write_tree_output(
                dir, "run", tracked_result
            )
            @test !isfile(joinpath(dir, "run.substitutions.csv"))
        end
    end
end

@testset "log file shared by the CLI and mcmc_sample" begin
    g = PottsGraph(5, 21; init=:rand)
    params = SamplingParameters(; Teq=2, burnin=1)
    mktempdir() do dir
        logpath, io = CLI.setup_logfile(dir, "run", "sample-chain", ["graph.dat", "-p"])
        @test logpath == joinpath(dir, "run.log")

        # one sink, hence one handle and one file position, shared by both loggers
        sink = PottsEvolver.file_sink(io, 1)
        with_logger(PottsEvolver.get_logger(-1, nothing, 1; extra_sinks=(sink,))) do
            @info "CLI PREAMBLE"
            mcmc_sample(
                g,
                3,
                params;
                init=:random_aa,
                verbose=-1,
                extra_sinks=(sink,),
                progress_meter=false,
            )
            @info "CLI TAIL"
        end
        close(io)

        content = read(logpath, String)
        @test startswith(content, "# pottsevolver ")
        positions = map(
            s -> first(something(findfirst(s, content), 0:0)),
            ["CLI PREAMBLE", "Sampling done", "CLI TAIL"],
        )
        @test all(>(0), positions) # nothing was lost
        @test issorted(positions)  # nothing was overwritten out of order
    end
end

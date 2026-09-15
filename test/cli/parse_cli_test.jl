@testset "parse_cli" begin
    @testset "required --params enforced" begin
        @test_throws ArgParse.ArgParseError PottsEvolver.CLI.parse_cli(
            ["sample-chain", "model.dat"]
        )
        @test_throws ArgParse.ArgParseError PottsEvolver.CLI.parse_cli(
            ["sample-tree", "model.dat", "tree.nwk"]
        )
    end

    @testset "required positionals enforced" begin
        @test_throws ArgParse.ArgParseError PottsEvolver.CLI.parse_cli(
            ["sample-chain", "-p", "params.toml"]
        )
        @test_throws ArgParse.ArgParseError PottsEvolver.CLI.parse_cli(
            ["sample-tree", "model.dat", "-p", "params.toml"]
        )
    end

    @testset "sample-chain: defaults" begin
        r = PottsEvolver.CLI.parse_cli(["sample-chain", "model.dat", "-p", "params.toml"])
        @test r["%COMMAND%"] == "sample-chain"
        sub = r["sample-chain"]
        @test sub["model"] == "model.dat"
        @test sub["params"] == "params.toml"
        @test sub["outdir"] == pwd()
        @test sub["prefix"] === nothing
        @test sub["init"] == "random_aa"
        @test sub["init-type"] === nothing
        @test sub["seed"] === nothing
        @test sub["verbose"] == 0
        @test sub["log-verbose"] == 1
        @test sub["compute-omega"] == false
        @test sub["translate"] == false
        @test sub["n-samples"] == 1
        @test sub["tvals"] == Float64[]
        @test sub["tvals-file"] === nothing
    end

    @testset "sample-tree: defaults" begin
        r = PottsEvolver.CLI.parse_cli(
            ["sample-tree", "model.dat", "tree.nwk", "-p", "params.toml"]
        )
        @test r["%COMMAND%"] == "sample-tree"
        sub = r["sample-tree"]
        @test sub["model"] == "model.dat"
        @test sub["tree"] == "tree.nwk"
        @test sub["params"] == "params.toml"
        @test sub["internals"] == false
    end

    @testset "sample-chain: --tvals collects a vector, overrides given" begin
        r = PottsEvolver.CLI.parse_cli([
            "sample-chain",
            "model.dat",
            "-p",
            "params.toml",
            "--tvals",
            "0.0",
            "1.5",
            "3.0",
            "-M",
            "7",
            "--seed",
            "42",
            "-v",
            "2",
            "--compute-omega",
            "--translate",
        ])
        sub = r["sample-chain"]
        @test sub["tvals"] == [0.0, 1.5, 3.0]
        @test sub["n-samples"] == 7
        @test sub["seed"] == 42
        @test sub["verbose"] == 2
        @test sub["compute-omega"] == true
        @test sub["translate"] == true
    end

    @testset "sample-tree: --internals, --init with fasta,seqname syntax" begin
        r = PottsEvolver.CLI.parse_cli([
            "sample-tree",
            "model.dat",
            "tree.nwk",
            "-p",
            "params.toml",
            "--internals",
            "--init",
            "root.fasta,seq3",
            "--init-type",
            "aa",
        ])
        sub = r["sample-tree"]
        @test sub["internals"] == true
        @test sub["init"] == "root.fasta,seq3"
        @test sub["init-type"] == "aa"
    end
end

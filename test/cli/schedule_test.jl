@testset "schedule_times" begin
    @testset "linear" begin
        spec = (; kind=:linear, burnin=100, spacing=50, n_samples=4)
        @test CLI.schedule_times(spec) == [100, 150, 200, 250]
        # a single sample is just the burnin
        @test CLI.schedule_times((; kind=:linear, burnin=7, spacing=50, n_samples=1)) == [7]
        # spacing 0 collapses, which the library warns about rather than refusing
        @test CLI.schedule_times((; kind=:linear, burnin=3, spacing=0, n_samples=3)) ==
            [3, 3, 3]
    end

    @testset "linear, as start/stop" begin
        # the same schedule as burnin=0, spacing=25, n_samples=5, said the other way
        spec = (; kind=:linear, start=0, stop=100, n_samples=5)
        @test CLI.schedule_times(spec) == [0, 25, 50, 75, 100]
        @test CLI.schedule_times(spec) ==
            CLI.schedule_times((; kind=:linear, burnin=0, spacing=25, n_samples=5))
        # `range` rejects a single point unless start == stop, so it is taken directly
        @test CLI.schedule_times((; kind=:linear, start=7, stop=99, n_samples=1)) == [7]
    end

    @testset "log" begin
        spec = (; kind=:log, start=1, stop=1000, n_samples=4)
        @test CLI.schedule_times(spec) ≈ [1, 10, 100, 1000]
        # geometric spacing is undefined from 0
        @test_throws ArgumentError CLI.schedule_times(
            (; kind=:log, start=0, stop=10, n_samples=3)
        )
        # a single point needs no ratio, so `start = 0` is fine there
        @test CLI.schedule_times((; kind=:log, start=0, stop=10, n_samples=1)) == [0]
    end

    @testset "explicit" begin
        @test CLI.schedule_times((; kind=:explicit, tvals=[0, 2, 7])) == [0, 2, 7]
    end

    @test_throws ArgumentError CLI.schedule_times((; kind=:quadratic, n_samples=3))
end

@testset "parse_tvals" begin
    # a linear table is read leniently: the flags may complete it, so unset stays `nothing`
    @test CLI.parse_tvals(Dict("spacing" => 50, "n_samples" => 3)) ==
        (; kind=:linear, burnin=nothing, spacing=50, n_samples=3)
    @test CLI.normalize_schedule(CLI.parse_tvals(Dict("spacing" => 50, "n_samples" => 3))) ==
        (; kind=:linear, burnin=250, spacing=50, n_samples=3)
    # and an empty one normalises to a single sample
    @test CLI.normalize_schedule(CLI.parse_tvals(Dict{String,Any}())) ==
        (; kind=:linear, burnin=0, spacing=0, n_samples=1)
    @test CLI.parse_tvals(
        Dict("kind" => "linear", "burnin" => 10, "spacing" => 5, "n_samples" => 2)
    ) == (; kind=:linear, burnin=10, spacing=5, n_samples=2)
    @test CLI.parse_tvals(
        Dict("kind" => "log", "start" => 1, "stop" => 100, "n_samples" => 3)
    ) == (; kind=:log, start=1, stop=100, n_samples=3)
    # `stop` selects the start/stop shape, with `linear` still the default kind
    @test CLI.parse_tvals(Dict("start" => 0, "stop" => 100, "n_samples" => 5)) ==
        (; kind=:linear, start=0, stop=100, n_samples=5)
    @test CLI.parse_tvals(
        Dict("kind" => "linear", "start" => 0, "stop" => 100, "n_samples" => 5)
    ) == (; kind=:linear, start=0, stop=100, n_samples=5)
    # it is a whole shape: no flag completes it, so a missing entry is an error
    @test_throws ArgumentError CLI.parse_tvals(Dict("start" => 0, "stop" => 100))
    @test_throws ArgumentError CLI.parse_tvals(Dict("stop" => 100, "n_samples" => 5))
    @test CLI.parse_tvals(Dict("tvals" => [1, 2, 4])).tvals == [1, 2, 4]

    # a log table cannot be completed from the command line, so it must be whole, and a
    # missing entry names itself rather than surfacing as a KeyError
    @test_throws ArgumentError CLI.parse_tvals(
        Dict("kind" => "log", "start" => 1, "n_samples" => 3)
    ) # no stop
    @test_throws ArgumentError CLI.parse_tvals(Dict("kind" => "quadratic"))

    @testset "tvals_file, resolved relative to the parameter file" begin
        mktempdir() do dir
            write(joinpath(dir, "sched.txt"), "0\n5\n10\n")
            spec = CLI.parse_tvals(Dict("tvals_file" => "sched.txt"); dir)
            @test spec.tvals == [0.0, 5.0, 10.0]
        end
    end
end

@testset "write_tvals" begin
    @testset "a structured schedule stays a handful of numbers" begin
        # the point of storing a spec: this is a 100_000-sample schedule
        spec = (; kind=:log, start=1.0, stop=1e5, n_samples=100_000)
        table, sidecar = CLI.write_tvals(spec)
        @test isnothing(sidecar)
        @test table == Dict(
            "kind" => "log", "start" => 1.0, "stop" => 1e5, "n_samples" => 100_000
        )
        # and it round-trips
        @test CLI.parse_tvals(table) == spec

        spec = (; kind=:linear, burnin=100, spacing=50, n_samples=4)
        table, _ = CLI.write_tvals(spec)
        @test CLI.parse_tvals(table) == spec

        # a start/stop schedule is written back as start/stop, not restated as a spacing,
        # which also avoids a float wobble when the spacing does not divide evenly
        spec = (; kind=:linear, start=0, stop=100, n_samples=4)
        table, _ = CLI.write_tvals(spec)
        @test table == Dict(
            "kind" => "linear", "start" => 0, "stop" => 100, "n_samples" => 4
        )
        @test CLI.parse_tvals(table) == spec
    end

    @testset "an explicit schedule keeps the form it was given in" begin
        # typed out on the command line: short by construction, so inline it
        spec = (; kind=:explicit, tvals=[0, 5, 10])
        table, sidecar = CLI.write_tvals(spec; prefix_path=joinpath("out", "run"))
        @test isnothing(sidecar)
        @test table == Dict("tvals" => [0, 5, 10])

        # read from a file: stored as a file again, so the parameter file stays readable
        spec = (; kind=:explicit, tvals=[0.0, 5.0], from_file=true)
        table, sidecar = CLI.write_tvals(spec; prefix_path=joinpath("out", "run"))
        @test table == Dict("tvals_file" => "run.tvals.txt")
        @test sidecar == (("out", "run.tvals.txt") => [0.0, 5.0])
    end
end

@testset "resolve_schedule" begin
    # the parsed-args shape `resolve_schedule` reads
    function cli_args(; n_samples=nothing, spacing=nothing, burnin=nothing,
                      tvals=Float64[], tvals_file=nothing)
        return Dict{String,Any}(
            "n-samples" => n_samples,
            "spacing" => spacing,
            "burnin" => burnin,
            "tvals" => tvals,
            "tvals-file" => tvals_file,
        )
    end

    function with_params(content)
        path, io = mktemp()
        write(io, content)
        close(io)
        return path
    end

    @testset "flags alone" begin
        empty_file = with_params("")
        spec = CLI.resolve_schedule(
            cli_args(; n_samples=100, spacing=50, burnin=10), empty_file
        )
        @test spec == (; kind=:linear, burnin=10.0, spacing=50.0, n_samples=100)

        # burnin defaults to 5 * spacing, as the library does
        spec = CLI.resolve_schedule(cli_args(; n_samples=2, spacing=4), empty_file)
        @test spec.burnin == 20.0

        # nothing at all: a single sample
        @test CLI.resolve_schedule(cli_args(), empty_file) ==
            (; kind=:linear, burnin=0, spacing=0, n_samples=1)
        # the spec's field order is canonical, so equal schedules compare equal
        @test keys(CLI.resolve_schedule(cli_args(; n_samples=2, spacing=1), empty_file)) ==
            (:kind, :burnin, :spacing, :n_samples)

        # explicit times
        @test CLI.resolve_schedule(cli_args(; tvals=[0.0, 2.0]), empty_file).tvals ==
            [0.0, 2.0]
    end

    @testset "the file supplies what the flags leave out" begin
        file = with_params("[schedule]\nspacing = 50\nburnin = 100\nn_samples = 7\n")

        # no flags: the file is the schedule
        @test CLI.resolve_schedule(cli_args(), file) ==
            (; kind=:linear, burnin=100, spacing=50, n_samples=7)

        # a flag overrides its own entry and leaves the others alone
        spec = CLI.resolve_schedule(cli_args(; n_samples=3), file)
        @test spec == (; kind=:linear, burnin=100, spacing=50, n_samples=3)
        spec = CLI.resolve_schedule(cli_args(; spacing=5), file)
        @test spec == (; kind=:linear, burnin=100, spacing=5.0, n_samples=7)
    end

    @testset "a start/stop table is restated so a flag can override one entry" begin
        file = with_params("[schedule]\nstart = 0\nstop = 100\nn_samples = 5\n")
        # untouched without flags: the shape is preserved
        @test CLI.resolve_schedule(cli_args(), file) ==
            (; kind=:linear, start=0, stop=100, n_samples=5)
        # with a flag: same times, restated as burnin/spacing, then overridden
        @test CLI.resolve_schedule(cli_args(; n_samples=3), file) ==
            (; kind=:linear, burnin=0, spacing=25.0, n_samples=3)
        @test CLI.resolve_schedule(cli_args(; spacing=10), file) ==
            (; kind=:linear, burnin=0, spacing=10.0, n_samples=5)
    end

    @testset "a log schedule comes from the file" begin
        file = with_params("[schedule]\nkind = \"log\"\nstart = 1\nstop = 100\nn_samples = 3\n")
        @test CLI.resolve_schedule(cli_args(), file) ==
            (; kind=:log, start=1, stop=100, n_samples=3)
        # linear flags cannot complete a log table: it is ignored, with a warning
        spec = @test_logs (:warn,) match_mode = :any CLI.resolve_schedule(
            cli_args(; n_samples=4, spacing=2), file
        )
        @test spec == (; kind=:linear, burnin=10.0, spacing=2.0, n_samples=4)
    end

    @testset "conflicts" begin
        empty_file = with_params("")
        # giving the times and describing them at once
        @test_throws ArgumentError CLI.resolve_schedule(
            cli_args(; n_samples=3, tvals=[0.0, 1.0]), empty_file
        )
        @test_throws ArgumentError CLI.resolve_schedule(
            cli_args(; spacing=5, tvals=[0.0, 1.0]), empty_file
        )
        # two ways of giving the times
        mktempdir() do dir
            f = joinpath(dir, "s.txt")
            write(f, "0\n1\n")
            @test_throws ArgumentError CLI.resolve_schedule(
                cli_args(; tvals=[0.0, 1.0], tvals_file=f), empty_file
            )
        end
        # explicit times on the command line win over a table, with a warning
        file = with_params("[schedule]\nspacing = 50\nn_samples = 7\n")
        spec = @test_logs (:warn,) match_mode = :any CLI.resolve_schedule(
            cli_args(; tvals=[0.0, 1.0]), file
        )
        @test spec.kind == :explicit
    end
end

@testset "load_init" begin
    @testset "keyword" begin
        g = PottsGraph(6, 21; init=:rand)
        s = PottsEvolver.CLI.load_init("random_aa", nothing, g)
        @test s isa AASequence
        @test length(s) == 6

        g_codon = PottsGraph(6, 21; init=:rand)
        s = PottsEvolver.CLI.load_init("random_codon", nothing, g_codon)
        @test s isa CodonSequence

        g_rna = PottsGraph(6, 5; init=:rand)
        s = PottsEvolver.CLI.load_init("random_rna", nothing, g_rna)
        @test s isa RNASequence

        g_num = PottsGraph(6, 7; init=:rand)
        s = PottsEvolver.CLI.load_init("random_num", nothing, g_num)
        @test s isa PottsEvolver.NumSequence
    end

    @testset "fasta: single-record file, no seqname" begin
        mktempdir() do dir
            file = joinpath(dir, "root.fasta")
            seq = AASequence("ACDEFG")
            write_fasta(file, [seq]; labels=["root"])

            g = PottsGraph(6, 21; init=:rand)
            s = PottsEvolver.CLI.load_init(file, "aa", g)
            @test s isa AASequence
            @test s.seq == seq.seq

            # multi-record file without a name should fail
            file2 = joinpath(dir, "many.fasta")
            write_fasta(file2, [seq, seq]; labels=["a", "b"])
            @test_throws ArgumentError PottsEvolver.CLI.load_init(file2, "aa", g)
        end
    end

    @testset "fasta: file,seqname (IQ-TREE style)" begin
        mktempdir() do dir
            file = joinpath(dir, "many.fasta")
            seqs = [AASequence("ACDEFG"), AASequence("GFEDCA")]
            write_fasta(file, seqs; labels=["seq1", "seq2"])

            g = PottsGraph(6, 21; init=:rand)
            s = PottsEvolver.CLI.load_init("$file,seq2", "aa", g)
            @test s.seq == seqs[2].seq

            @test_throws ArgumentError PottsEvolver.CLI.load_init("$file,nope", "aa", g)
        end
    end

    @testset "no such file" begin
        g = PottsGraph(6, 21; init=:rand)
        @test_throws ArgumentError PottsEvolver.CLI.load_init("does_not_exist.fasta", "aa", g)
    end

    @testset "explicit init-type overrides inference" begin
        mktempdir() do dir
            file = joinpath(dir, "root.fasta")
            seq = AASequence("ACDEFG")
            write_fasta(file, [seq]; labels=["root"])
            g = PottsGraph(6, 21; init=:rand)
            @test PottsEvolver.CLI.load_init(file, "aa", g) isa AASequence
            @test_throws ArgumentError PottsEvolver.CLI.load_init(file, "bogus", g)
        end
    end

    @testset "inference: unambiguous by content" begin
        g21 = PottsGraph(6, 21; init=:rand)
        g5 = PottsGraph(6, 5; init=:rand)

        mktempdir() do dir
            # amino-acid-exclusive letter forces AASequence
            f_aa = joinpath(dir, "aa.fasta")
            write_fasta(f_aa, [AASequence("ACDEFG")]; labels=["s"])
            @test PottsEvolver.CLI.load_init(f_aa, nothing, g21) isa AASequence

            # 'U' forces RNASequence
            f_rna = joinpath(dir, "rna.fasta")
            write_fasta(f_rna, [RNASequence("ACGUAC")]; labels=["s"])
            @test PottsEvolver.CLI.load_init(f_rna, nothing, g5) isa RNASequence

            # amino-acid-exclusive letter but incompatible graph -> error
            @test_throws ArgumentError PottsEvolver.CLI.load_init(f_aa, nothing, g5)
        end
    end

    @testset "inference: ambiguous content resolved by length/q" begin
        mktempdir() do dir
            # "ACGACG" is valid RNA (no exclusive letters) and, at L=6/q=21, also a valid
            # (if biologically odd) AASequence; codon would need length 18 (3*6).
            f = joinpath(dir, "ambiguous.fasta")
            write_fasta(f, [RNASequence("ACGACG")]; labels=["s"])

            g_aa = PottsGraph(6, 21; init=:rand)
            @test PottsEvolver.CLI.load_init(f, nothing, g_aa) isa AASequence

            g_rna = PottsGraph(6, 5; init=:rand)
            @test PottsEvolver.CLI.load_init(f, nothing, g_rna) isa RNASequence

            # length 18 over a q=21 graph must resolve to CodonSequence (6 codons)
            f_codon = joinpath(dir, "codon.fasta")
            open(f_codon, "w") do io
                write(io, ">s\n", "ACGACGACGACGACGACG\n")
            end
            @test PottsEvolver.CLI.load_init(f_codon, nothing, g_aa) isa CodonSequence
        end
    end

    @testset "inference: 'T' rules out RNA but not aa/codon" begin
        mktempdir() do dir
            f_t = joinpath(dir, "t2.fasta")
            open(f_t, "w") do io
                write(io, ">s\n", "ACGTAC\n")
            end
            g_aa = PottsGraph(6, 21; init=:rand)
            @test PottsEvolver.CLI.load_init(f_t, nothing, g_aa) isa AASequence

            g_rna = PottsGraph(6, 5; init=:rand)
            @test_throws ArgumentError PottsEvolver.CLI.load_init(f_t, nothing, g_rna)
        end
    end
end

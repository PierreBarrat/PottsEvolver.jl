@testset "GenCode: Symbol to Int to Symbol" begin
    for i in 1:length(codon_alphabet)
        aa_char = genetic_code(codon_alphabet[i])
        aa_int = genetic_code(i)
        @test (isnothing(aa_int) && aa_char == '*') || (aa_alphabet(aa_int) == aa_char)
    end
end

@testset "Genetic code" begin
    # This mainly tests self-consistency
    for aa in 1:length(aa_alphabet)
        @test aa in map(genetic_code, PottsEvolver.reverse_code(aa))
        @test genetic_code(PottsEvolver.reverse_code_rand(aa)) == aa

        aa_c = aa_alphabet(aa)
        @test aa_c in map(genetic_code, PottsEvolver.reverse_code(aa_c))
        @test genetic_code(PottsEvolver.reverse_code_rand(aa_c)) == aa_c
    end
end

@testset "Valid/Gap/Invalid codons" begin
    stop_codons = ["TAA", "TAG", "TGA"]
    for codon in symbols(codon_alphabet)
        # Check PottsEvolver.isstop
        @test in(prod(PottsEvolver.bases(codon)), stop_codons) == PottsEvolver.isstop(codon)
        # Check PottsEvolver.isvalid
        all_nt = all(in(['A', 'C', 'G', 'T']), PottsEvolver.bases(codon))
        all_gaps = all(==('-'), PottsEvolver.bases(codon))
        @test PottsEvolver.isvalid(codon) == (all_nt || all_gaps)
    end
end

@testset "Accessible codons" begin
    gap_codon = codon_alphabet(PottsEvolver.Codon("---"))
    codons = filter(PottsEvolver.iscoding, 1:length(codon_alphabet)) # Ints - only coding

    # Test the one argument version of accessible codons
    for codon in codons
        acc_codons = PottsEvolver.accessible_codons(codon)[1]
        acc_aas = PottsEvolver.accessible_codons(codon)[2]
        @test codon ∉ acc_codons
        @test gap_codon in acc_codons
        @test codon in PottsEvolver.accessible_codons(gap_codon)[1]
        @test genetic_code.(acc_codons) == acc_aas
    end
    @test gap_codon ∉ PottsEvolver.accessible_codons(gap_codon)[1]

    # Test the two argument version of accessible codons
    for codon in codons
        acc_codons = mapreduce(b -> PottsEvolver.accessible_codons(codon,b)[1], vcat, 1:3)
        acc_aas = mapreduce(b -> PottsEvolver.accessible_codons(codon,b)[2], vcat, 1:3)
        @test codon in acc_codons
        @test gap_codon ∉ acc_codons
        @test genetic_code.(acc_codons) == acc_aas
    end
    for b in 1:3
        @test isnothing(PottsEvolver.accessible_codons(gap_codon, b)[1])
    end

    # Check consistency between base-specific and total accessibility
    # only expected to work for coding (non-gap) codons
    for codon in codons
        # two argument version - filter self out
        acc_cod_1 = mapreduce(vcat, 1:3) do b
            PottsEvolver.accessible_codons(codon, b)[1]
        end
        acc_cod_1 = filter(!=(codon), acc_cod_1)
        acc_cod_1 = sort(acc_cod_1)

        # one argument version - filter gap out
        acc_cod_2 = filter(PottsEvolver.accessible_codons(codon)[1]) do c
            !PottsEvolver.isgap(codon_alphabet(c)) # need to filter gap out
        end |> sort

        @test acc_cod_1 == acc_cod_2
    end
end

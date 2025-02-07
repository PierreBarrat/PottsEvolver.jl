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
        @test in(prod(bases(codon)), stop_codons) == PottsEvolver.isstop(codon)
        # Codons made of only nucleotides are always valid
        @test PottsEvolver.isvalid(codon) ==
            (all(in(['A', 'C', 'G', 'T']), bases(codon)) || all(==('-'), bases(codon)))
    end
end

@testset "Accessible codons" begin
    # Check consistency between base-specific and total accessibility
    # only expected to work for coding (non-gap) codons
    codons = map(codon_alphabet, filter(PottsEvolver.iscoding, symbols(codon_alphabet)))
    for codon in codons
        acc_cod_1 = mapreduce(vcat, 1:3) do b
            PottsEvolver.accessible_codons(codon, b)[1]
        end |> sort
        acc_cod_2 = filter!(PottsEvolver.accessible_codons(codon)[1]) do c
            !PottsEvolver.isgap(codon_alphabet(c)) # need to filter gap out
        end |> sort
        @test acc_cod_1 == acc_cod_2

        acc_aa_1 = mapreduce(vcat, 1:3) do b
            PottsEvolver.accessible_codons(codon, b)[2]
        end |> sort
        acc_aa_2 = filter!(PottsEvolver.accessible_codons(codon)[2]) do c
            !PottsEvolver.isgap(c, aa_alphabet)
        end |> sort
        @test acc_aa_1 == acc_aa_2
    end
end

#====================================#
############# AASequence #############
#====================================#

@testset "AASequence Tests" begin
    # Test constructing the object
    aa_seq = AASequence([1, 2, 3, 4, 5])
    @test aa_seq isa AASequence
    @test length(aa_seq) == 5

    # Test indexing and setting the index in a sequence
    @test aa_seq[2] == 2
    aa_seq[2] = 10
    @test aa_seq[2] == 10

    # Test copying the object
    aa_seq_copy = copy(aa_seq)
    @test aa_seq_copy == aa_seq
    @test pointer(aa_seq.seq) != pointer(aa_seq_copy.seq)

    # Test in place copy
    aa_seq_2 = AASequence([5, 4, 3, 2, 1])
    copy!(aa_seq_copy, aa_seq_2)
    @test aa_seq_copy == AASequence([5, 4, 3, 2, 1])
    @test pointer(aa_seq_copy.seq) != pointer(aa_seq_2.seq)
end

#=============================================#
################ CodonSequence ################
#=============================================#

@testset "CodonSequence" begin
    @testset "From aa" begin
        aas = collect(1:6) # Int
        codon_seq = CodonSequence(aas; source=:aa)

        @test codon_seq isa CodonSequence{Int}
        @test PottsEvolver.translate(codon_seq).seq == codon_seq.aaseq
    end

    @testset "From codons" begin
        # Create an example codon sequence with matching amino acids
        codons = Int8.([1, 2, 3, 6, 7])
        aas = map(genetic_code, codons)
        codon_seq = CodonSequence(codons; source=:codon)

        # Test constructing the object
        @test codon_seq isa CodonSequence{Int8}
        @test length(codon_seq) == 5

        # Test indexing and setting the index in a sequence
        @test codon_seq[3] == 3
        codon_seq[3] = 7
        @test codon_seq[3] == 7
        @test codon_seq.aaseq[3] == genetic_code(7)

        # Test copying the object
        codon_seq_copy = copy(codon_seq)
        @test codon_seq_copy == codon_seq
        @test pointer(codon_seq.seq) != pointer(codon_seq_copy.seq)
        @test pointer(codon_seq.aaseq) != pointer(codon_seq_copy.aaseq)

        # Test in place copy
        codon_seq_2 = CodonSequence(Int8.([7, 6, 3, 2, 1]); source=:codon)
        copy!(codon_seq_copy, codon_seq_2)
        @test codon_seq_2.seq == Int8.([7, 6, 3, 2, 1])
        @test codon_seq_copy.seq == codon_seq_2.seq
        @test codon_seq_copy.aaseq == codon_seq_2.aaseq
        @test pointer(codon_seq_copy.seq) != pointer(codon_seq_2.seq)

        # Test fail in place copy when different lengths
        codon_seq_3 = CodonSequence(Int8.([7, 6, 3, 2]); source=:codon)
        @test_throws ArgumentError copy!(codon_seq_copy, codon_seq_3)
    end
end

#============================================================#
##################### Numerical sequence #####################
#============================================================#

@testset "NumSequence Tests" begin
    # Test constructing the object
    L, q = (10, 5)
    T = Int32

    num_seq = NumSequence(L, q; T)
    @test num_seq isa NumSequence{T,T(q)}
    @test length(num_seq) == L

    # Test indexing and setting the index in a sequence
    @test num_seq[4] >= 1 && num_seq[4] <= q
    num_seq[4] = 2
    @test num_seq[4] == 2

    # Test copying the object
    num_seq_copy = copy(num_seq)
    @test num_seq_copy == num_seq
    @test pointer(num_seq.seq) != pointer(num_seq_copy.seq)
    num_seq_copy[4] = 1 # num_seq[4] == 2 from test above
    @test num_seq_copy != num_seq

    # Test manually copying
    num_seq_copy = NumSequence{T,q}(copy(PottsEvolver.sequence(num_seq)))
    @test num_seq_copy == num_seq

    # Test in place copy
    num_seq_2 = NumSequence{T, T(q)}(ones(Int, L))
    num_seq_3 = NumSequence(ones(Int, L), 4) # will be NumSequence{Int, 4}, but q above is 5
    copy!(num_seq_copy, num_seq_2)
    @test num_seq_copy == num_seq_2
    @test pointer(num_seq_copy.seq) != pointer(num_seq_2.seq)
    @test_throws MethodError copy!(num_seq_3, num_seq_2)

    # Test failures
    @test_throws ArgumentError NumSequence([1, 2, 3]) # need to provide q
    @test_throws ArgumentError NumSequence([1, 2, 3], 2)
    @test_throws ArgumentError NumSequence([0, 1, 2], 21)
end

@testset "Sequences to alignent" begin
    M, L, q = (3, 10, 5)
    numseqs = [NumSequence(L, q; T=Int32) for _ in 1:M]
    aaseqs = [AASequence(L) for _ in 1:M] # default integer type should be IntType
    codonseqs = [CodonSequence(L) for _ in 1:M]

    A = Alignment(numseqs)
    @test A isa Alignment{Nothing,Int32}
    @test isnothing(A.alphabet)

    A = Alignment(aaseqs)
    @test A isa Alignment{Char,PottsEvolver.IntType}
    @test A.alphabet == aa_alphabet

    A = Alignment(codonseqs; as_codons=true)
    @test A isa Alignment{PottsEvolver.Codon,PottsEvolver.IntType}
    @test A.alphabet == codon_alphabet

    A = Alignment(codonseqs; as_codons=false)
    @test A isa Alignment{Char,PottsEvolver.IntType}
    @test A.alphabet == aa_alphabet
end

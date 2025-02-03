@testset "Energy differences" begin
    L = 5
    @testset "AA sequences" begin
        q = 21
        g = PottsGraph(L, q; init=:rand)

        a = AASequence(L)
        a.seq[1] = 1
        ref_ΔE = PottsEvolver.compute_energy_differences(a, g)

        b = copy(a)
        i = 1
        x = rand(2:q)
        b.seq[i] = x
        ΔE_true = PottsEvolver.compute_energy_differences(b, g)

        ΔE_buffer = similar(ref_ΔE)
        ΔE_test = PottsEvolver.compute_energy_differences(ref_ΔE, a, i, x, g; ΔE=ΔE_buffer)

        @test all(ΔE_test .≈ ΔE_true)
    end

    @testset "Codon sequences" begin
        q = length(codon_alphabet)
        g = PottsGraph(L, q; init=:rand)

        a = CodonSequence(L)
        a.seq[1] = 1
        a.aaseq[1] = 1
        ref_ΔE = PottsEvolver.compute_energy_differences(a, g)

        b = copy(a)
        i = 1
        x = rand(2:q)
        b.seq[i] = x
        b.aaseq[i] = genetic_code(x)
        ΔE_true = PottsEvolver.compute_energy_differences(b, g)

        ΔE_buffer = copy(ref_ΔE)
        ΔE_test = PottsEvolver.compute_energy_differences(ref_ΔE, a, i, x, g; ΔE=ΔE_buffer)

        @test all(ΔE_test .≈ ΔE_true)
    end
end







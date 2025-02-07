import TreeTools.Generate: balanced_binary_tree

@testset "prepare_tree" begin
    tree = balanced_binary_tree(8, 1.0)
    rootseq = CodonSequence(5)
    ptree = PottsEvolver.prepare_tree(tree, rootseq)

    @test TreeTools.share_labels(tree, ptree)
    @test ptree isa Tree{PottsEvolver.Sequence{typeof(rootseq)}}
    labels = map(label, nodes(tree; skiproot=true))
    for l in labels
        @test branch_length(tree[l]) ≈ branch_length(ptree[l])
    end
    @test ismissing(branch_length(root(ptree)))
    @test data(root(ptree)).seq == rootseq
    @test data(root(ptree)).seq !== rootseq
end

@testset "branch length meaning" begin
    L = 4
    b1 = 2
    b2 = 2.1 # rtol is 1e-6 for int
    b3 = 2.5

    p = BranchLengthMeaning(:step, :exact)
    @test PottsEvolver.steps_from_branchlength(b1, p, L) == b1
    @test_throws ArgumentError PottsEvolver.steps_from_branchlength(b2, p, L)
    @test_throws ArgumentError PottsEvolver.steps_from_branchlength(b3, p, L)

    p = BranchLengthMeaning(:step, :round)
    @test PottsEvolver.steps_from_branchlength(b1, p, L) == b1
    @test PottsEvolver.steps_from_branchlength(b2, p, L) == round(Int, b2)
    @test PottsEvolver.steps_from_branchlength(b3, p, L) == round(Int, b3)

    p = BranchLengthMeaning(:step, :poisson)
    @test PottsEvolver.steps_from_branchlength(b1, p, L) isa Int
    @test PottsEvolver.steps_from_branchlength(b2, p, L) isa Int
    @test PottsEvolver.steps_from_branchlength(b3, p, L) isa Int

    p = BranchLengthMeaning(:sweep, :exact)
    @test PottsEvolver.steps_from_branchlength(b1, p, L) == L * b1
    @test_throws ArgumentError PottsEvolver.steps_from_branchlength(b2, p, L)
    @test PottsEvolver.steps_from_branchlength(b3, p, L) ≈ L * b3

    p = BranchLengthMeaning(:sweep, :round)
    @test PottsEvolver.steps_from_branchlength(b1, p, L) == L * b1
    @test PottsEvolver.steps_from_branchlength(b2, p, L) == L * round(Int, b2)
    @test PottsEvolver.steps_from_branchlength(b3, p, L) ≈ L * b3
end

@testset "main function" begin
    L = 3
    g = PottsGraph(L, 21; init=:rand)
    # tree with branches of length one mut
    rootseq = AASequence(L)
    tree = PottsEvolver.prepare_tree(balanced_binary_tree(8, 1.0), rootseq)
    params = SamplingParameters(; step_meaning=:changed, Teq=5, burnin=100) # Teq & burnin should not matter!

    tree = PottsEvolver.mcmc_sample_tree!(g, tree, params)
    for node in nodes(tree; skiproot=true)
        nseq = data(node).seq
        aseq = data(ancestor(node)).seq
        @test hamming(nseq, aseq; normalize=false) == 1
    end
end

@testset "access functions" begin
    L, q = 10, 21
    g = PottsGraph(L, q; init=:rand)
    tree = TreeTools.Generate.balanced_binary_tree(8, 1.0)
    params = SamplingParameters(; step_meaning=:changed, Teq=1, burnin=100) # burnin should not matter!

    # providing root sequence
    rootseq = CodonSequence(L)
    sampled_tree = PottsEvolver.mcmc_sample_tree(g, tree, rootseq, params)
    @test tree !== sampled_tree # should be a copy
    @test data(root(sampled_tree)).seq == rootseq
    @test data(root(sampled_tree)).seq !== rootseq

    # providing init kwarg and burnin
    rootseq = PottsEvolver.NumSequence(L, q)
    sampled_tree = PottsEvolver.mcmc_sample_tree(g, tree, params; init=rootseq)
    @test tree !== sampled_tree # should be a copy
    @test data(root(sampled_tree)).seq != rootseq # finite burnin: with probability ≈1 they should be different

    # providing init kwarg and no burnin
    params = SamplingParameters(; step_meaning=:accepted, Teq=0, burnin=0) # burnin should not matter!
    sampled_tree = @test_logs min_level = Logging.Warn PottsEvolver.mcmc_sample_tree(
        g, tree, params; init=:random_aa
    )
end

@testset "index swap" begin
    L, q = 10, 21
    M = 5
    g = PottsGraph(L, q; init=:rand)
    tree = balanced_binary_tree(8, 1.0)
    params = SamplingParameters(; step_meaning=:changed, Teq=0)

    # Test for 'raw' sequence output
    data_pertree = mcmc_sample(g, tree, M, params; init=:random_aa, alignment_output=false)
    data_pernode = PottsEvolver.pernode_alignment(data_pertree)
    for nlabel in map(label, nodes(tree)), m in 1:M
        # m represents the tree, *i.e.* the sampling realization
        @test if isleaf(tree[nlabel])
            data_pertree[m].leaf_sequences[nlabel] == data_pernode.sequences[nlabel][m]
        else
            data_pertree[m].internal_sequences[nlabel] == data_pernode.sequences[nlabel][m]
        end
    end

    # for alignment output
    data_pertree = mcmc_sample(
        g, tree, M, params; init=:random_codon, alignment_output=true
    )
    data_pernode = PottsEvolver.pernode_alignment(data_pertree)
    for nlabel in map(label, nodes(tree)), m in 1:M
        # m represents the tree, *i.e.* the sampling realization
        seq_left = if isleaf(tree[nlabel])
            find_sequence(nlabel, data_pertree[m].leaf_sequences)[2]
        else
            find_sequence(nlabel, data_pertree[m].internal_sequences)[2]
        end
        @test seq_left == data_pernode.sequences[nlabel][m]
    end
end

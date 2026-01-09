```@meta 
DocTestSetup = quote
    using PottsEvolver
    L = 20
    q = 21
    g_null = PottsGraph(L, q, Float64; init=:null); # null initialization
    g_rand = PottsGraph(L, q, Float64; init=:rand); # random initialization
end
```

# PottsGraph

The `PottsGraph` structure represents a Potts model, a statistical model used to describe interactions in a system with multiple states (e.g., amino acids in a protein sequence). It is defined by the following fields:

- `J`: A 4-dimensional array of dimensions `q x q x L x L` representing the coupling parameters between states at different positions.
- `h`: A 2-dimensional array of dimensions `q x L` representing the local field parameters for each state at each position.
- `β`: The inverse temperature parameter, used during sampling. 
- `alphabet`: An optional `Alphabet` object to map states to characters (e.g., amino acids), see the [BioSequenceMappings.jl](https://github.com/PierreBarrat/BioSequenceMappings.jl) package.

`q` represents the number of states (e.g., 21 for amino acids + gap), and `L` is the length of the sequence.

## Structure and Fields

The `PottsGraph` is defined as follows:

```julia
@kwdef mutable struct PottsGraph{T<:AbstractFloat}
    J::Array{T,4}
    h::Array{T,2}
    β::T = 1.0
    alphabet::Union{Nothing, BioAequenceMappings.Alphabet{Char,<:Integer}} = aa_alphabet
end
```

### Creating a PottsGraph

You can create a `PottsGraph` with null or random initialization:

```@repl pottsgraph_1
using PottsEvolver
L, q = (10, 21);
g_null = PottsGraph(L, q, Float64; init=:null); # null initialization
g_null.J[:, :, 1, 2]
g_rand = PottsGraph(L, q, Float64; init=:rand); # random initialization
size(g_rand)
g_rand.J[:, :, 1, 2] # coupling matrix between sites 1 and 2
```
It is sometimes easier to view the couplings `J` as a square matrix. 
```jldoctest pottsgraph_1
julia> J_matrix = PottsEvolver.couplings_as_matrix(g_rand.J);

julia> size(J_matrix) # (L*q x L*q)
(420, 420)

julia> J_matrix ≈ J_matrix' # The coupling matrix should be symetric
true

julia> J_tensor = PottsEvolver.matrix_as_couplings(J_matrix, L, q); # back to tensor

julia> J_tensor ≈ g_rand.J
true
```
!!! note
    There is no check that the coupling matrix is symetric, although it should normally be the case for Potts models. 

## Reading and writing `PottsGraph`

### Reading from a File

The `read_graph` function reads a `PottsGraph` from a file. The file should follow the format 
```
J i j a b value # coupling lines
...
h i a value # field lines
```
where `i` and `j` are sequence positions and `a`, `b` are integer states (_e.g._ representing amino acids). 
```@repl pottsgraph_1
# Read the graph in the example folder
g = read_graph("../../example/toy_potts.dat)
```

### Writing to a File

The `write` function writes the parameters of a `PottsGraph` to a file.
```julia
write("path/to/output.txt", g; sigdigits=5, index_style=1)
```

## Basic Functions

### Energy of a Sequence

The `energy` function calculates the statistical energy of a sequence. 

```@repl pottsgraph_1
seq = AASequence(map(aa_alphabet, collect("ACDEFGHIKLMNPQRSTVWY")))
E = energy(seq, g_rand)
```

### Gauge Change

The `set_gauge!` function allows you to change the gauge of the `PottsGraph`. Supported gauges include `:zero_sum` and `:lattice_gas`.

```jldoctest pottsgraph_1
julia> PottsEvolver.set_gauge!(g_rand, :zero_sum);

julia> all(≈(0; atol=1e-14), sum(g_rand.h, dims=1)) # sum of h values along each sequence position is zero
true

julia> all(≈(0; atol=1e-14), sum(g_rand.J, dims=[1,2])) #same goes for J
true
```

## Additional Utilities

### Profile Model

The `profile_model` function creates a `PottsGraph` with only fields that fit the single-site frequencies in a given matrix.

```jldoctest pottsgraph_1
julia> f1 = rand(21, 10);  # Example frequency matrix - does not have to be normalized

julia> profile = PottsEvolver.profile_model(f1; pc=1e-2);

julia> all(≈(0; atol=1e-14), profile.J)
true
```

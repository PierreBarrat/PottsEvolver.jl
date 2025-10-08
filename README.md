# PottsEvolver

[![Build Status](https://github.com/PierreBarrat/PottsEvolver.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/PierreBarrat/PottsEvolver.jl/actions/workflows/CI.yml?query=branch%3Amaster)

In progress.

Simulate protein evolution using a Potts model and MCMC sampling. 
Can simulate either a linear chain, or along branches of a tree. 

To learn more, check out the example Pluto notebook at `examples/example_notebook.jl`. 

## Running the example notebook 

- Install [julia](https://julialang.org/install/)
- Launch julia and install [Pluto](https://plutojl.org/) by running
  ```
  using Pkg
  Pkg.add("Pluto")
  ```
- start Pluto from within julia by running
  ```
  using Pluto
  Pluto.run()
  ```
  Pluto should start in a browser window.
- Open the notebook from within Pluto. 

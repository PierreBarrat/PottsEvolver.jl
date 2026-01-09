# PottsEvolver

[![Build Status](https://github.com/PierreBarrat/PottsEvolver.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/PierreBarrat/PottsEvolver.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PierreBarrat.github.io/PottsEvolver.jl/dev/)
<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PierreBarrat.github.io/PottsEvolver.jl/stable/) -->

Simulate protein evolution using a Potts model and MCMC sampling. 
Can simulate either a single lineage sampled at given times, or along branches of a tree. 

This is the code used for the paper. 
```
Generative continuous time model reveals epistatic signatures in protein evolution
Andrea Pagnani & Pierre Barrat-Charlaix
biorXiv, https://doi.org/10.1101/2025.09.17.676821
```
Please cite it if you use the package in your work. 

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

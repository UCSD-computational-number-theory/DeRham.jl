# DeRham.jl

[![CI](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml)

This package aims to implement cutting-edge-speed algorithms that calculate the cohomology of algebraic varieties and the extra structure on them. Longer term goals include having fast versions of all the various forms of Kedlaya's algorithm, Harvey's algorithm, Moonen's algorithm, etc.

**This package is a work-in-progress.** PRs and contributions are welcome. 

## Current functionality

The implementation that is currently available is Kedlaya's algorithm that computes the zeta functions of projective hypersurfaces using a method known as controlled reduction (developed by Costa, Harvey, and Kedlaya), without the nondegeneracy optimization (see Proposition 1.15 [here](https://edgarcosta.org/assets/articles/EdgarCosta-PhDthesis.pdf)).

### Getting Started

Dependencies: Need to install Julia and Oscar [ToDo: Add installation guides]. Within Julia, require the BigIntegers, LinearAlgebra, and Combinatorics package [ToDo: Add more dependencies].

### Sample code

Here is an example for how to use the code. Execute the following lines in Julia REPL:  
```
include("DeRham.jl")
p = 7
R, (x,y,z) = polynomial_ring(GF(p), ["x","y","z"])
f = y^2*z - x^3 - x*z^2 - z^3
DeRham.zeta_function(f,fastevaluation=true,algorithm=:naive)
```

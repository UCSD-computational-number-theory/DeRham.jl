# DeRham.jl

This package aims to implement cutting-edge-speed algorithms that calculate the cohomology of algebraic varieties and the extra structure on them. Longer term goals include having fast versions of all the various forms of Kedlaya's algorithm, Harvey's algorithm, Moonen's algorithm, etc.

**This package is a work-in-progress.** PRs and contributions are welcome. 

## Current functionality

All that is currently implemented is Kedlaya's algorithm with controlled reduction, without the nondegeneracy optimization (see Proposition 1.15 [here](https://edgarcosta.org/assets/articles/EdgarCosta-PhDthesis.pdf)).

### Getting Started

Dependencies: Need to install Julia and Oscar [ToDo: Add installation guides]. Within Julia, require the BigIntegers, LinearAlgebra, and Combinatorics package.

### Sample code

Here is an example for how to use the code. Execute the following lines in Julia REPL:  
```
include("ZetaFunction.jl")
using Oscar
n = 2
d = 3
p = 7
R = GF(p,1)
PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
x,y,z = Vars
f = y^2*z - x^3 - x*z^2 - z^3
ZetaFunction.computeAll(n,d,f,7,p,R,PR,Vars)
```

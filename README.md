# Griffiths-Dwork construction implementation

Given a homogenous, irreducible polynomial in $\mathbb{Q}_p[x_1, ... x_n]$, this script will generate the basis vectors for the $n$-th de Rham cohomology of $(X, Z)/\mathbb{Q}_p$ by following the Griffiths-Dwork construction algorithm.

Dependencies: Need to install Julia and Oscar [ToDo: Add installation guides]. Within Julia, require the BigIntegers, LinearAlgebra, and Combinatorics package.

Here is an example for how to use the code. Execute the following lines in command line:  
```
julia
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
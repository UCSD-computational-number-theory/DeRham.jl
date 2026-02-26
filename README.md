
# DeRham.jl

[![CI](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/UCSD-computational-number-theory/DeRham.jl/actions/workflows/CI.yml)

This package aims to implement cutting-edge algorithms that calculate

* point counts of equations mod $p$ (algebraic varieties over $\mathbb{F}_p$)
* the cohomology of algebraic varieties
* the extra structure on this cohomology

Long term goals include having fast versions of all the various forms of Kedlaya's algorithm, Harvey's algorithm, Harvey's Trace Formula, Moonen's algorithm, etc.

**This package is work-in-progress research software.** It may have rough edges. PRs and contributions are welcome.

## Current functionality

We currently implement a variant of Kedlaya's algorithm that computes the zeta functions (and consequently the newton polygons) of smooth projective hypersurfaces of degree $$d$$ over $$\mathbb{F}_p$$ (for $$p\nmid d$$) using a method known as controlled reduction, see Proposition 1.15 [here](https://edgarcosta.org/assets/articles/EdgarCosta-PhDthesis.pdf).

### Getting Started

#### Installing for the first time

First, install Julia 1.12 and Oscar. If you're new to Julia and Oscar, you can find a tutorial for non-experts [here](https://jjgarzella.github.io/blog/how-to-install-julia/). Then, clone this repo with `git clone https://github.com/UCSD-computational-number-theory/DeRham.jl.git`. 
If you'd like to use CUDA, you'll need to install the CUDA driver, see [the instructions in the CUDA.jl docs](https://cuda.juliagpu.org/stable/installation/overview/).

`Revise.jl` is recommended for hot reloading. From a new Julia REPL, run
`using Pkg; Pkg.add("Revise")` (or, using Julia's package prompt, run `] add Revise`) once, and then run `using Revise` upon opening every new REPL.

The first time you install DeRham.jl, it will take a while for Julia to download and compile all of the dependincies. To do this, make sure your Julia REPL is in the DeRham.jl project folder. Then, in the package prompt (i.e. after pressing `]`), run

```
activate .
add https://github.com/UCSD-computational-number-theory/GPUFiniteFieldMatrices.jl.git
instantiate
```

The `instantiate` step will take a while as Julia's package manager downloads and compiles all dependencies.

Then, press backspace to go back to a julia prompt and run

```
using DeRham, Oscar
```

#### Starting a new session after installation is complete

After starting a new Julia REPL, run

```
using Revise
] activate .
using DeRham, Oscar
```

### Sample code

#### Zeta Functions and Newton polygons

Execute the following lines in Julia REPL:

```julia
p = 7
R, (x,y,z) = GF(p)[:x,:y,:z]
f = y^2*z - x^3 - x*z^2 - z^3
DeRham.zeta_function(f)
DeRham.newton_polygon(f)

g = DeRham.random_hypersurface(4,3,p)
DeRham.zeta_function(g,givefrobmat=true)

h = DeRham.fermat_hypersurface(4,4,p)
DeRham.newton_polygon(g)
```

#### Nondegeneracy Conditions (advanced)

We support using the nondegeneracy condition in [Costa's Thesis](https://edgarcosta.org/assets/articles/EdgarCosta-PhDthesis.pdf) called "S-smoothness", where "S" is a subset of `{0, ..., n-1}` (where `n` is the number of variables in the polynomial ring R).

By default, S is the full set `{0, ..., n-1}`, which is equivalent to the smoothness of the hypersufaces. Most users will not have
to worry about this.

There are two situations (currently) which require the use of S. First, if the degree `d` of the hypersurface is less than the number of varibles `n`. In this case, S can have size at most 3.

```julia
p = 7
R, (x,y,z,w) = GF(p)[:x,:y,:z]
f = x^4 + y^4 + z^4 + w^4 + x*y*z*w

DeRham.zeta_function(f,S=[0,1,2])
```

Secondly, if one uses that `:varbyvar` reduction policy, then S must be `{n}`

```julia
p = 11
R, (x,y,z,w) = GF(p)[:x,:y,:z]
f = x^4 + y^4 + z^4 + w^4 + x*y*z*w

DeRham.zeta_function(f,S=[n], algorithm=:varbyvar)
```

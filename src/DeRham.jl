module DeRham

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

verbose = false 

include("Utils.jl")
include("FindMonomialBasis.jl")
include("PrecisionEstimate.jl")
include("SmallestSubsetSmooth.jl")
include("PolynomialWithPole.jl")

include("StandardReduction.jl")

include("ControlledReduction.jl")
include("Frobenius.jl")
include("FinalReduction.jl")

include("ZetaFunction.jl")

# TODO: export Zeta Function functions, but first we probably want to rename compute_all

end

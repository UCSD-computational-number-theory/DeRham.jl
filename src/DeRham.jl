module DeRham

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

verbose = true

include("Utils.jl")
include("FindMonomialBasis.jl")
include("SlopesPolygon.jl")
include("PolynomialWithPole.jl")
include("PrecisionEstimate.jl")
include("SmallestSubsetSmooth.jl")

include("StandardReduction.jl")

include("ControlledReduction.jl")
include("Frobenius.jl")
include("FinalReduction.jl")

include("CharPolyFrob.jl")


include("ZetaFunction.jl")

# TODO: export Zeta Function functions

end

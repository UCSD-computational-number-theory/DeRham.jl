module DeRham

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics
using Memoize

#include("NemoAdditions.jl")

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

include("ExamplePolynomials.jl")

include("ZetaFunction.jl")
include("PointCounts.jl")

# TODO: export Zeta Function functions

end

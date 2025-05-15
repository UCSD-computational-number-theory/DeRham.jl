module DeRham

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics
using Memoize
using LRUCache
using LFUDACache
using OhMyThreads

using CUDA
using GPUFiniteFieldMatrices
# Pkg.add(url="https://github.com/UCSD-computational-number-theory/GPUFiniteFieldMatrices.jl")

# comment this out when not debugging
#using Debugger

#include("NemoAdditions.jl")

include("Utils.jl")
include("GradedExpCache.jl")
include("LinearAlgebraWrappers.jl")
include("FindMonomialBasis.jl")
include("SlopesPolygon.jl")
include("PolynomialWithPole.jl")
include("PrecisionEstimate.jl")
include("SmoothNondegenerate.jl")

include("StandardReduction.jl")

include("PEP.jl")
include("EvaluatePEP.jl")
include("ComputeRuv.jl")
include("ControlledReduction.jl")

include("Frobenius.jl")
include("FinalReduction.jl")
include("CharPolyFrob.jl")

include("ExamplePolynomials.jl")

include("ZetaFunction.jl")
include("PointCounts.jl")
include("NewtonPolygon.jl")
include("K3OverF3.jl")

# TODO: export Zeta Function functions

end

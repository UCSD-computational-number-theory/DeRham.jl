module DeRham

using Oscar:
    Oscar,
    GF,
    MPolyBuildCtx,
    PadicField,
    QQ,
    RingElem,
    ZZ,
    ZZMatrix,
    ZZModMatrix,
    ZZRingElem,
    base_ring,
    basis,
    change_base_ring,
    characteristic,
    charpoly,
    coeff,
    coefficient_ring,
    coefficients,
    data,
    derivative,
    divexact,
    evaluate,
    exponent_vector,
    factor,
    finish,
    gens,
    grade,
    hom,
    homogenizer,
    ideal,
    is_prime,
    is_smooth,
    lift,
    map_coefficients,
    map_entries,
    matrix,
    matrix_space,
    modulus,
    ncols,
    nrows,
    number_of_rows,
    nvars,
    polynomial_ring,
    power_series_ring,
    proj,
    push_term!,
    quo,
    residue_ring,
    rref,
    special_linear_group,
    terms,
    total_degree,
    valuation,
    vars,
    zero!,
    zero_matrix,
    zzModMatrix
using BitIntegers: BitIntegers
using LinearAlgebra: LinearAlgebra, I, nullspace, rank
using Combinatorics: Combinatorics
using Memoize: Memoize
using LRUCache: LRUCache, LRU
using LFUDACache: LFUDACache, LFUDA
using OhMyThreads: OhMyThreads

using CUDA: CUDA, @cuda, blockDim, blockIdx, threadIdx, totalmem
using GPUFiniteFieldMatrices:
    GPUFiniteFieldMatrices,
    CuModArray,
    CuModMatrix,
    CuModVector,
    KaratsubaArray,
    KaratsubaMatrix,
    KaratsubaVector
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
include("ControlledReductionCache.jl")

include("ZetaFunctionParams.jl")
include("PEP.jl")
include("EvaluatePEP.jl")
include("ComputeRuv.jl")
include("ControlledReduction.jl")
include("PChunk.jl")
include("VarByVar.jl")
include("DepthFirst.jl")
include("AKR.jl")

include("Frobenius.jl")
include("FinalReduction.jl")
include("CharPolyFrob.jl")

include("ExamplePolynomials.jl")

include("PrecisionInformation.jl")
include("K3OverF3.jl")

include("./invariants/ZetaCoefficients.jl")
include("./invariants/HodgePolygon.jl")
include("./invariants/Smoothness.jl")
include("./invariants/NewtonPolygon.jl")
include("./invariants/HasseWittMatrix.jl")
include("./invariants/aNumber.jl")
include("./invariants/FrobeniusMatrix.jl")
include("./invariants/CohomologyBasis.jl")
include("./invariants/PointCounts.jl")
include("./invariants/LPolynomial.jl")
include("./invariants/ArtinMazurHeight.jl")


# TODO: export Zeta Function functions

end

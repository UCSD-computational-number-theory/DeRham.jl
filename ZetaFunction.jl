module ZetaFunction 

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("AutomatedScript.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")

"""
    computeFrobeniusMatrix(n,d,Reductions,T)

Computes Frobenius Matrix

INPUTS: 
* "n" -- dimension of ambient projective space 
* "d" -- degree of the polynomial f2
* "Reductions" -- output of computeReductionOfTransformLA
* "T" -- output of computeT
"""
function computeFrobeniusMatrix(n,d,Reductions,T)
    FrobMatTemp = []
    denomArray = []
    for i in 1:length(Reductions)
        push!(denomArray, QQ(lift(ZZ,Reductions[i][3])))
        #push!(denomArray,lift(ZZ,Factorial(PrecisionRing(p*(Basis[i][2]+N-1)-1),PrecisionRing(1))/(p^(n-1))))
        push!(FrobMatTemp,T*transpose(AutomatedScript.convert_p_to_m([Reductions[i][1]],AutomatedScript.gen_exp_vec(n+1,d*n-n-1))))
        #(p^(n-1)/Factorial(PrecisionRing(p*(Basis[i][2]+N-1)-1),PrecisionRing(1)))
    end
    FrobMat = hcat(FrobMatTemp...)
    #MS = matrix_space(ZZ,nrows(FrobMat),ncols(FrobMat))
    #FM = MS()
   
    @assert nrows(FrobMat) == ncols(FrobMat) "Frobenius matrix is not a square matrix."
    R = matrix_space(QQ, nrows(FrobMat), ncols(FrobMat))
    FM = Array{QQFieldElem}(undef, nrows(FrobMat), ncols(FrobMat))
    for i in axes(FrobMat, 1)
        for j in axes(FrobMat, 2)
            FM[i,j] = lift(ZZ,FrobMat[i,j])/denomArray[i]
        end
    end

    return R(FM)
end

"""
    LPolynomial(FM, q)

Given the Frobenius matrix, computes the corresponding L-polynomial det(1-tq^{-1}FM)

INPUT: 
* "FM" -- Frobenius matrix 
* "q" -- size of the base field 
"""

function LPolynomial(FM, q)
    @assert size(FM, 1) == size(FM, 2) "FM is not a square matrix"

    P, T = polynomial_ring(QQ, "T")
    f = charpoly(P, FM)

    return reverse(f)
end 


"""
    computeAll(n, d, f, precision, p, R, PR, vars)

Wrapper function that outputs the Frobenius Matrix

INPUTS: 
* "n" -- number of variables - 1
* "d" -- degree of f
* "f" -- Oscar polynomial
* "precision" -- not in use
* "p" -- prime number
* "R" -- basefield(parent(f))
* "PR" -- parent(f)
* "vars" -- generators of PR
"""
function computeAll(n, d, f, precision, p, R, PR, var, verbose=false)
    #Nm = PrecisionEstimate.compute_precisions_each(p,precision,n)
    #N = max(Nm...)
    N = 6
    s = N + n - 1
    #M = Int(precision + floor((p*s-1)/(p-1) + 1))
    M = 15
    #PrecisionRing = PadicField(p,M)
    PrecisionRing = residue_ring(ZZ, p^M)
    println(typeof(PrecisionRing))
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f, R, PR)
    fLift = ControlledReduction.liftCoefficients(PrecisionRing, PrecisionRingPoly, f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp, ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    T = ControlledReduction.computeT(BasisTLift, fLift, n, d, PrecisionRing, PrecisionRingPoly)
    #S = SmallestSubsetSmooth.smallest_subset_s_smooth(fLift,n)
    S = [0,1,2]
    #S = []
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = ControlledReduction.applyFrobeniusToBasis(Basis, n, d, fLift, N, p, PrecisionRing, PrecisionRingPoly)
    psuedoInverseMatTemp = CopiedFindMonomialBasis.psuedo_inverse_controlled(f,S,R,PR)
    psuedoInverseMat = zeros(PrecisionRing, nrows(psuedoInverseMatTemp), ncols(psuedoInverseMatTemp))
    for i in 1:nrows(psuedoInverseMat)
        for j in 1:ncols(psuedoInverseMat)
            psuedoInverseMat[i,j] = PrecisionRing(lift(ZZ, psuedoInverseMatTemp[i,j]))
        end
    end
    Reductions = ControlledReduction.computeReductionOfTransformLA(FBasis, n, d, p, N, S, fLift, psuedoInverseMat, PrecisionRing, PrecisionRingPoly)
    FM = computeFrobeniusMatrix(n, d, Reductions, T) 

    if verbose
        println("The Frobenius matrix is $FM")
    end

    return LPolynomial(FM, p)
end

end 

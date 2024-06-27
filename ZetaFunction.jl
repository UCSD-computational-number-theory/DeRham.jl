module ZetaFunction 

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")
include("Frobenius.jl")
include("FinalReduction.jl")

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
    println(Reductions)
    for i in 1:length(Reductions)
        #push!(denomArray, QQ(lift(ZZ,Reductions[i][3])))
        #push!(denomArray,lift(ZZ,Factorial(PrecisionRing(p*(Basis[i][2]+N-1)-1),PrecisionRing(1))/(p^(n-1))))
        push!(FrobMatTemp,T*transpose(Utils.convert_p_to_m([Reductions[i][1]],Utils.gen_exp_vec(n+1,d*n-n-1))))
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
            println(lift(ZZ,FrobMat[i,j]))
            #println(denomArray[i])
            #FM[i,j] = lift(ZZ,FrobMat[i,j])/denomArray[i]
            FM[i,j] = lift(ZZ,FrobMat[i,j])
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
    if verbose
        println("Working with a degree $d hypersurface in P^$n")
    end 

    N = 2 # series precision 
    s = N + n - 1
    #M = Int(precision + floor((p*s-1)/(p-1) + 1))
    M = 3 # Absolute precision

    if verbose
        println("We work modulo $p^$M, and compute up to the $N-th term of the Frobenius power series")
    end 

    #PrecisionRing = PadicField(p,M)
    PrecisionRing, = residue_ring(ZZ, p^M)
    println(typeof(PrecisionRing))
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f, R, PR)
    num_BasisT = length(BasisT)

    if verbose
        println("There are $num_BasisT basis elements in H^$n")
    end 

    fLift = Utils.liftCoefficients(PrecisionRing, PrecisionRingPoly, f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp, Utils.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    T = FinalReduction.computeT(BasisTLift, fLift, n, d, PrecisionRing, PrecisionRingPoly)
    #S = SmallestSubsetSmooth.smallest_subset_s_smooth(fLift,n)
    S = [0,1,2]
    #S = []
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = Frobenius.applyFrobeniusToBasis(Basis, n, d, fLift, N, p, PrecisionRing, PrecisionRingPoly)
    pseudoInverseMatTemp = CopiedFindMonomialBasis.pseudo_inverse_controlled_lifted(f,S,R,PR,M)
    pseudoInverseMat = zeros(PrecisionRing, nrows(pseudoInverseMatTemp), ncols(pseudoInverseMatTemp))

    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    fLift = Utils.liftCoefficients(ZZ,PRZZ,f)
    controlledMatrixZZ = CopiedFindMonomialBasis.compute_controlled_matrix(fLift, d * n - n + d - length(S), S, ZZ, PRZZ)
    pseudoInverseMatModP = matrix(ZZ, [lift(ZZ,x) for x in Array(pseudoInverseMatTemp)])
    pseudoInverseMatNew = Utils.henselLift(p,M,controlledMatrixZZ, pseudoInverseMatModP)
    
    #for i in 1:nrows(pseudoInverseMat)
    #    for j in 1:ncols(pseudoInverseMat)
    #        pseudoInverseMat[i,j] = PrecisionRing(lift(ZZ, pseudoInverseMatTemp[i,j]))
    #    end
    #end
    Reductions = ControlledReduction.reducetransform_LA_descending(FBasis, n, d, p, N, S, fLift, pseudoInverseMat, PrecisionRing, PrecisionRingPoly)
    println(Reductions)
    FM = computeFrobeniusMatrix(n, d, Reductions, T)
    println(FM)

    if verbose
        println("The Frobenius matrix is $FM")
    end

    return LPolynomial(FM, p)
end

end 

#=
include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")
include("ZetaFunction.jl")
include("TestControlledReduction.jl")
n = 2
d = 3
p = 7
R = GF(p,1)
PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
x0,x1,x2 = Vars
f = x1^2*x2 - x0^3 - x0*x2^2 - x2^3
S = [0,1,2]
Test = ZetaFunction.computeAll(n,d,f,1,p,R,PR,var)
=#

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
    compute_frobenius_matrix(n,d,Reductions,T)

Computes Frobenius Matrix

INPUTS: 
* "n" -- dimension of ambient projective space 
* "d" -- degree of the polynomial f2
* "Reductions" -- output of computeReductionOfTransformLA
* "T" -- output of computeT
"""
function compute_frobenius_matrix(n, d, Reductions, T)
    R = parent(T[1,1])
    FrobMatTemp = []
    denomArray = []
    ev = Utils.gen_exp_vec(n+1,d*n-n-1,:invlex)
    VS = matrix_space(R,length(ev),1)
    for i in 1:length(Reductions)
        temp = VS()
        temp2 = Utils.convert_p_to_m([Reductions[i][1][1]],ev)
        for i in 1:length(ev)
            temp[i,1] = R(temp2[i])
        end
        #push!(denomArray, QQ(lift(ZZ,Reductions[i][3])))
        #push!(denomArray,lift(ZZ,Factorial(PrecisionRing(p*(Basis[i][2]+N-1)-1),PrecisionRing(1))/(p^(n-1))))
        push!(FrobMatTemp,T*temp)
        #(p^(n-1)/Factorial(PrecisionRing(p*(Basis[i][2]+N-1)-1),PrecisionRing(1)))
    end
    FrobMat = hcat(FrobMatTemp...)
    #MS = matrix_space(ZZ,nrows(FrobMat),ncols(FrobMat))
    #FM = MS()
    println(FrobMat)
   
    @assert nrows(FrobMat) == ncols(FrobMat) "Frobenius matrix is not a square matrix."
    R = matrix_space(R, nrows(FrobMat), ncols(FrobMat))
    FM = Array{ZZModRingElem}(undef, nrows(FrobMat), ncols(FrobMat))
    for i in axes(FrobMat, 1)
        for j in axes(FrobMat, 2)
            #println(lift(ZZ,FrobMat[i,j]))
            #println(denomArray[i])
            #FM[i,j] = lift(ZZ,FrobMat[i,j])/denomArray[i]
            #FM[i,j] = lift(ZZ,FrobMat[i,j])
            FM = FrobMat[i,j]
        end
    end

    return R(FM)
end

"""
    LPolynomial(FM, q)

Given the Frobenius matrix, computes the corresponding L-polynomial det(1-tq^{-1}FM)

INPUT: 
* "FM" -- Frobenius matrix 
"""

function LPolynomial(FM)
    @assert size(FM, 1) == size(FM, 2) "FM is not a square matrix"

    P, T = polynomial_ring(parent(FM[1,1]), "T")
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
function compute_all(f, precision, verbose=false)
    p = Int64(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    d = degree(f,1)
    PR = parent(f)
    R = coefficient_ring(parent(f))
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

    precisionring, = residue_ring(ZZ, p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])
    basis = CopiedFindMonomialBasis.compute_monomial_bases(f, R, PR) # basis of cohomology 
    num_BasisT = length(basis)

    if verbose
        println("There are $num_BasisT basis elements in H^$n")
    end 

    T = FinalReduction.computeT(f, basis, M)
    #S = SmallestSubsetSmooth.smallest_subset_s_smooth(fLift,n)
    S = [0,1,2]
    #S = []

    BasisTLift = []
    for i in basis
        temp = []
        for j in i
            push!(temp, Utils.liftCoefficients(precisionring,precisionringpoly,j))
        end
        push!(BasisTLift,temp)
    end

    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end

    fLift = Utils.liftCoefficients(precisionring, precisionringpoly, f)
    FBasis = Frobenius.applyFrobeniusToBasis(Basis, n, d, fLift, N, p, precisionring, precisionringpoly)
    l = d * n - n + d - length(S)
    pseudo_inverse_mat_new = CopiedFindMonomialBasis.pseudo_inverse_controlled_lifted(f,S,l,M)
    #pseudoInverseMat = zeros(PrecisionRing, nrows(pseudoInverseMatTemp), ncols(pseudoInverseMatTemp))

    #PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    #fLift = Utils.liftCoefficients(ZZ,PRZZ,f)
    #controlledMatrixZZ = CopiedFindMonomialBasis.compute_controlled_matrix(fLift, d * n - n + d - length(S), S, ZZ, PRZZ)
    #pseudoInverseMatModP = matrix(ZZ, [lift(ZZ,x) for x in Array(pseudoInverseMatTemp)])
    #pseudo_inverse_mat_new = Utils.henselLift(p,M,controlledMatrixZZ, pseudoInverseMatModP)
    
    #for i in 1:nrows(pseudoInverseMat)
    #    for j in 1:ncols(pseudoInverseMat)
    #        pseudoInverseMat[i,j] = PrecisionRing(lift(ZZ, pseudoInverseMatTemp[i,j]))
    #    end
    #end
    Reductions = ControlledReduction.reducetransform_LA_descending(FBasis, n, d, p, N, S, fLift, pseudo_inverse_mat_new, precisionring, precisionringpoly)
    println(Reductions)
    ev = Utils.gen_exp_vec(n+1,n*d-n-1,:invlex)
    println(Utils.convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev))
    FM = compute_frobenius_matrix(n, d, Reductions, T)
    println(FM)

    if verbose
        println("The Frobenius matrix is $FM")
    end

    return LPolynomial(FM)
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
Test = ZetaFunction.compute_all(f,3)
=#

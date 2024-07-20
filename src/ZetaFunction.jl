#module ZetaFunction 
#
#using Oscar
#using BitIntegers
#using LinearAlgebra
#using Combinatorics
#
#
#include("Utils.jl")
#include("FindMonomialBasis.jl")
#include("PrecisionEstimate.jl")
#include("SmallestSubsetSmooth.jl")
#include("PolynomialWithPole.jl")
#
#include("StandardReduction.jl")
#
#include("ControlledReduction.jl")
#include("Frobenius.jl")
#include("FinalReduction.jl")
#
#verbose = false

"""
    compute_frobenius_matrix(n,d,Reductions,T)

Computes Frobenius Matrix

INPUTS: 
* "n" -- dimension of ambient projective space 
* "d" -- degree of the polynomial f2
* "N" -- integer, series precision
* "Reductions" -- output of computeReductionOfTransformLA
* "T" -- output of computeT
"""
function compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, termorder)
    verbose && println("Terms after controlled reduction: $Reductions")
    R = parent(T[1,1])
    frob_mat_temp = []
    denomArray = []
    ev = gen_exp_vec(n+1,d*n-n-1,termorder)
    VS = matrix_space(R,length(ev),1)
    for i in 1:length(Reductions)
        e = Basis[i][2] # pole order of basis element 
        N = N_m[e]
        verbose && println("e: $e")

        verbose && println(p*(e+N-1)-1)

        ff = factorial(ZZ(p*(e+N-1)-1)) 
        val_ff = valuation(ff,p)
        final_val = (n-1) - val_ff  
        ff_invertible = ff / p^val_ff

        inverse_ff = inv(R(ff_invertible))


        temp = VS()
        temp2 = convert_p_to_m([Reductions[i][1][1]],ev)
        for i in 1:length(ev)
            temp[i,1] = R(temp2[i])
        end
        temp = T * temp
        verbose && println("temp: $temp")
        for i in 1:length(temp)
#            println(temp[i])
            ele = inverse_ff * temp[i]
            if 0 ≤ final_val
                ele *= p^final_val
            else
                lifted = lift(ZZ,ele)
                ele = divexact(ele,p^(-final_val)) 
            end
            temp[i] = ele
        end 
        
        push!(frob_mat_temp, temp)
        #push!(denomArray, QQ(lift(ZZ,Reductions[i][3])))
        #push!(denomArray,lift(ZZ,Factorial(PrecisionRing(p*(Basis[i][2]+N-1)-1),PrecisionRing(1))/(p^(n-1))))
        #push!(FrobMatTemp,T*temp)
        #(p^(n-1)/Factorial(PrecisionRing(p*(Basis[i][2]+N-1)-1),PrecisionRing(1)))
    end
    frob_mat = hcat(frob_mat_temp...)
    #MS = matrix_space(ZZ,nrows(FrobMat),ncols(FrobMat))
    #FM = MS()
    #println(FrobMat)
    
    #=
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
    =#

    return frob_mat
end

"""
    LPolynomial(FM, q)

Given the Frobenius matrix, computes the corresponding L-polynomial det(1-tq^{-1}FM)

INPUT: 
* "FM" -- Frobenius matrix 
"""

function LPolynomial(FM,q,n)
    @assert size(FM, 1) == size(FM, 2) "FM is not a square matrix"

    P, T = polynomial_ring(parent(FM[1,1]), "T")
    f = charpoly(P, FM)

    k = degree(f)
    bound = binomial(k,Int(ceil(k/2)))*(q^(n/2))*ceil(k/2)
    R, t = polynomial_ring(ZZ,"t")
    result = R(0)
    for i in 0:k
        if abs(ZZ(coeff(f,i))) > bound
            result = result + (ZZ(coeff(f,i)) - characteristic(P))*t^(k-i)
        else
            result = result + (ZZ(coeff(f,i)))*t^(k-i)
        end
    end


    return reverse(f), result
end 


"""
    computeAll(f, verbose=false, givefrobmat=false)

Wrapper function that outputs the Frobenius Matrix

INPUTS: 
* "f" -- Oscar polynomial

KEYWORD ARGUMENTS:
verboselevel -- print statements at various levels of depth
givefrobmat -- should the funciton also output the appoximated frobenius matrix
algorithm -- the algorithm used for controlled reduction
termorder -- the term ordering that should be used in vector representations
>>>if you don't know what this is, ignore it.
vars_reversed -- reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.

"""
function zeta_function(f; verbose=false, givefrobmat=false, algorithm=:costachunks, termorder=:invlex, vars_reversed=true)
    p = Int64(characteristic(parent(f)))
    q = p
    n = nvars(parent(f)) - 1
    d = degree(f,1)
    PR = parent(f)
    R = coefficient_ring(parent(f))

    verbose && println("Working with a degree $d hypersurface in P^$n")

    basis = compute_monomial_bases(f, R, PR, termorder) # basis of cohomology 

    verbose && println("Basis of cohomology is $basis")

    k = sum([length(tmp) for tmp in basis]) # dimension of H^n

    verbose && println("There are $k basis elements in H^$n")

    r_m = relative_precision(k, p)
    N_m = series_precision(r_m, p, n) # series precision 
    M = algorithm_precision(r_m, N_m, p)

    verbose && println("We work modulo $p^$M, and compute up to the $N_m-th term of the Frobenius power series")
     
    precisionring, = residue_ring(ZZ, p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])

    T = computeT(f, basis, M, termorder)
    verbose && println("T matrix is $T")
    #S = SmallestSubsetSmooth.smallest_subset_s_smooth(fLift,n)
    S = collect(0:n)

    BasisTLift = []
    for i in basis
        temp = []
        for j in i
            push!(temp, liftCoefficients(precisionring,precisionringpoly,j))
        end
        push!(BasisTLift,temp)
    end

    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end

    verbose && println(Basis)

    fLift = liftCoefficients(precisionring, precisionringpoly, f)
    FBasis = applyFrobeniusToBasis(Basis,fLift, N_m, p, termorder)
    l = d * n - n + d - length(S)
    pseudo_inverse_mat_new = pseudo_inverse_controlled_lifted(f,S,l,M,termorder)
    MS = matrix_space(precisionring, nrows(pseudo_inverse_mat_new), ncols(pseudo_inverse_mat_new))
    pseudo_inverse_mat = MS()
    for i in 1:nrows(pseudo_inverse_mat_new)
        for j in 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i,j] = precisionring(ZZ(pseudo_inverse_mat_new[i,j]))
        end
    end
    #=
    pseudo_inverse_mat = zeros(Int, nrows(pseudo_inverse_mat_new),ncols(pseudo_inverse_mat_new))
    for i in 1:nrows(pseudo_inverse_mat_new)
        for j in 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i,j] = ZZ(pseudo_inverse_mat_new[i,j])
        end
    end
    printMat(pseudo_inverse_mat)
    =#
    #pseudoInverseMat = zeros(PrecisionRing, nrows(pseudoInverseMatTemp), ncols(pseudoInverseMatTemp))

    #PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    #fLift = liftCoefficients(ZZ,PRZZ,f)
    #controlledMatrixZZ = compute_controlled_matrix(fLift, d * n - n + d - length(S), S, ZZ, PRZZ)
    #pseudoInverseMatModP = matrix(ZZ, [lift(ZZ,x) for x in Array(pseudoInverseMatTemp)])
    #pseudo_inverse_mat_new = henselLift(p,M,controlledMatrixZZ, pseudoInverseMatModP)
    
    #for i in 1:nrows(pseudoInverseMat)
    #    for j in 1:ncols(pseudoInverseMat)
    #        pseudoInverseMat[i,j] = PrecisionRing(lift(ZZ, pseudoInverseMatTemp[i,j]))
    #    end
    #end
    Reductions = reducetransform_LA_descending(FBasis, N_m, S, fLift, pseudo_inverse_mat, p, termorder)
    verbose && println(Reductions)
    ev = gen_exp_vec(n+1,n*d-n-1,termorder)
    verbose && println(convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev))
    FM = compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, termorder)
    verbose && println(FM)

    if verbose
        println("The Frobenius matrix is $FM")
    end

    if givefrobmat
        (FM,LPolynomial(FM,q,n))
    else
        LPolynomial(FM,q,n)
    end
end

#end 

#=
include("DeRham.jl")
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
verbose = false
n = 2
d = 3
p = 7
R = GF(p,1)
PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
x0,x1,x2 = Vars
f = x1^2*x2 - x0^3 - x0*x2^2 - x2^3
S = [0,1,2]
Test = compute_all(f)
=#

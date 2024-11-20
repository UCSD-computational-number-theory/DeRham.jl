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
function compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, termorder, verbose)
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
    zeta_function(f; verbose=false, givefrobmat=false, algorithm=:costachunks, termorder=:invlex, vars_reversed=true)


Wrapper function that outputs the Frobenius Matrix

INPUTS: 
* "f" -- Oscar polynomial

KEYWORD ARGUMENTS:
TODO:verboselevel -- print statements at various levels of depth
RIGHTNOW: verbose -- print stuff
givefrobmat -- should the funciton also output the appoximated frobenius matrix
algorithm -- the algorithm used for controlled reduction
termorder -- the term ordering that should be used in vector representations
>>>if you don't know what this is, ignore it.
vars_reversed -- reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.

"""
function zeta_function(f; verbose=false, givefrobmat=false, algorithm=:costachunks, termorder=:invlex, vars_reversed=true, fastevaluation=false)
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

    hp = hodgepolygon(f; basis=basis)
    r_m = calculate_relative_precision(hp, n-1, p)
    #N_m = series_precision(r_m, p, n) # series precision 
    #M = algorithm_precision(r_m, N_m, p)
    N_m = series_precision(p,n,d,r_m)
    M = algorithm_precision(p,n,d,r_m,N_m)

    verbose && println("We work modulo $p^$M, and compute up to the $N_m-th term of the Frobenius power series")
     
    precisionring, = residue_ring(ZZ, p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])

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
    verbose && println("Basis of cohomology is $Basis")


    fLift = liftCoefficients(precisionring, precisionringpoly, f)
    FBasis = applyFrobeniusToBasis(Basis,fLift, N_m, p, termorder,vars_reversed,verbose=verbose)
    l = d * n - n + d - length(S)
    pseudo_inverse_mat_new = pseudo_inverse_controlled_lifted(f,S,l,M,termorder,vars_reversed)
    MS = matrix_space(precisionring, nrows(pseudo_inverse_mat_new), ncols(pseudo_inverse_mat_new))
    pseudo_inverse_mat = MS()

    T = computeT(f, basis, M, termorder, vars_reversed)
    verbose && println("T matrix is $T")

    for i in 1:nrows(pseudo_inverse_mat_new)
        for j in 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i,j] = ZZ(pseudo_inverse_mat_new[i,j])
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
    #TODO: check which algorithm we're using
    Reductions = reducetransform(FBasis, N_m, S, fLift, pseudo_inverse_mat, p, termorder,algorithm,vars_reversed,fastevaluation,verbose=verbose)
    if verbose
        for i in 1:length(Basis)
            basis_elt = Basis[i]
            after_reduction = Reductions[i]
            println("Basis element $basis_elt becomes $after_reduction after controlled reduction")
        end 
    end 
    ev = gen_exp_vec(n+1,n*d-n-1,termorder)
    verbose && println(convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev))
    FM = compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, termorder, verbose)
    display(FM)
    verbose && println("The Frobenius matrix is $FM")

    #reductions_verbose = convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev)

    #verbose && println("convert_p_to_m is $reductions_verbose")

    #FM = compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis)

   # verbose && println("The Frobenius matrix is $FM")

    if givefrobmat
        (FM,LPolynomial(FM,q,n))
    else
        LPolynomial(FM,q,n)
    end
end

"""
    hodgepolygon(basis::Array)

Calculates the hodge polygon of the cohomology module with 
griffiths-dwork basis basis

basis -- an array of "polynomials with pole" as descirbed in PolynomialWithPole.jl
"""
function hodgepolygon(basis::Array,n)
    #WRONG: n = highestpoleorder(basis)
    hodgenumbers = zeros(Int,n)
    for i in 0:n-1
        h = length(termsoforder(basis,n-i))
        hodgenumbers[i+1] = h
    end

    SlopesPolygon(hodgenumbers)
end

"""
    hodgepolygon(f; termorder=:invlex)

Calculates the hodge polygon of f

f - the polynomial to get the hodge polygon of
"""
function hodgepolygon(f; termorder=:invlex, basis=nothing)
    n = nvars(parent(f)) - 1
    PR = parent(f)
    R = coefficient_ring(parent(f))

    if basis == nothing
        basis = compute_monomial_bases(f, R, PR, termorder) # basis of cohomology 
    end
    
    Basis = []
    for i in 1:n
        for j in basis[i]
            push!(Basis,[j,i])
        end
    end

    hodgepolygon(Basis,n)
end


"""
    check_smoothness(f)

Using Oscar, check whether f defines a smooth hypersurface.

f is assumed to be homogeneous already. Otherwise Oscar
will throw an error.

note: the name issmooth / is_smooth is already taken by oscar,
and really we're just wrapping that method.
"""
function check_smoothness(f)
    p = characteristic(parent(f))
    nVars = length(gens(parent(f)))

    graded, _ = grade(parent(f))

    R, _ = quo(graded, ideal(graded, [graded(f)]))

    V = proj(R)

    is_smooth(V)
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
n = 4
p = 13
F = GF(p)
R, (x1,x2,x3,x4,x5) = polynomial_ring(F, ["x$i" for i in 0:n])

f = x1^3 + x2^3 + x3^3 + x4^3 + x5^3
@time DeRham.zeta_function(f)
=#

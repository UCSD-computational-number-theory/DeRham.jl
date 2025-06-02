#module StandardReduction

#using Oscar
#using BitIntegers
#using LinearAlgebra

#include("PrecisionEstimate.jl")
#include("CopiedFindMonomialBasis.jl")
#include("Utils.jl")
#include("CopiedFindMonomialBasis.jl")

"""
       monomial_change_basis(f, l, basis)
Returns the integer-coefficient change-of-basis matrix of Homog_l from the 
{cohomology relations,cohomology basis} basis to the standard monomial basis 

INPUTS: 
* "f" -- polynomial
* "l" -- integer, degree of homogeneous polynomials in question
* "basis" -- list, basis elements in the cohomology basis, assumed to be monomials of degree l 
* "termorder" -- the monomial ordering of vectors
"""
function monomial_change_basis(f, l, basis, params, cache)
       #println(basis)
       p = characteristic(parent(f))
       n = nvars(parent(f)) - 1
       S = [i for i in 0:n]
       d = total_degree(f)
       PR = parent(f)
       R = coefficient_ring(parent(f))
       PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
       f_lift = liftCoefficients(ZZ,PRZZ,f)

       # Lift coefficients of elements in basis to integers
       basis_lift = []
       for i in basis
              push!(basis_lift, liftCoefficients(ZZ,PRZZ,i,false))
       end

       exp_vec = gen_exp_vec(n+1, l, params.termorder)

       # matrix for the map (\mu_0, \dots, \mu_n) \mapsto \sum_{i\in n} \mu_i \partial_i f 
       change_basis_matrix = compute_controlled_matrix(f_lift,l,S,ZZ,PRZZ,params,cache)

       if length(basis) == 0
              return change_basis_matrix
       else 
              # column vectors corresponding to monomials in the basis of cohomology 
              basis_columns = transpose(convert_p_to_m(basis_lift,exp_vec))
              
              #change_basis_matrix_aug = hcat(change_basis_matrix,basis_columns); 
              change_basis_matrix_aug = hcat(basis_columns,change_basis_matrix); 
              return change_basis_matrix_aug
       end
end 

"""
       monomial_change_basis_inverse(f,l,basis)  
Computes the mod p inverse of the matrix output by monomial_change_basis(f,l,basis)

INPUTS: 
* "f" -- polynomial
* "l" -- integer, degree of homogeneous polynomials in question
* "basis" -- list, basis elements in the cohomology basis, assumed to be monomials of degree l 
* "termorder" -- the monomial ordering of vectors
"""
function monomial_change_basis_inverse(f,l,basis,params,cache)    
       PR = parent(f)
       R = coefficient_ring(PR)  
       A = monomial_change_basis(f,l,basis,params,cache)
       
       flag, B = Nemo.is_invertible_with_inverse(matrix(R,[R(x) for x in Array(A)]), side=:right)
       
       if flag
           return (A,B)
       else 
           throw(ArgumentError("matrix from f is not invertible"))
       end
   end
   
   """
       monomial_change_basis_inverse_lifted(f,l,basis,M)  
Computes the mod p^M Hensel-lifted inverse of the matrix output by monomial_change_basis(f,l,basis)

INPUTS: 
* "f" -- polynomial
* "l" -- integer, degree of homogeneous polynomials in question
* "basis" -- list, basis elements in the cohomology basis, assumed to be monomials of degree l 
* "M" -- integer, desired mod p^M precision
* "termorder" -- the monomial ordering of vectors
"""
function monomial_change_basis_inverse_lifted(f, l, basis, M, params, cache)
    PR = parent(f)
    R = coefficient_ring(PR)
    (A, Sol_fp) = monomial_change_basis_inverse(f,l,basis,params,cache)
    lift_to_int(s) = map(x -> lift(ZZ,x),s)

    Sol_mod_p_int = lift_to_int(Sol_fp)

    p = characteristic(parent(f))
    return A, henselLift(p, M, A, Sol_mod_p_int, params)
end


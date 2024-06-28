module StandardReduction

using Oscar
using BitIntegers
using LinearAlgebra

include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("Utils.jl")
include("CopiedFindMonomialBasis.jl")

"""
       monomial_change_basis(f, l, basis)
Returns the integer-coefficient change-of-basis matrix of Homog_l from the 
{cohomology relations,cohomology basis} basis to the standard monomial basis 

INPUTS: 
* "f" -- polynomial
* "l" -- integer, degree of homogeneous polynomials in question
* "basis" -- list, basis elements in the cohomology basis, assumed to be monomials of degree l 
"""
function monomial_change_basis(f, l, basis)
       p = characteristic(parent(f))
       n = nvars(parent(f)) - 1
       S = [i for i in 0:n]
       d = degree(f,1)
       PR = parent(f)
       R = coefficient_ring(parent(f))
       PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
       f_lift = Utils.liftCoefficients(ZZ,PRZZ,f)

       # Lift coefficients of elements in basis to integers
       basis_lift = []
       for i in basis
              push!(basis_lift, Utils.liftCoefficients(ZZ,PRZZ,i,false))
       end

       exp_vec = Utils.gen_exp_vec(n+1, l,:invlex)

       # matrix for the map (\mu_0, \dots, \mu_n) \mapsto \sum_{i\in n} \mu_i \partial_i f 
       change_basis_matrix = CopiedFindMonomialBasis.compute_controlled_matrix(f_lift,l,S,ZZ,PRZZ)
       
       # column vectors corresponding to monomials in the basis of cohomology 
       basis_columns = transpose(matrix(ZZ,Utils.convert_p_to_m(basis_lift,exp_vec)))
       
       #change_basis_matrix_aug = hcat(change_basis_matrix,basis_columns); 
       change_basis_matrix_aug = hcat(basis_columns,change_basis_matrix); 
       return change_basis_matrix_aug
end 

"""
       monomial_change_basis_inverse(f,l,basis)  
Computes the mod p inverse of the matrix output by monomial_change_basis(f,l,basis)

INPUTS: 
* "f" -- polynomial
* "l" -- integer, degree of homogeneous polynomials in question
* "basis" -- list, basis elements in the cohomology basis, assumed to be monomials of degree l 
"""
function monomial_change_basis_inverse(f,l,basis)    
       PR = parent(f)
       R = coefficient_ring(PR)  
       A = monomial_change_basis(f,l,basis)
       
       flag, B = is_invertible_with_inverse(matrix(R,[R(x) for x in Array(A)]), side=:left)
       
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
"""
   function monomial_change_basis_inverse_lifted(f, l, basis, M)
       PR = parent(f)
       R = coefficient_ring(PR)
       (A, Sol_fp) = monomial_change_basis_inverse(f,l,basis)
       lift_to_int64(s) = Int64.(map(x -> lift(ZZ,x),s))
   
       Sol_mod_p_int = lift_to_int64(Sol_fp)
   
       p = characteristic(parent(f))
       return A,Utils.henselLift(p,M,A,Sol_mod_p_int)
   end

end
#=
"""
       stdRed_step(f, h, standard_right_inverse, e, M)

Use Griffiths-Dwork relation to reduce the pole order of h Omega/f^m by 1 

INPUTS: 
* "f" -- polynomial, defining equation of hypersurface 
* "h" -- polynomial
* "standard_right_inverse" -- matrix, solution to linear algebra problem
* "e" -- integer, pole order of differential h Omega/f^e
* "M" -- integer, mod p^M precision
"""
function stdRed_step(f, h, standard_right_inverse, e, target_exp_vec, domain_exp_vec, M)
       # Number of variables, n
       n = nvars(parent(poly)) - 1
       S = [i for i in 0:n]
       d = degree(f,1)
       l = d*n - n - 1
       PR = parent(f)
       R = coefficient_ring(parent(f))
       PrecisionRing, = residue_ring(ZZ,p^M)
       PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])

       h_vec = Utils.convert_p_to_m(h, mono_vec)
       mu_vec = standard_right_inverse * transpose(h_vec)
       mu = []
       for i in 1:(n+1)
              push!(mu, Utils.convert_m_to_p(transpose(mu_vec[Int((i-1)*(length(gJS)/(n+1))+1):Int(i*(length(gJS)/(n+1)))]),domain_exp_vec,PrecisionRing,PrecisionRingPoly)[1])
       end
       muPartials = [derivative(mu[i], i) for i in 1:(n+1)]
       return sum(muPartils)
end 

end
=#
       #=
       polyMat = Utils.convert_p_to_m([poly],Utils.gen_exp_vec(n+1,d*f_tilde_exp - n - 1))

       # Groebner basis G_i
       pseudoInverseMat = CopiedFindMonomialBasis.pseudo_inverse_controlled(f,[],R,PR)[2]

       polycTemp = pseudoInverseMat*transpose(polyMat)
       polyc = []
       for i in 1:(n+1)
              push!(polyc, Utils.convert_m_to_p(transpose(polycTemp[Int((i-1)*(length(polycTemp)/(n+1))+1):Int(i*(length(polycTemp)/(n+1)))]),Utils.gen_exp_vec(n+1,n*d-n),R,PR)[1])
       end
       partials = [ derivative(polyc[i], i) for i in 1:(n+1) ]


       # partial derivatives of G_i, and \sum_i{\partial_i{G_i}}
       result = sum(partials)

       return (result, f_tilde_exp - 1, denom*(f_tilde_exp-1))
       
       fPartials = [ derivative(f, i) for i in 1:(n+1) ]
       fPartials = reverse(fPartials)
       polyc, t = reduce_with_quotients(poly,fPartials)
       polyPartials = [ derivative(polyc[i], i) for i in 1:(n+1) ]
       result = sum(polyPartials)
       return (result, f_tilde_exp - 1, denom*(f_tilde_exp-1))
       =#

# reduction: exponent downto #variables
#=
function stdandardReduction(f, poly, f_tilde_exp, denom)
       # Number of variables, n
       num_vars = nvars(parent(poly)) - 1

       while f_tilde_exp > num_vars
          (poly, f_tilde_exp, denom) = stdRed_step(f, poly, f_tilde_exp, denom)
       end

       return (poly, f_tilde_exp, denom)
end
=#
#
# pseudo_inverse_classical(f, R, PR), where R is the field and PR associated poly ring
#CopiedFindMonomialBasis.pseudo_inverse_classical(f, R, PR)




#=
# Examples
p = 7
R = GF(p)
PR, vars = polynomial_ring(R, ["x", "y", "z"])
x, y, z = vars
f2 = x^2 * z - x^3 - x * z^2 - z^3

p = 7
R = GF(p)
# R = residue_ring(ZZ, 7)
PR, vars = polynomial_ring(R, ["w", "x", "y", "z"])
w, x, y, z = vars
f3 = x^4 + y^4 + z^4 + w^4
=#

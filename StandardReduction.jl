module StandardReduction

using Oscar
using BitIntegers
using LinearAlgebra

include("PrecisionEstimate.jl")
include("FindMonomialBasis.jl")
include("CopiedFindMonomialBasis.jl")
include("AutomatedScript.jl")
include("CopiedFindMonomialBasis.jl")

# One step reduction (formula)
# 
# @param poly: polynomial f
# @param f_tilde_exp: exponent of \tilde{f} in the denominator
# @param denom: the coefficient in the denominator (p-denominator)
function stdRed_step(poly, f_tilde_exp, denom)
       # Number of variables, n
       n = nvars(parent(poly)) - 1

       # Groebner basis G_i
       groeb_basis = CopiedFindMonomialBasis.compute_monomial_bases(poly, base_ring(parent(poly)), parent(poly))

       # partial derivatives of G_i, and \sum_i{\partial_i{G_i}}
       partials = [ derivative(groeb_basis, i) for i in 1:n+1 ]
       result = sum(partials)

       return (result, f_tilde_exp - 1, denom / (f_tilde_exp - 1))
end

# reduction: exponent downto #variables
function stdandardReduction(poly, f_tilde_exp, denom)
       # Number of variables, n
       num_vars = nvars(parent(poly)) - 1

       while f_tilde_exp > num_vars
          (poly, f_tilde_exp, denom) = stdRed_step(poly, f_tilde_exp, denom)
       end

       return (poly, f_tilde_exp, denom)
end
#
# psuedo_inverse_classical(f, R, PR), where R is the field and PR associated poly ring
#CopiedFindMonomialBasis.psuedo_inverse_classical(f, R, PR)

end


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
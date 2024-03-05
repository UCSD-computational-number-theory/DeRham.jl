module StandardReduction

using Oscar
using BitIntegers
using LinearAlgebra

include("PrecisionEstimate.jl")
include("FindMonomialBasis.jl")
include("AutomatedScript.jl")
include("CopiedFindMonomialBasis.jl")

# One step reduction (formula)
# 
# @poly: polynomial f
# @f_tilde_exp: exponent of \tilde{f} in the denominator
# @denom: the coefficient produced with the formula, in the denominator (p-denominator)
function stdRed_step(poly, f_tilde_exp, denom)
       # Number of variables, n
       num_vars = nvars(parent(poly)) - 1

       # Groebner basis G_i
       groeb_basis = FindMonomialBasis.find_monomial_basis(poly, n, false)

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

# psuedo_inverse_classical(f, R, PR), where R is the field and PR associated poly ring
CopiedFindMonomialBasis.psuedo_inverse_classical(f, R, PR)

end

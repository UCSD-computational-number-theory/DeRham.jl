module StandardReduction

using Oscar
using BitIntegers
using LinearAlgebra

include("PrecisionEstimate.jl")
include("FindMonomialBasis.jl")
include("CopiedFindMonomialBasis.jl")
include("Utils.jl")
include("CopiedFindMonomialBasis.jl")

"""
       stdRed_step(f, poly, f_tilde_exp, denom)

Standard reduction step 

INPUTS: 
* "f" -- polynomial, defining equation of hypersurface 
* "poly" -- ? 
* "f_tilde_exp" -- vector?, exponents of f? 
* "denom" -- integer?, coefficient in the denominator of what? 
"""
function stdRed_step(f, poly, f_tilde_exp, denom)
       # Number of variables, n
       n = nvars(parent(poly)) - 1
       d = degree(f,1)
       PR = parent(f)
       R = coefficient_ring(parent(f))
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
       =#
       fPartials = [ derivative(f, i) for i in 1:(n+1) ]
       polyc, t = reduce_with_quotients(poly,fPartials)
       polyPartials = [ derivative(polyc[i], i) for i in 1:(n+1) ]
       result = sum(polyPartials)
       return (result, f_tilde_exp - 1, denom*(f_tilde_exp-1))
end

# reduction: exponent downto #variables
function stdandardReduction(f, poly, f_tilde_exp, denom)
       # Number of variables, n
       num_vars = nvars(parent(poly)) - 1

       while f_tilde_exp > num_vars
          (poly, f_tilde_exp, denom) = stdRed_step(f, poly, f_tilde_exp, denom)
       end

       return (poly, f_tilde_exp, denom)
end
#
# pseudo_inverse_classical(f, R, PR), where R is the field and PR associated poly ring
#CopiedFindMonomialBasis.pseudo_inverse_classical(f, R, PR)

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

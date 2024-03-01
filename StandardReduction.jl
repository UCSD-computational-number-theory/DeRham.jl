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
# @F_poly: polynomial f
# @f_tilde_exp: exponent of \tilde{f} in the denominator
# @p_denom: the coefficient produced with the formula, in the denominator
function stdRed_step(F_poly, f_tilde_exp, p_denom)
       # Number of variables, n
       num_vars = nvars(parent(F_poly)) - 1

       # Groebner basis G_i
       groeb_basis = FindMonomialBasis.find_monomial_basis(F_poly, n, false)

       # partial derivatives of G_i, and \sum_i{\partial_i{G_i}}
       partials = [ derivative(groeb_basis, i) for i in 1:n+1 ]
       result = sum(partials)

       return (result, f_tilde_exp - 1, p_denom / (f_tilde_exp - 1))
end

# reduction: exponent downto #variables
function stdandardReduction(F_poly, f_tilde_exp, p_denom)
       # Number of variables, n
       num_vars = nvars(parent(F_poly)) - 1

       while f_tilde_exp > num_vars
          (F_poly, f_tilde_exp, p_denom) = stdRed_step(F_poly, f_tilde_exp, p_denom)
       end

       return (F_poly, f_tilde_exp, p_denom)
end

# psuedo_inverse_classical(f, R, PR), where R is the field and PR associated poly ring

#=
# (From Thomas' code)
# Finds a monomial basis for polynomial de Rham cohomology with Z(f)
# for each case of `m = 1 ... n`.
function find_monomial_basis(f, m)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    if m * d - n - 1 <= 0 || m * d - n - d < 0
        return []
    end

    rows = binomial(m*d-1, n)
    section = binomial(m*d - d, n)
    cols = section * (n + 1)
    U = matrix_space(R, rows, cols)
    M = U()

    mon = compute_monomials(n + 1, m*d - n - d)
    partials = [ derivative(f, i) for i in 1:n+1 ]

    # Computes (md - 1 \choose n) x (n + 1)(md - d \choose n) matrix corresponding
    # to linear transformation
    #   F:  (Homog_{md-n-d})^{n+1}  ->  Homog_{md-n-1}
    #       (m_0, ..., m_n)        |->  \sum_{i=0}^n m_i \partial_i f
    for i in 1:n+1
        for monomial in 1:length(mon)
            current_term = polynomial_to_lex_vector(mon[monomial] * partials[i], n)
            M[:, section * (i-1) + monomial] = current_term
        end
    end

    row_monomials = compute_monomials(n + 1, m*d - n - 1)
    non_pivot_rows = find_non_pivot_rows(M)
    return map((i) -> row_monomials[i], non_pivot_rows)
end
=#

end

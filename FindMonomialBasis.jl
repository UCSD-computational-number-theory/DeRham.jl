module FindMonomialBasis
using Oscar


# Coefficient Ring
# p = 7
R = ZZ # GF(p)

# Polynomial Setup

# n = 2 # number of variables - 1
# d = 3 # degree of homogenous polynomial
# S, vars = PolynomialRing(R, ["x$i" for i in 0:n])
# x, y,z = vars
# f = -x^3 - x*z^2 + y^2*z - z^3

n = 2
d = 5
S, vars = PolynomialRing(R, ["x$i" for i in 0:n])
x, y, z = vars
f = x^5 + y^5 + z^5 + x * y^3 * z

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

# Find all monomial bases
function find_monomial_bases(f)
    n = nvars(parent(f)) - 1
    res = []
    for m in 1:n
        push!(res, (m, find_monomial_basis(f, m)))
    end
    return res
end

# Given a matrix M, find rows number corresponding to non-pivot
# rows.
function find_non_pivot_rows(M)
    res = []
    N = rref(M)[2]
    for i in 1:size(N, 1)
        if all(N[i, :] .== 0)
            push!(res, i)
        end
    end
    return res
end

# Computes all monomials of degree `d` in `n` variables.
function compute_monomials(n, d)
    result = []
    
    function backtrack(start, current)
        if length(current) == d
            push!(result, prod(S(var) for var in current))
            return
        end

        for i in start:n
            backtrack(i, [current..., vars[i]])
        end
    end

    backtrack(1, [])

    return result
end

# Given a homogenous polynomial `f` of degree `d`` in `n` variables, computes
# the coefficient vector
function polynomial_to_lex_vector(f, n)
    d = total_degree(f)

    mon = compute_monomials(n + 1, d)
    res = fill(0, length(mon))
    for i in 1:length(mon)
        res[i] = coeff(f, mon[i])
    end
    return res
end

# Given a lexicographically ordered vector of monomial coefficients, returns
# the associated polynomial
function lex_vec_to_polynomial(vect, n, d)
    res = S()
    mon = compute_monomials(n + 1, d)
    for i in 1:length(vect)
        res += S(vect[i]) * mon[i]
    end
    return res
end
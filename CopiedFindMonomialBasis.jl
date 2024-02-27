module CopiedFindMonomialBasis
using Oscar

macro assert(ex)
    return :( $ex ? nothing : throw(AssertionError($(string(ex)))) )
end

# Polynomial Setup

# n = 2 # number of variables - 1
# d = 3 # degree of homogenous polynomial
# S, vars = PolynomialRing(R, ["x$i" for i in 0:n])
# x, y,z = vars
# f = -x^3 - x*z^2 + y^2*z - z^3

# Examples:
# p = 7,11,13
# y^2 - x^3 - x
# y^2 - x^3 - x - 1
# x^4 + y^4 + z^4 + w^4
# x^4 + y^4 + z^4 + w^4 + 2xyzw
# x^5 + y^5 + z^5 + xy^3z

p = 41
R = GF(p)
n = 2
d = 5
PR, vars = PolynomialRing(R, ["x$i" for i in 0:n])
x, y, z = vars
f = x^5 + y^5 + z^5 + x * y^3 * z

function pivot_columns(M)
    rank, M = rref(M)
    res = fill(0, rank)
    ncols = size(M, 2)
    j = 1
    for i in 1:rank
        while j <= ncols && M[i, j] == 0
            j += 1
        end
        res[i] = j
        j+=1
    end
    return res
end

# Computes the basis
function main_find_basis(f, p, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    basis, _ = find_monomial_basis(f, 2, false, p, R, PR)

    return basis
end

function psuedo_inverse_classical(f, p, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    _, M = find_monomial_basis(f, 2, false, p, R, PR)

    flag, B = is_invertible_with_inverse(M, side=:right)

    if flag
        return Array(B)
    end
end

function psuedo_inverse_controlled(f, p, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    _, M = find_monomial_basis(f, 2, true, p, R, PR)

    flag, B = is_invertible_with_inverse(M, side=:right)

    if flag
        return Array(B)
    end
end

# Find all monomial bases
# res = find_monomial_bases(f)
# basis = res[m][2]
# uncontrolled_matrix = res[m][3]
# controlled_matrix = res[m][4]

# (m, basis, matrix, matrix)
function find_monomial_bases(f)
    n = nvars(parent(f)) - 1
    res = []
    for m in 1:n
        if m*d - n - d <= 0
            push!(res, (m, compute_monomials(m*d - n - 1, d), nothing, nothing))
        else
            uncontrolled = find_monomial_basis(f, m, false)
            controlled = find_monomial_basis(f, m, true)
            push!(res, (m, uncontrolled[1], uncontrolled[2], controlled[2]))
        end
    end
    return res
end

# Computes all monomials of degree `d` in `n` variables.
function compute_monomials(n, d)
    result = []
    
    function backtrack(start, current)
        if length(current) == d
            push!(result, prod(PR(var) for var in current))
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
    res = fill(R(0), length(mon))
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
        res += PR(vect[i]) * mon[i]
    end
    return res
end

end
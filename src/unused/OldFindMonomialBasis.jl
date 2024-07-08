module FindMonomialBasis
using Oscar


# M_l g = (mu0, ... mn)


macro assert(ex)
    return :( $ex ? nothing : throw(AssertionError($(string(ex)))) )
end

# Coefficient Ring
# p = 7
p = 41
R = GF(p) # GF(p)

# Polynomial Setup

# n = 2 # number of variables - 1
# d = 3 # degree of homogenous polynomial
# S, vars = polynomial_ring(R, ["x$i" for i in 0:n])
# x, y,z = vars
# f = -x^3 - x*z^2 + y^2*z - z^3


# 1. take nullspace of A^T
# 2. Add to original matrix these as cols
# 3. Now we have a square matrix, so the columns of this forms a basis for
#   for the monomials we were given
# 3. take inverse of this mat
# 4. Now we multiple the inverse A^-1 by the elements of the nullspace we found vi
#   to get A vi the associated representation in the matrix

n = 2
d = 5
S, vars = polynomial_ring(R, ["x$i" for i in 0:n])
x, y, z = vars
f = x^5 + y^5 + z^5 + x * y^3 * z

# Finds a monomial basis for polynomial de Rham cohomology with Z(f)
# for each case of 1 < m < n.

function find_monomial_basis(f, m, controlled)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    if controlled
        rows = binomial(m*d - 1, n)
        sec = binomial(m*d - d - 1, n)
    else
        rows = binomial(m*d-1, n)
        sec = binomial(m*d - d, n)
    end
    cols = sec * (n + 1)
    U = matrix_space(R, rows, cols)
    M = U()

    # Uncontrolled: rows = binomial(m*d - n - d, n)
    # if m*d - n - d - 1 > 0:
    if controlled
        mon = compute_monomials(n + 1, m*d - n - d - 1)
    else
        mon = compute_monomials(n + 1, m*d - n - d)
    end
    partials = [ derivative(f, i) for i in 1:n+1 ]

    # Computes (md - 1 \choose n) x (n + 1)(md - d \choose n) matrix corresponding
    # to linear transformation
    #   F:  (Homog_{md-n-d})^{n+1}  ->  Homog_{md-n-1}
    #       (m_0, ..., m_n)        |->  \sum_{i=0}^n m_i [x_i] \partial_i f
    for i in 1:n+1
        for monomial in 1:length(mon)
            if controlled
                current_term = polynomial_to_lex_vector(vars[i] * mon[monomial] * partials[i], n)
            else
                current_term = polynomial_to_lex_vector(mon[monomial] * partials[i], n)
            end
            # Uncontrolled: current_term = polynomial_to_lex_vector(mon[monomial] * partials[i], n)
            M[:, sec * (i-1) + monomial] = current_term
        end
    end

    nrows = size(M, 1)
    ncols = size(M, 2)
    pivot_rows = pivot_columns(transpose(M))
    non_pivot_rows = setdiff([1:nrows;], pivot_rows)

    # Uncontrolled: row_monomials = compute_monomials(n + 1, m*d - n - 1)
    row_monomials = compute_monomials(n + 1, m*d - n - 1)

    # non_pivot_rows = find_non_pivot_rows(M)
    return map((i) -> row_monomials[i], non_pivot_rows), M
end

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
        res += S(vect[i]) * mon[i]
    end
    return res
end

end
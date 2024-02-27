module FindMonomialBasis
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
# x^5 + y^5 + z^5 + x

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

# Given $f, m$, over the ring $\texttt{PR} = R[x_0, \dots, x_n]$, computes the matrix for the map
# $$(\mu_0, \dots, \mu_n) \mapsto \sum_{i=0}^n \mu_i \partial_i f$$
# where each $\mu_i$ are monomials of degree $dm - n - 1 - (d-1)$.
#
# Usage: "Linear Algebra Problem", l = d*m - n - 1. (See Section 1.2.2)
#
function compute_classical_matrix(f, l, m, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    @assert(0 <= m && m <= n)

    section = binomial(n + l - (d-1), n)
    domain_mons = compute_monomials(n+1, l - (d - 1), PR)

    if length(domain_mons) <= 0
        return []
    end
    
    U = matrix_space(R, binomial(n + l, n), (n+1) * section)
    M = U()

    partials = [ derivative(f, i) for i in 1:n+1 ]

    for i in 1:n+1
        for monomial in eachindex(domain_mons)
            M[:, section * (i-1) + monomial] = polynomial_to_vector(domain_mons[monomial] * partials[i], n+1, R, PR, order=:lex)
        end
    end
    
    return M
end

# Given $f, m, l$, over the ring $\texttt{PR} = R[x_0, \dots, x_n]$, and the set $S\subseteq \{0, ..., n\}$, computes
# the matrix for the map
# $$(\mu_0, \dots, \mu_n) \mapsto \sum_{i\in S} \mu_i \partial_i f + \sum_{i\not\in S} \mu_i x_i \partial_i f$$
# where $\mu_i$ for $0\leq i < |S|$ are of degree $l - (d - 1)$ and $\mu_i$ for
# $|S| \leq i \leq n$ are of degree $l - d$.
#
# Usage: If $f$ is nondegenerate, you can use this for cases $l=0,d,\dots, dn$. (See Prop 1.18)
# Usage: Use this for case $l = dn - n + d - |S|$. (See Prop 1.15)
#
function compute_controlled_matrix(f, l, S, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    len_S = length(S)
    
    @assert(len_S >= 0 && len_S <= n+1)

    in_set_mons = compute_monomials(n+1, l - (d - 1), PR)
    not_in_set_mons = compute_monomials(n+1, l - d, PR)

    in_set_section = binomial(n + l - (d-1), n)
    not_in_set_section =  binomial(n + l - d, n)
    cols = len_S * in_set_section + (n + 1 - len_S) * not_in_set_section

    if len_S > 0
        @assert(length(in_set_mons) > 0)
    end
    if len_S < n+1
        @assert(length(not_in_set_mons) > 0)
    end
    
    U = matrix_space(R, binomial(n + l, n), cols)
    M = U()

    partials = [ derivative(f, i) for i in 1:n+1 ]

    for i in 1:len_S
        for monomial in eachindex(in_set_mons)
            M[:, in_set_section * (i-1) + monomial] = polynomial_to_vector(in_set_mons[monomial] * partials[i], n+1, R, PR, order=:lex)
        end
    end

    for i in (len_S+1):n+1
        for monomial in eachindex(not_in_set_mons)
            M[:, not_in_set_section * (i-1) + monomial] = polynomial_to_vector(not_in_set_mons[monomial] * vars[i] * partials[i], n+1, R, PR, order=:lex)
        end
    end
    return M
end

# Computes the monomial basis of $H_{dR}^n(U_{\QQ_p})$. In particular, we find the monomials
# of degree $l = dm - n - 1$ in $F_p[x_0, \dots, x_n]$ that project onto a basis of the cokernel
# of the map computed in `compute_classical_mat()`.
function compute_monomial_basis(f, m, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    row_monomials = compute_monomials(n + 1, m*d - n - 1, PR)

    M = compute_classical_mat(f, d*m - n - 1, m, R, PR)
    if isempty(M)
        return row_monomials
    end

    rows = size(M)[1]

    pivot_rows = pivot_columns(transpose(M))
    non_pivot_rows = setdiff([1:rows;], pivot_rows)
    return map((i) -> row_monomials[i], non_pivot_rows)
end

# Computes the the monomial bases for different `m`. That is,
# `compute_monomial_bases(f, R, PR)[m]` will give the `m`-th case.
function compute_monomial_bases(f, R, PR)
    n = nvars(parent(f)) - 1

    res = []

    for m in 1:n
        push!(res, compute_monomial_basis(f, m, R, PR))
    end
    return res
end

function psuedo_inverse_classical(f, m, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    
    M = compute_classical_matrix(f, d*m - n - 1, m, R, PR)

    flag, B = is_invertible_with_inverse(M, side=:left)

    if flag
        return Array(B)
    end
end

function psuedo_inverse_controlled(f, S, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    
    len_S = length(S)

    M = compute_controlled_matrix(f, d * n - n + d - len_S, S, R, PR)

    flag, B = is_invertible_with_inverse(M, side=:right)

    if flag
        return Array(B)
    end
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
function compute_monomials(n, d, PR)
    if d <= 0
        return []
    end

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
function polynomial_to_vector(f, n, R, PR; order=:lex)

    if order == :lex
        d = total_degree(f)

        mon = compute_monomials(n, d, PR)
        res = fill(R(0), length(mon))
        for i in eachindex(mon)
            res[i] = coeff(f, mon[i])
        end
        return res
    else
        throw(ArgumentError("Invalid option '$order'"))
    end
end

# Given a lexicographically ordered vector of monomial coefficients, returns
# the associated polynomial
function lex_vec_to_polynomial(vect, n, d, PR)
    res = PR()
    mon = compute_monomials(n + 1, d, PR)
    for i in 1:length(vect)
        res += PR(vect[i]) * mon[i]
    end
    return res
end

end
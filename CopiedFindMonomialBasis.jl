module CopiedFindMonomialBasis

using Oscar

include("Utils.jl")

# Polynomial Setup

# n = 2 # number of variables - 1
# d = 3 # degree of homogenous polynomial
# S, vars = polynomial_ring(R, ["x$i" for i in 0:n])
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
PR, vars = polynomial_ring(R, ["x$i" for i in 0:n])
x, y, z = vars
f = x^5 + y^5 + z^5 + x * y^3 * z

# Given $f, m$, over the ring $\texttt{PR} = R[x_0, \dots, x_n]$, computes the matrix for the map
# $$(\mu_0, \dots, \mu_n) \mapsto \sum_{i=0}^n \mu_i \partial_i f$$
# where each $\mu_i$ are monomials of degree $dm - n - 1 - (d-1)$.
#
function compute_basis_matrix(f, l, m, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    vars = gens(PR)

    @assert(0 <= m && m <= n)

    section = binomial(n + l - (d-1), n)
    domain_mons = Utils.compute_monomials(n+1, l - (d - 1), PR,:invlex)

    if length(domain_mons) <= 0
        return []
    end
    
    U = matrix_space(R, binomial(n + l, n), (n+1) * section)
    M = U()

    partials = [ derivative(f, i) for i in 1:n+1 ]

    for i in 1:n+1
        for monomial in eachindex(domain_mons)
            M[:, section * (i-1) + monomial] = Utils.polynomial_to_vector(domain_mons[monomial] * partials[i], n+1, R, PR, :invlex)
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
function compute_controlled_matrix(f, l, S, R, PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    vars = gens(PR)

    len_S = length(S)
    
    @assert(len_S >= 0 && len_S <= n+1)

    in_set_mons = Utils.compute_monomials(n+1, l - (d - 1), PR,:invlex)
    not_in_set_mons = Utils.compute_monomials(n+1, l - d, PR,:invlex)

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
    partials = reverse(partials)

    for i in 1:len_S
        for monomial in eachindex(in_set_mons)
            M[:, in_set_section * (i-1) + monomial] = Utils.polynomial_to_vector(in_set_mons[monomial] * partials[i], n+1, R, PR, :invlex)
        end
    end

    for i in (len_S+1):n+1
        for monomial in eachindex(not_in_set_mons)
            M[:, not_in_set_section * (i-1) + monomial] = Utils.polynomial_to_vector(not_in_set_mons[monomial] * vars[i] * partials[i], n+1, R, PR, :invlex)
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
    vars = gens(PR)

    ev = Utils.gen_exp_vec(n + 1, m*d - n - 1)
    row_monomials = Utils.gen_mon(ev,R,PR)

    M = compute_basis_matrix(f, d*m - n - 1, m, R, PR)
    if isempty(M)
        return row_monomials
    end

    rows = size(M)[1]

    pivot_rows = Utils.pivot_columns(transpose(M))
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


# Computes the pseudo_inverse for the controlled case.
"""
Solves the linear algebra problem in section 1.5.2 of Costa's thesis, 
top of page 23. 
In other words, finds a pseudo-inverse to the linear map described
there.
The map is constructed as a matrix from the polynomial f and the set S.

f - the polynomial defining the hypersurface
S - the set in [0..n] to be used for the linear algebra problem

I think these are correct: (TODO)
R - coefficient_ring(parent(f))
PR- paren(f)
"""
function pseudo_inverse_controlled(f, S, l, R, PR)
    n = nvars(parent(f)) - 1
    
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    fLift = Utils.liftCoefficients(ZZ,PRZZ,f)
    
    U = compute_controlled_matrix(fLift, l, S, ZZ, PRZZ)
    
    flag, B = is_invertible_with_inverse(matrix(R,[R(x) for x in Array(U)]), side=:right)
    
    if flag
        return (U,B)
    else 
        throw(ArgumentError("matrix from f is not right invertible"))
    end
end

function pseudo_inverse_controlled_lifted(f,S,l,M)
    PR = parent(f)
    R = coefficient_ring(PR)
    (U, Sol_fp) = pseudo_inverse_controlled(f,S,l,R,PR)
    lift_to_int64(s) = Int64.(map(x -> lift(ZZ,x),s))

    Sol_mod_p_int = lift_to_int64(Sol_fp)

    println("Solution mod p: $Sol_fp")

    p = characteristic(parent(f))
    return Utils.henselLift(p,M,U,Sol_mod_p_int)
end

# Computes the pseudo_inverse for the classical case.
#TODO: update this to reflex changes to pseudo_inverse_controlled
#it's used in standard reduction, I'll plan to take care of it then
function pseudo_inverse_classical(f, R, PR)
    return pseudo_inverse_controlled(f, [i for i in 1:n+1], R, PR)
end

function pseudo_inverse_classicalm(f, m, R, PR)
    return pseudo_inverse_controlled(f, [i for i in 1:n+1], R, PR)
end

end

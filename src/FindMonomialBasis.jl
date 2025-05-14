#module FindMonomialBasis

#using Oscar


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

#p = 41
#R = GF(p)
#n = 2
#d = 5
#PR, vars = polynomial_ring(R, ["x$i" for i in 0:n])
#x, y, z = vars
#f = x^5 + y^5 + z^5 + x * y^3 * z

# Given $f, m$, over the ring $\texttt{PR} = R[x_0, \dots, x_n]$, computes the matrix for the map
# $$(\mu_0, \dots, \mu_n) \mapsto \sum_{i=0}^n \mu_i \partial_i f$$
# where each $\mu_i$ are monomials of degree $dm - n - 1 - (d-1)$.
#
function compute_basis_matrix(f, l, m, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    vars = gens(PR)
    R = base_ring(PR)

    @assert(0 <= m && m <= n)

    section = binomial(n + l - (d-1), n)
    domain_mons = compute_monomials(n+1, l - (d - 1), PR, params.termorder, cache, vars_reversed=cache.vars_reversed)

    if length(domain_mons) <= 0
        return []
    end
    
    M = zero_matrix(R, binomial(n + l, n), (n+1) * section)

    partials = [ derivative(f, i) for i in 1:n+1 ]

    #TODO: if the user makes a mistake and puts in a surface that
    #has two few variables, the error happens here.
    #We should probably double check for this in ZetaFunction.jl
    for i in 1:n+1
        for monomial in eachindex(domain_mons)
            M[:, section * (i-1) + monomial] = polynomial_to_vector(domain_mons[monomial] * partials[i], n+1, params.termorder, cache, vars_reversed=cache.vars_reversed)
        end
    end
    
    return M
end


"""
    compute_controlled_matrix(f, l, S, R, PR, params)

Given f, m, l, over the ring PR = R[x_0, ..., x_n], and the set S subseteq {0, ..., n}, computes
the matrix for the map
(mu_0, ..., mu_n) |--> sum_{i in S} mu_i partial_i f + sum_{i notin S} mu_i x_i partial_i f
 where mu_i for 0<= i < |S| are of degree l - (d - 1) and mu_i for
 |S| <= i <= n are of degree l - d.

"""
function compute_controlled_matrix(f, l, S, R, PR, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    vars = gens(PR)

    len_S = length(S)
    notS = setdiff(collect(0:n),S)
    # notS = [2]
    # Stilda = [1,1,0]
    len_notS = length(notS)
    
    @assert(len_S >= 0 && len_S <= n+1)

    Stilda = zeros(Int, n+1)
    for i in S
        Stilda[i+1] = 1
    end

    # The code that was formerly here looked as follows:
    #
#    in_S_mons_vec = gen_exp_vec(n+1, l-(d-1), params.termorder)
#    not_in_S_mons_vec = gen_exp_vec(n+1, l-d, params.termorder)
    #
    # Thus, I figured that when we get them from the cache, which 
    # is calculated using vars_reversed=true in gen_exp_vec,
    # I figured we'd need to reverse here:
    
    in_S_mons_vec = cache[l-(d-1)]
    not_in_S_mons_vec = cache[l-d]


    if cache.vars_reversed
        reverse!.(in_S_mons_vec)
        reverse!.(not_in_S_mons_vec)
    end

    in_S_mons = gen_mon(in_S_mons_vec, PR)
    not_in_S_mons = gen_mon(not_in_S_mons_vec, PR)

    if cache.vars_reversed
        reverse!.(in_S_mons_vec)
        reverse!.(not_in_S_mons_vec)
    end

    in_set_section = binomial(n + l - (d-1), n)
    not_in_set_section =  binomial(n + l - d, n)
    cols = len_S * in_set_section + (n + 1 - len_S) * not_in_set_section

    if len_S > 0
        @assert(length(in_S_mons_vec) > 0)
    end
    if len_S < n+1
        @assert(length(not_in_S_mons_vec) > 0)
    end
    
    U = matrix_space(R, binomial(n + l, n), cols)
    M = U()

    partials = zeros(parent(f),n+1)
    for i in 1:n+1
        partials[i] = derivative(f, i)
    end
    #partials = [ derivative(f, i) for i in 1:n+1 ] # ∂_0, ∂_1, ∂_2
    
    if params.vars_reversed == true
        partials = reverse(partials)  # one needs to be quite careful with the ordering of partials 
        Stilda = reverse(Stilda)
        vars = reverse(vars)        
    end
    #println("partials = $partials")
    #println("Stilda = $Stilda")
    #println("vars = $vars")
    #println("in_S_mons = $in_S_mons")
    #println("in_S_mons_vec = $in_S_mons_vec")

    col_idx = 1
    for i in 1:(n+1)
        if Stilda[i] == 1
            for monomial in eachindex(in_S_mons)
                #M[:, col_idx] = polynomial_to_vector(in_S_mons[monomial] * partials[i], n+1, params.termorder)
                #println(in_S_mons[monomial] * partials[i])
                M[:, col_idx] = polynomial_to_vector(in_S_mons[monomial] * partials[i], n+1, params.termorder, cache, vars_reversed=cache.vars_reversed)
                col_idx = col_idx + 1
            end
        else
            for monomial in eachindex(not_in_S_mons)
                #M[:, col_idx] = polynomial_to_vector(not_in_S_mons[monomial] * vars[i] * partials[i], n+1, params.termorder)
                M[:, col_idx] = polynomial_to_vector(not_in_S_mons[monomial] * vars[i] * partials[i], n+1, params.termorder, cache, vars_reversed=cache.vars_reversed)
                col_idx = col_idx + 1
            end
        end 
    end 

    return M
end

# Computes the monomial basis of $H_{dR}^n(U_{\QQ_p})$. In particular, we find the monomials
# of degree $l = dm - n - 1$ in $F_p[x_0, \dots, x_n]$ that project onto a basis of the cokernel
# of the map computed in `compute_classical_mat()`.
function compute_monomial_basis(f, m, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    vars = gens(PR)

    ev = cache[m*d - n - 1]#gen_exp_vec(n + 1, m*d - n - 1, params.termorder)

    if cache.vars_reversed
        reverse!.(ev)
    end

    row_monomials = gen_mon(ev,PR)

    if cache.vars_reversed
        reverse!.(ev)
    end

    M = compute_basis_matrix(f, d*m - n - 1, m, params, cache)
    if isempty(M)
        return row_monomials
    end

    rows = size(M)[1]

    temp = size(M)

    if (0 < params.verbose)
        println("Finding pivot rows of matrix of size $(size(M))")
        @time pivot_rows = pivot_columns(transpose(M))
    else
        pivot_rows = pivot_columns(transpose(M))
    end
    non_pivot_rows = setdiff([1:rows;], pivot_rows)
    return map((i) -> row_monomials[i], non_pivot_rows)
end

# Computes the the monomial bases for different `m`. That is,
# `compute_monomial_bases(f, R, PR)[m]` will give the `m`-th case.
function compute_monomial_bases(f, params, cache)
    n = nvars(parent(f)) - 1

    res = []

    for m in 1:n
        push!(res, compute_monomial_basis(f, m, params, cache))
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
termorder - the order of the monomials used for vectors

I think these are correct: (TODO)
R - coefficient_ring(parent(f))
PR- paren(f)
"""
function pseudo_inverse_controlled(f, S, l, R, PR, params, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    fLift = liftCoefficients(ZZ,PRZZ,f)
    #println("fLift=$fLift")

    if (0 < params.verbose)
        println("Computing controlled matrix...")
        @time U = compute_controlled_matrix(fLift, l, S, ZZ, PRZZ, params, cache)
    else
        U = compute_controlled_matrix(fLift, l, S, ZZ, PRZZ, params, cache)
    end
    #U = compute_controlled_matrix(f, l, S, R, PR, params)
    
    (6 < params.verbose) && println("controlled matrix: \n$U")
  
    #temp = size(U)
    if (0 < params.verbose)
        println("Computing pseudoinverse of matrix of size $(size(U))")
        @time flag, B = is_invertible_with_inverse(matrix(R,[R(x) for x in Array(U)]), side=:right)
    else
        flag, B = is_invertible_with_inverse(matrix(R,[R(x) for x in Array(U)]), side=:right)
    end
    
    (6 < params.verbose) && println("pinv mod p: \n$B")

    if flag
        return (U,B)
    else 
        if S == collect(0:n) && l == ((d-2)*(n+1)+1)
            throw(ArgumentError("f is not smooth"))
        else
            throw(ArgumentError("matrix from f is not right invertible"))
        end 
    end
end

"""
    pseudo_inverse_controlled_lifted(f,S,l,M,params,cache)

Solves the linear algebra problem as in
`pseudo_inverse_controlled`, but then 
hensel lifts to Z/p^MZ

f - the polynomial definitng the hypersurface
S - the set in [0..n] to be used for the linear algebra problem
l - ???? we need to document this, it's something used by compute_contolled_matrix
M - the absolute precision to lift to.
params - 
cache - 
"""
function pseudo_inverse_controlled_lifted(f,S,l,M,params,cache)
    PR = parent(f)
    R = coefficient_ring(PR)
    (U, Sol_fp) = pseudo_inverse_controlled(f,S,l,R,PR,params,cache)
    lift_to_int(s) = map(x -> lift(ZZ,x),s)

    if (0 < params.verbose)
        println("Lifting mod p solution to the integers")
        @time Sol_mod_p_int = lift_to_int(Sol_fp)
    else
        Sol_mod_p_int = lift_to_int(Sol_fp)
    end
    #U_int = lift_to_int64(U)

    #println("Solution mod p: $Sol_fp")
    #println("U lifted: $U_int")

    p = characteristic(PR)

    if (0 < params.verbose)
        println("Hensel lifting matrix of size $(size(Sol_fp))")
        @time henselLift(p,M,U,Sol_mod_p_int)
    else
        henselLift(p,M,U,Sol_mod_p_int)
    end
end

## Computes the pseudo_inverse for the classical case.
##TODO: update this to reflex changes to pseudo_inverse_controlled
##it's used in standard reduction, I'll plan to take care of it then
#function pseudo_inverse_classical(f, R, PR)
#    return pseudo_inverse_controlled(f, [i for i in 1:n+1], R, PR)
#end
#
#function pseudo_inverse_classicalm(f, m, R, PR)
#    return pseudo_inverse_controlled(f, [i for i in 1:n+1], R, PR)
#end

#end

module AutomatedScript

using Oscar

##################################################################
# Given a homogenous, irreducible polynomial in Q_p[x1, ... xn],
# this script will generate the basis vectors for the n-th de
# Rham cohomology of (X, Z)/Q_p by following the Griffiths-Dwork
# construction algorithm.
##################################################################

# TODO: Modify all the below to work in Fp.

#= Example
n = 2
d = 5
p = 41
Fp = GF(p)

R = Fp
PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
polynomial = Vars[1]^5 + Vars[2]^5 + Vars[3]^5 + Vars[1]*(Vars[2]^3)*Vars[3]
partials = [ derivative(polynomial, i) for i in 1:(n+1) ]
=#

# TODO: Check that it is homogenous. Let d = degree.

# TODO: Ensure partials are not all zero.

"""
    gen_exp_vec(n, d, order)

Returns all nonnegative integer lists of length n who entires sum to d

These are the exponent vectors for all the homogeneous monomials of
degree d, in n variables.

TODO: give this function @memoize. perhaps for some big examples the 
storage required to store the result isn't worth it. However, for 
tests where we're running many similar examples of medium size,
I think it'll be nicer.

INPUTS: 
* "n" -- integer
* "d" -- integer
* "order" -- string, monomial ordering, defaulted to lexicographic ordering. Also supports neglex 
"""
function gen_exp_vec(n, d, order="lex")
    @assert (n >= 0) && (d >= 0) "n and d need to be non-negative"
    result = Any[]
    if n == 1
        return [[d]]
    end

    if order == "lex"
        if d == 1
            for i in 1:n
                s = zeros(Int64,n)
                s[end-i+1] = 1
                push!(result,s)
            end
            return result
        end

        for i in 0:d
            y = gen_exp_vec(n-1,d-i,order)
            for j in axes(y,1)
                prepend!(y[j],i)
            end
            append!(result,y)
        end

    elseif order == "neglex"
        if d == 1
            for i in 1:n
                s = zeros(Int64,n)
                s[i] = 1
                push!(result,s)
            end
            return result
        end

        for i in 0:d
            y = gen_exp_vec(n-1,d-i,order)
            for j in axes(y,1)
                prepend!(y[j],i)
            end
            prepend!(result,y)
        end

    elseif order == "invlex"
        vecs = gen_exp_vec(n,d,"lex")
        result = reverse.(vecs)
    end 

    result
end

function gen_mon(exp_vec, R, PR)
    result = []
    for i in axes(exp_vec,1)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), exp_vec[i])
        monomial = finish(B)
        push!(result,monomial)
    end
    result
end

# Computes all monomials of degree `d` in `n` variables.
function compute_monomials(n, d)
    S, vars = polynomial_ring(ZZ, ["x$i" for i in 1:n])

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

    result
end

# Computes the relations
function compute_relations(monomials, partials)
    result = []
    for i in axes(monomials,1)
        for j in axes(partials,1)
            push!(result, monomials[i]*partials[j])
        end
    end
    result
end

# Converts vector of homogeneous polynomials to a matrix of their coefficents
function convert_p_to_m(polys, expvec)
    result = []
    for i in axes(polys,1)
        temp = []
        for j in axes(expvec,1)
            push!(temp, coeff(polys[i], expvec[j]))
        end
        push!(result, transpose(temp))
    end
    vcat(result...)
end

# Converts Matrix of coefficents to vector of polynomials, each row is one polynomial
function convert_m_to_p(mat, expvec, R, PR)
    result = []
    for i in axes(mat,1)
        B = MPolyBuildCtx(PR)
        for j in axes(expvec,1)
            push_term!(B, mat[i,j], expvec[j])
        end
        push!(result,finish(B))
    end
    result
end

# Computes the basis vectors associated with case h. Columns of
# returned matrix will be linearly independent vectors.
function basis_vectors(n, d, p, precision, polynomial, R, PR, order="lex")
    Qp = PadicField(p,precision)
    result = []
    partials = [ derivative(polynomial, i) for i in 1:(n+1) ]
    # If number of monomials is too small, just use
    # constant as basis vector.
    for h in 1:n
        if h*d - n - 1 <= 0
            push!(result,[PR(1), h]) 
        else
            # compute all monomials of degree `hd - n - 1`
            expvec = gen_exp_vec(n+1, h*d - n - 1, order)

            if h*d - n - d <= 0
                B = gen_mon(expvec,R,PR)
            else
                # TODO: compute all distinct products between the partial derivatives
                # and monomials(n+1, hd - n - d). These are our relations.
                #rmonomials = compute_monomials(n+1, h*d - n - d)
                rexpvec = gen_exp_vec(n+1,h*d - n - d, order)
                rmonomials = gen_mon(rexpvec,R,PR)
                relations = compute_relations(rmonomials, partials)

                # TODO: Check that the number of relations is <= the number of
                # monomials (which is = n + d - 1 choose d)

                # TODO: Let the i-th vector element of `monomials` be the i-th basis
                # element in the monomial basis. Convert all relations into that form.
                M = matrix(R, convert_p_to_m(relations, expvec))
                v, N = nullspace(M)
                B = convert_m_to_p(transpose(N),expvec, R, PR)
            end
            Bh = []
            for i in B
                #push!(Bh, [change_base_ring(Qp,map_coefficients(lift,i)),h])
                push!(Bh, [map_coefficients(lift,i),h])
            end
            append!(result, Bh)
        end
    end
    result
end

# TODO: Compute Omega in terms of differential symbol and exterior
# products (perhaps as symbols).

# TODO: For each linearly independent vector m, compute $$m * Omega / P^h$$.
# Add these to the list of basis vectors.

# Future TODO: Parallelize these algorithms either by hand or by trying to use
# more library functions (for example the Combinatorics.jl library seems promising
# and likely is automatically parallel).
end

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
PR, Vars = PolynomialRing(R, ["x$i" for i in 0:n])
polynomial = Vars[1]^5 + Vars[2]^5 + Vars[3]^5 + Vars[1]*(Vars[2]^3)*Vars[3]
partials = [ derivative(polynomial, i) for i in 1:(n+1) ]
=#

# TODO: Check that it is homogenous. Let d = degree.

# TODO: Ensure partials are not all zero.

# Generates vector of exponentials
function gen_exp_vec(n,d)
    result = Any[]
    if n == 1
        return [[d]]
    end
    if d == 1
        for i in 1:n
            s = zeros(Int64,n)
            s[i] = 1
            push!(result,s)
        end
        return result
    end
    for i in 0:d
        y = gen_exp_vec(n-1,d-i)
        for j in axes(y,1)
            append!(y[j],i)
        end
        append!(result,y)
    end
    return result
end

function gen_mon(exp_vec, R, PR)
    result = []
    for i in axes(exp_vec,1)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), exp_vec[i])
        monomial = finish(B)
        push!(result,monomial)
    end
    return result
end

# Computes all monomials of degree `d` in `n` variables.
function compute_monomials(n, d)
    S, vars = PolynomialRing(ZZ, ["x$i" for i in 1:n])

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

# Computes the relations
function compute_relations(monomials, partials)
    result = []
    for i in axes(monomials,1)
        for j in axes(partials,1)
            push!(result, monomials[i]*partials[j])
        end
    end
    return result
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
    return vcat(result...)
end

# Converts Matrix of coefficents to vector of polynomials, each row is one polynomial
function convert_m_to_p(mat, expvec, R, PR)
    result = []
    for i in 1:nrows(mat)
        B = MPolyBuildCtx(PR)
        for j in axes(expvec,1)
            push_term!(B, mat[i,j], expvec[j])
        end
        push!(result,finish(B))
    end
    return result
end

# Computes the basis vectors associated with case h. Columns of
# returned matrix will be linearly independent vectors.
function basis_vectors(n, d, polynomial, R, PR)
    result = []
    partials = [ derivative(polynomial, i) for i in 1:(n+1) ]
    # If number of monomials is too small, just use
    # constant as basis vector.
    for h in 1:n
        if h*d - n - 1 <= 0
            push!(result,[PR(1), h]) 
        else
            # compute all monomials of degree `hd - n - 1`
            expvec = gen_exp_vec(n+1, h*d - n - 1)

            if h*d - n - d <= 0
                B = gen_mon(expvec,R,PR)
            else
                # TODO: compute all distinct products between the partial derivatives
                # and monomials(n+1, hd - n - d). These are our relations.
                #rmonomials = compute_monomials(n+1, h*d - n - d)
                rexpvec = gen_exp_vec(n+1,h*d - n - d)
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
                push!(Bh, [map_coefficients(lift,i),h])
            end
            append!(result, Bh)
        end
    end
    return result
end

# TODO: Compute Omega in terms of differential symbol and exterior
# products (perhaps as symbols).

# TODO: For each linearly independent vector m, compute $$m * Omega / P^h$$.
# Add these to the list of basis vectors.

# Future TODO: Parallelize these algorithms either by hand or by trying to use
# more library functions (for example the Combinatorics.jl library seems promising
# and likely is automatically parallel).
end
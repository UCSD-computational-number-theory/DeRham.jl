using Oscar

##################################################################
# Given a homogenous, irreducible polynomial in Q_p[x1, ... xn],
# this script will generate the basis vectors for the n-th de
# Rham cohomology of (X, Z)/Q_p by following the Griffiths-Dwork
# construction algorithm.
##################################################################

# TODO: Modify all the below to work in Fp.

n = 2
p = 2
Fp = GF(p)

R, (x, y, z) = PolynomialRing(ZZ, ["x$i" for i in 1:n])
polynomial = y^2*z - x^3 - x*z^2 - z^3

# TODO: Check that it is homogenous. Let d = degree.

d = 3

partials = [ derivative(polynomial, i) for i in 1:n ]

# TODO: Ensure partials are not all zero.

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

# Computes the basis vectors associated with case h. Columns of
# returned matrix will be linearly independent vectors.
function basis_vectors(h)
    # If number of monomials is too small, just use
    # constant as basis vector.
    if h*d - n - 1 <= 0
        return R(1)

    # compute all monomials of degree `hd - n - 1`
    monomials = compute_monomials(n + 1, h*d - n - 1)

    # TODO: compute all distinct products between the partial derivatives
    # and monomials(n+1, hd - n - d - 1). These are our relations.

    # TODO: Check that the number of relations is <= the number of
    # monomials (which is = n + d - 1 choose d)

    # TODO: Let the i-th vector element of `monomials` be the i-th basis
    # element in the monomial basis. Convert all relations into that form.
    vectors = []

    # Compute the nullspace of the matrix created via the above row vectors
    M = vcat(vectors...)

    return nullspace(M)

# TODO: Compute Omega in terms of differential symbol and exterior
# products (perhaps as symbols).

# TODO: For each linearly independent vector m, compute $$m * Omega / P^h$$.
# Add these to the list of basis vectors.

# Future TODO: Parallelize these algorithms either by hand or by trying to use
# more library functions (for example the Combinatorics.jl library seems promising
# and likely is automatically parallel).
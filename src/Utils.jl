
"""
   gen_exp_vec(n, d, order; vars_reversed=false)

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
function gen_exp_vec(n, d, order::Symbol=:lex; vars_reversed=false)
    result = Vector{Vector{Int64}}(undef, binomial(n+d-1,d))
    for i in 1:binomial(n+d-1,d)
        result[i] = zeros(Int64,n)
    end
    if order == :lex
        for i in 1:n
            dtemp = copy(d)
            k = 0
            while k <= (length(result) - 1)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[length(result)-k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[length(result)-k][i] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[length(result)-(k+j-i)][j] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[length(result)-(j+k-1)][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[length(result)-(j+k-1)][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    elseif order == :neglex
        for i in 1:n
            dtemp = copy(d)
            k = 1
            while k <= length(result)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[k][i] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[k+j-i][j] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[j+k-1][i] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    elseif order == :invlex
        for i in 1:n
            dtemp = copy(d)
            k = 0
            while k <= (length(result) - 1)
                if i > 1
                    dtemp = copy(d)
                    for j in 1:n
                        dtemp = dtemp - result[length(result)-k][j]
                    end
                end
                if i == n && dtemp > 0
                    result[length(result)-k][n-i+1] = dtemp
                    k = k + 1
                    continue
                end
                if dtemp == 0
                    k = k + 1
                    continue
                end
                if dtemp == 1 && i > 1
                    for j in i:n
                        result[length(result)-(k+j-i)][n-j+1] = 1
                    end
                    k = k + n - i + 1
                    continue
                end
                dtemp2 = copy(d - dtemp)
                if i == 1 || dtemp2 == 0
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-1,d-dtemp)
                            result[length(result)-(j+k-1)][n-i+1] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-1,d-dtemp)
                        dtemp = dtemp - 1
                    end
                else
                    while dtemp >= 0
                        for j in 1:binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                            result[length(result)-(j+k-1)][n-i+1] = dtemp
                        end
                        k = k + binomial(n-i+d-dtemp-dtemp2-1,d-dtemp-dtemp2)
                        dtemp = dtemp - 1
                    end
                end
            end
        end
    else
        throw(ArgumentError("Unsupported order '$order'"))
    end
    
    if vars_reversed
        for i in 1:binomial(n+d-1,d)
            reverse!(result[i])
        end
    end
    return result
end

"""
FIXME/DOCUMENTME: exp_vec seems to be a 2d array here?

It will also work if exp_vec is an array of arrays.
"""
function gen_mon(exp_vecs, PR)
    result = zeros(PR,length(exp_vecs))
    R = base_ring(PR)
    B = MPolyBuildCtx(PR)
    for i in axes(exp_vecs,1)
        push_term!(B, one(R), exp_vecs[i])
        result[i] = finish(B)
    end
    result
end

"""
Computes the monomials in n variables, of degree d, in the
variable order order, as Oscar polynonmials.
"""
function compute_monomials(n,d,PR,order=:lex,cache=nothing; vars_reversed=false)
    if n < 0 || d < 0
        return []
    end
    if cache != nothing
        exp_vecs = cache[d]
    else
        exp_vecs = gen_exp_vec(n,d,order)
    end

    if vars_reversed
        reverse!.(exp_vecs)
        res = gen_mon(exp_vecs,PR)
        reverse!.(exp_vecs)
        res
    else
        gen_mon(exp_vecs,PR)
    end
end

"""
    polynomial_to_vector(f, n, R, PR; order=:lex)

Convert (homogeneous) polynomial to vector form with specified order. Default is lexicographic.

f - an oscar polynomial
n - the number of variables (minus 1???)
order - a symbol which denotes the term order
cache - a GradedExpCache which gives us the variables
vars_reversed - should we reverse the variables?
"""
function polynomial_to_vector(f, n, order=:lex,cache=nothing; vars_reversed=false)

    #TODO: remove vars_reversed and set it based on the cache's property

    #TODO: if f turns out to be zero, we don't know what the degree should be.
    #
    #How best to fix this?
    R = base_ring(parent(f))
    d = total_degree(f)

    if cache != nothing
        mon_vec = cache[d]
        #println(cache[d])
    else
        mon_vec = gen_exp_vec(n,d,order,vars_reversed=vars_reversed)
    end
    #if vars_reversed 
    #    mon = gen_mon([reverse(tmp) for tmp in mon_vec], R, PR)
    #else 
    #    mon = compute_monomials(n, d, PR, order)
    #end 

    #mon = compute_monomials(n, d,PR,order)
    res = fill(R(0), length(mon_vec))
    for i in eachindex(mon_vec)
    #    println("interpreting vector $(mon_vec[i])")
        if vars_reversed
            reverse!(mon_vec[i])
    #        println("changed to $(mon_vec[i])")
            res[i] = coeff(f, mon_vec[i])
    #        println("coeff: $(res[i])")
            reverse!(mon_vec[i])
        else
            res[i] = coeff(f, mon_vec[i])
        end
    end

    #error()
    #println(res)
    res
end

# 
"""
vector_to_polynomial(vect, n, d, PR, order=:lex)

Convert vector to polynomial with specified order. Default is lexicographic.

vect - vector of coefficients
n - number of variables-1
d - homogeneous degree
PR - polynomial ring to be the parent of the return value
order - a symbol which denotes the term order
"""
function vector_to_polynomial(vect, n, d, PR, order=:lex, vars_reversed=false)
    
    C = MPolyBuildCtx(PR)
    R = base_ring(PR)
    exp_vecs = gen_exp_vec(n+1,d,order)
    @assert length(vect) == length(exp_vecs) "vector has incorrect length for the specified degree"
    for i in eachindex(vect)
        #res += PR(vect[i]) * mon[i]
        v = vect[i]
        if v isa AbstractFloat
            v = Integer(vect[i])
        end

        if vars_reversed
            push_term!(C, R(v), reverse(exp_vecs[i]))
        else 
            push_term!(C, R(v), exp_vecs[i])
        end
    end
    res = finish(C)

    res
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
function convert_p_to_m(polys, expvecs; vars_reversed=false)
    R = coefficient_ring(parent(polys[1]))
    MS = matrix_space(R, length(polys), length(expvecs))
    result = MS()
    for i in axes(polys,1)
        for j in axes(expvecs,1)
            exp_vec = expvecs[j]

            if vars_reversed
                reverse!(exp_vec)
            end
                
            result[i,j] = coeff(polys[i], exp_vec)

            if vars_reversed
                reverse!(exp_vec)
            end
        end
    end
    return result
    #=
    result = []
    for i in axes(polys,1)
        temp = []
        for j in axes(expvec,1)
            push!(temp, coeff(polys[i], expvec[j]))
        end
        push!(result, transpose(temp))
    end
    vcat(result...)
    =#
end

# Converts Matrix of coefficents to vector of polynomials, each row is one polynomial
function convert_m_to_p(mat, expvec, R, PR)
    result = []
    B = MPolyBuildCtx(PR)
    for i in axes(mat,1)
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

# MARK - stuff copied from the old utils file

# Given an Oscar Matrix M, return a the indices of those columns whcih are pivots after being particular
# into reduced row echelon form.
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

"""
    liftCoefficients(R, PR, f)

Lifts the coefficeints of f to the ring R.

Works by lifting coefficients to ZZ and then
converting elements of ZZ to elements of R.
Thus, this method is mostly useful for
lifting something mod p^m to p^n,
for m < n.

INPUTS:
* "f" -- the polynomial to be lifted
* "R" -- the ring for the coefficients to end up in
* "PR" -- the polynomial ring (over R) for the result to end up in
"""
function liftCoefficients(R, PR, f, positiveLift=true)
    t = terms(f)
    sum = 0
    for i in t
        ev = exponent_vector(i,1)
        c = coeff(i,1)
        B = MPolyBuildCtx(PR)
        charBaseField = characteristic(parent(f))
        if positiveLift && (lift(ZZ,c) > div(charBaseField, 2))
            push_term!(B, R(lift(ZZ,c)-charBaseField), ev)
        else
            push_term!(B, R(lift(ZZ,c)), ev)
        end
        sum = sum + finish(B)
    end
    return sum
end

function Factorial(x,y)
    fact = 1
    if x == 0
        return 1
    end
    #BEWARE could be infinite
    while x != y
        fact = fact*x
        x = x - 1
    end
    return fact
end

"""
Performs a mathematical mod (in julia if you mod a negative number, it stays negative). 
Was using for debugging.

"""
function julia_signed_mod(x,n)
    ((x % n) + n) % n
end

"""
    henselLift(p, precision, A, T)
Hensel lifts mod p solution T to the linear system AX-I=0 to mod p^precision

INPUTS:
* "p" -- integer, a prime number 
* "precision" -- integer 
* "A" -- matrix, integer coefficients
* "T" -- matrix, integer coefficients, satisfies AT-I=0 mod p
"""
function henselLift(p, precision, A, T)
    i = 1

    temp_T = zero_matrix(base_ring(T),size(T)...)
    copy!(temp_T,T)
    #add!(temp_T,temp_T,T)
    temp_2T = zero_matrix(base_ring(T),size(T)...)
    temp_AT = zero_matrix(base_ring(T),size(A,1),size(T,2))
    temp_TAT = zero_matrix(base_ring(T),size(T)...)
    while i < precision
        Oscar.Nemo.mul!(temp_2T,2,temp_T)
        #println("Computed: $temp_2T")
        #actual = 2*T
        #println("Actual: $(2*T)")
        #println("Correct? $(temp_2T == actual)")
        
        Oscar.Nemo.mul!(temp_AT,A,temp_T)
        #println("Computed: $temp_AT")
        #actual = A*T
        #println("Actual: $(A*T)")
        #println("Correct? $(temp_AT == actual)")
        
        Oscar.Nemo.mul!(temp_TAT,temp_T,temp_AT)
        #println("Computed: $temp_TAT")
        #actual = T*(A*T)
        #println("Actual: $(T*(A*T))")
        #println("Correct? $(temp_TAT == actual)")
        
        Oscar.Nemo.sub!(temp_T,temp_2T,temp_TAT)
        #println("Computed: $temp_T")
        #actual = 2*T - T * (A*T)
        #println("Actual: $(2*T - T * (A*T))")
        #println("Correct? $(temp_T == actual)")

        #T = 2*T - T * (A*T)
        #println("After step $i: $(julia_signed_mod.(T,p^(i+1)))")
        #println("After step $i: $(julia_signed_mod.(temp_T,p^(i+1)))")

        i *= 2
    end
    R, pi = residue_ring(ZZ, p^precision)
    #stuff = [R(x) for x in Array(T)]
    return matrix(R,temp_T)
end

function findpivotcolumnsU(U)
    n = ncols(U)
    pivotcols = []
    for i in 1:nrows(U)
        for j in 1:ncols(U)
            if U[i,j] != 0
                push!(pivotcols,j)
                break
            end
        end
    end
    return pivotcols
end

function Linverse(L)
    MS = matrix_space(parent(L[1,1]),nrows(L),ncols(L))
    A = MS()
    for j in 1:ncols(L)
        for i in i:nrows(L)
            if i == j 
                A[i,j] = inverse(L[i,j])
            else
                A[i,j] = inverse(L[i,j])*sum(-L[i,k]*A[k,j] for k in 1:(i-1))
            end
        end
    end
end

function LUpseudoinverse(L,U,P1,P2)
    A = permutecolumns(Linverse(L),P1)
    MS = matrix_space(parent(U[1,1]),ncols(U),nrows(U))
    B = MS()
    for j in 1:ncols(B)
        for i in 1:nrows(U)
            B[nrows(U)-i+1,ncols(B)-j+1] = (A[[nrows(U)-i+1,ncols(B)-j+1] + sum(-U[nrows(U)-i+1,ncols(B)-j+1+k]*B[nrows(U)-i+1+k,ncols(B)-j+1] for i in 1:(nrows(B)-ncols(B)))*inverse(U[nrows(U)-i+1,ncols(B)-j+1])])
        end
    end
    return permuterows(B,inverseperm(P2))
end

function permuterows(M,P)
    A = copy(M)
    for i in 1:length(P)
        A[P[i],:] = M[i,:]
    end
    return A
end

function permutecolumns(M,P)
    A = copy(M)
    for i in 1:length(P)
        A[:,P[i]] = M[:,i]
    end
    return A
end

function inverseperm(P)
    S = copy(P)
    for i in 1:length(P)
        S[P[i]] = i
    end
    return S
end


#end

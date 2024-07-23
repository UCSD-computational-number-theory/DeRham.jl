module Utils

using Oscar

macro assert(ex)
    return :( $ex ? nothing : throw(AssertionError($(string(ex)))) )
end

"""
polynomial_to_vector(f, n, R, PR; order=:lex)

Convert polynomial to vector form with specified order. Default is lexicographic.

f is an oscar polynomial
n is the number of variables (minus 1???)
R is the base ring
PR is the polynomial ring, i.e. it is parent(f)
order is a symbol, which must be :lex or the function aborts
"""
function polynomial_to_vector(f, n, R, PR; order=:lex)
    vars = gens(PR)

    if order == :lex
        d = total_degree(f)

        mon = compute_monomials(n, d, PR,vars)
        res = fill(R(0), length(mon))
        for i in eachindex(mon)
            res[i] = coeff(f, mon[i])
        end
        return res
    else
        throw(ArgumentError("Invalid option '$order'"))
    end
end

# Convert vector to polynomial with specified order. Default is lexicographic.
function vector_to_polynomial(vect, n, d, PR, order=:lex)
    if order == :lex
        res = PR()
        mon = compute_monomials(n + 1, d, PR)
        for i in eachindex(vect)
            res += PR(vect[i]) * mon[i]
        end
    else
        throw(ArgumentError("Invalid option '$order'"))
    end
    return res
end

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

# Computes all monomials of degree `d` in `n` variables in the polynomial ring `PR`.
function compute_monomials(n, d, PR, vars, order=:lex)
    if d <= 0
        return []
    end

    if order == :lex
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
    else
        throw(ArgumentError("Invalid option '$order'"))
    end
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

function getTerms(poly)
    t = terms(poly[1])
    result = []
    for i in t
        push!(result,[i,poly[2]])
    end
    return result
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
    henselLift(precision, A, T)
Hensel lifts mod p solution T to the linear system AX-I=0 to mod p^precision

INPUTS: 
* "A" -- matrix, integer coefficients
* "T" -- matrix, integer coefficients, satisfies AT-I=0 mod p
"""
function henselLift(precision, A, T)
    i = 1 
    while i < precision
        T = 2*T - T * (A*T)
        i *= 2 
    end
    return T
end

end

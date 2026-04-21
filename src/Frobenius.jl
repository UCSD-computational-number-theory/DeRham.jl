
"""
    computeD(N, m)

Returns a list of length N where D_{j,m} = sum_{i=j}^{N-1} (-1)^{i+j}binom{-m}{i}binom{i}{j}

INPUTS: 
* "N" -- integer
* "m" -- integer
"""
function computeD(N, m)
    D = zeros(Int,N)
    for j in 0:(N-1)
        D[j+1] = sum((-1)^(i+j)*binomial(-m,i)*binomial(i,j) for i in j:(N-1))
    end
    return D
end

"""
    applyFrobeniusToMon(n,d,f,N,p,beta,m,R,PR)

Computes the power series expansion of p^{m-n-1}sigma(x^{beta}Omega/f^m) 
using formula (1.10) in Costa's thesis


INPUTS: 
* "n" -- integer, dimension of ambient projective space
* "d" -- integer, degree of the hypersurface f
* "f" -- polynomial, defining homogeneous equation of the hypersurface lifted to characteristic 0
* "N" -- integer, series precision
* "p" -- integer, a prime number that is the characteristic of the base field of the hypersurface
* "beta" -- vector, representing the exponents in the monomial of the basis element
* "m" -- integer, pole order of the basis element 
* "R" -- ring, precision ring 
* "PR" -- ring, polynomial ring with coefficients in R 
* "termorder" -- the term ordering that should be used in vector representations
* "vars_reversed" -- reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.
"""
function applyFrobeniusToMon(n, d, f, N, p, beta, m, R, PR, termorder,vars_reversed; verbose=false)
    (9 < verbose) && println("Applying Frobenius to $beta with pole order $m")
    (9 < verbose) && println("N=$N, m=$m")
    if N == 0
        return []
    end
    (9 < verbose) && println("Scaling by factorial of: ", p * (N + m - 1) - 1)
    Factorial = factorial(ZZ(p * (N + m - 1) - 1))
    o = ones(Int64, n+1)
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), o)
    X1 = finish(B)
    D = computeD(N,m)
    result = []
    for j in 0:(N-1)
        e = j + m
        factorial_e = R(ZZ(Factorial//factorial(ZZ(p * e - 1))))
        ev = gen_exp_vec(n+1,d*j,termorder)
        fj = f^j
        sum = 0
        for alpha in ev
            B = MPolyBuildCtx(PR)
            push_term!(B, R(1),Int64(p) * (beta + alpha + o))
            monomial = div(finish(B), X1)
            coefficient = R(factorial_e * (D[j+1] * (coeff(fj,alpha)))) # needs to act by Frobenius on coeff() if base field is not FF_p
            sum = sum + coefficient * monomial
            if (9 < verbose) && coefficient != 0
                Djm = D[j+1]
                C_jalpha = coeff(fj,alpha)
                println("Djm=$Djm, C_jalpha=$C_jalpha\n")
            end
        end
        push!(result, [sum, p*(m+j)])
    end
    return result
end

"""
    applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)

Applies the frobenius to all the elements of Basis

INPUTS: 
* "Basis" -- array of basis elmenets
* "n" -- number of variables minus 1
* "f" -- polynomial which is the denominator of poles (lifted version)
* "N_m" -- series precision
* "p" -- the prime
* "params" - the ControlledReductionParamaters
"""
function applyFrobeniusToBasis(Basis,f,N_m,p,params)
    termorder = params.termorder
    vars_reversed = params.vars_reversed
    verbose = params.verbose

    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    result = []
    for b in Basis
        Fmon = applyFrobeniusToMon(n,d,f,N_m[b[2]],p,exponent_vector(b[1],1),b[2],R,PR,termorder,vars_reversed,verbose=verbose)
        (9 < verbose) && println("Applying Frobenius to $b gives $Fmon")
        push!(result, Fmon)
    end
    return result
end
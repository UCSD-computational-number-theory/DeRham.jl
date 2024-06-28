module Frobenius

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
#include("FindMonomialBasis.jl")
include("Utils.jl")
#include("SmallestSubsetSmooth.jl")
include("StandardReduction.jl")
include("PolynomialWithPole.jl")

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
"""
function applyFrobeniusToMon(n, d, f, N, p, beta, m, R, PR)
    #FIXME reversed to match Costa's code
    beta = reverse(beta)
    println("N=$N, m=$m")
    Factorial = factorial(big(p * (N + m - 1) - 1))
    o = ones(Int64, n+1)
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), o)
    X1 = finish(B)
    D = computeD(N,m)
    result = []
    for j in 0:(N-1)
        e = j + m
        factorial_e = R(ZZ(Factorial/factorial(big(p * e - 1))))
        println("e=$e,factorial_e=$factorial_e")
        ev = Utils.gen_exp_vec(n+1,d*j,:invlex)
        fj = f^j
        sum = 0
        for alpha in ev
            B = MPolyBuildCtx(PR)
            push_term!(B, R(1), p * (beta + alpha + o))
            monomial = div(finish(B), X1)
            sum = sum + R(factorial_e * (D[j+1] * (coeff(fj,alpha)^p))) * monomial
            #println(typeof((D[j+1]*(coeff(map_coefficients(lift,fj),alpha)^p))*monomial))
        end
        push!(result, [sum, p*(m+j)])
    end
    return result
end


#=
function applyFrobenius(n,d,f,N,p,poly,R,PR)
    t = terms(poly)
    temp = []
    for i in t
        ev = exponent_vector(i[1],1)
        push!(temp, applyFrobeniusToMon(n,d,f,N,p,ev,poly[2],R,PR))
    end
    return temp
end
=#

"""
    applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)

Applies the frobenius to all the elements of Basis

INPUTS: 
* "Basis" -- array of basis elmenets
* "n" -- number of variables minus 1
* "d" -- degree
* "f" -- polynomial which is the denominator of poles (lifted version)
* "N" -- series precision
* "p" -- the prime
* "R" -- basering(parent(f))
* "PR" -- parent(f)
"""
function applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)
    result = []
    for b in Basis
        Fmon = applyFrobeniusToMon(n,d,f,N,p,exponent_vector(b[1],1),b[2],R,PR)
        println(Fmon)
        push!(result, Fmon)
    end
    return result
end
end

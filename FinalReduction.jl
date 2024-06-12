module FinalReduction

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

function computeT(Basis,f,n,d,R,PR)
    ev = Utils.gen_exp_vec(n+1,d*n-n-1)
    mons = Utils.gen_mon(ev,R,PR)
    T = []
    for m in mons
        temp = []
        for i in 0:(n-1)
            if i == 0
            else
                mterms = terms(m)
                sum = 0
                for t in mterms
                    sum = sum + StandardReduction.stdRed_step(f,t,n-i+1,1)[1]
                end
                m = inv(R(n-i))*sum
            end
            for j in 0:(length(Basis[n-i])-1)
                c = coeff(PR(m),exponent_vector(Basis[n-i][length(Basis[n-i])-j],1))
                push!(temp, c)
                m = m - c*Basis[n-i][length(Basis[n-i])-j]
            end
        end
        push!(T, transpose(temp))
    end
    return transpose(vcat(T...))
end

end

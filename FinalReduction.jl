module FinalReduction

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("Utils.jl")
#include("SmallestSubsetSmooth.jl")
include("StandardReduction.jl")

function ComputeT(f, Basis, M)
    p = characteristic(parent(f))
    n = nvars(parent(f)) - 1
    d = degree(f,1)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    partials_index = [(n+1)-i for i in 0:n]

    precisionring, = residue_ring(ZZ,p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])
    f_lift = Utils.liftCoefficients(precisionring,precisionringpoly,f)

    exp_vec = Utils.gen_exp_vec(n+1, d*n-n-1, :invlex)
    monomials = Utils.gen_mon(exp_vec, precisionring, precisionringpoly)
    #len = binomial(n+(d*n-n-1), n)
    len = length(monomials)

    partials = [ derivative(f_lift, i) for i in 1:n+1 ]
    partials = reverse(partials)
    println(partials)

    T = []
    for i in 0:n-1
        l = d*(n-i) - n - 1
        if l > 0
            exp_vec = Utils.gen_exp_vec(n+1,l,:invlex)
            if (l-(d-1)) >= 0
                monomials_domain = Utils.compute_monomials(n+1, l-(d-1), precisionringpoly,:invlex)
                len_domain = length(monomials_domain)
                basis = Basis[n-i]
                len_basis = length(basis)
                change_basis,change_basis_inverse = StandardReduction.monomial_change_basis_inverse_lifted(f,l,basis,M)
                change_basis = matrix(precisionring,[precisionring(x) for x in Array(change_basis)])
                tmp = []
                for j in 1:len # indexing over monomials 
                    # row vector for monomials[j] with respect to standard monomial basis
                    row_vec = matrix(precisionring, 1, length(exp_vec), Utils.convert_p_to_m([monomials[j]],exp_vec))
                    vec = change_basis_inverse * transpose(row_vec)
                    println(vec)
                    #print(vec)
                    for k in 1:length(basis)
                        #push!(tmp, precisionring(ZZ(factorial(n-1-i)))*vec[length(exp_vec)-len_basis+k,1])
                        push!(tmp, precisionring(ZZ(factorial(n-1-i)))*vec[k,1])
                    end
                    term = 0
                    vec_test = vec[1,1] * pvars[1]^3 
                    for t in 1:len_domain*(n+1)
                        t_partial = ceil(Int, t/(n+1))
                        t_domain = t%(len_domain)
                        if t_domain == 0
                            t_domain = len_domain
                        end 
                        vec_test = vec_test + vec[t+len_basis,1]*partials[t_partial]*monomials_domain[t_domain]

                        term = term + vec[t+len_basis,1]*derivative(monomials_domain[t_domain],partials_index[t_partial])
                    end 
                    println(vec_test)
                    @assert vec_test == monomials[j]
                    monomials[j] = term
                end 
                pushfirst!(T,tmp)
            end 
        end 
    end 
    pushfirst!(T,monomials)
    return T 
end 
end 

    #=
    T = []
    for mon in monomials
        temp = []
        for i in 0:n-1 
            l = d*(n-i) - n - 1
            if l >= 0
                for j in 0:(length(Basis[n-i])-1)
                    c = coeff(mon,exponent_vector(Basis[n-i][length(Basis[n-i])-j],1))
                    push!(temp, PrecisionRing(ZZ(c*factorial(n-1-i))))
                    mon = mon - c*Basis[n-i][length(Basis[n-i])-j]
                end
            end
            if i != n-1
                # apply standard reduction to terms in monomials[i] to reduce pole order by 1 
                standard_right_inverse = CopiedFindMonomialBasis.pseudo_inverse_controlled_lifted(f,S,l,M)
                target_exp_vec = Utils.gen_exp_vec(n+1, l, :invlex)
                domain_exp_vec = Utils.gen_exp_vec(n+1, l-(d-1), :invlex)
                sum = 0
                for t in mterms
                    sum = sum + StandardReduction.stdRed_step(f,t,standard_right_inverse,n-i,target_exp_vec,domain_exp_vec,M)
                end 
                m = sum  
            end
        push!(T,temp)    
        end
    end
    return transpose(matrix(PrecisionRing, length(monomials), sum([length(Basis)[i] for i in 1:length(Basis)], T)))
end  
=#



    #=
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
=#



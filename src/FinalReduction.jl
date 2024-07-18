#module FinalReduction

#using Oscar
#using BitIntegers
#using LinearAlgebra
#using Combinatorics

#include("Utils.jl")
#include("StandardReduction.jl")

function computeT(f, Basis, M)
    p = characteristic(parent(f))
    n = nvars(parent(f)) - 1
    d = degree(f,1)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    partials_index = [(n+1)-i for i in 0:n]

    precisionring, pi = residue_ring(ZZ,p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])
    f_lift = liftCoefficients(precisionring,precisionringpoly,f)

    exp_vec = gen_exp_vec(n+1, d*n-n-1, :invlex)
    monomials = gen_mon(exp_vec, precisionring, precisionringpoly)
    #len = binomial(n+(d*n-n-1), n)
    len = length(monomials)

    partials = [ derivative(f_lift, i) for i in 1:n+1 ]
    partials = reverse(partials)

    T = zero_matrix(precisionring, 0, len)

    for i in 0:n-1
        l = d*(n-i) - n - 1
        if l > 0
            exp_vec = gen_exp_vec(n+1,l,:invlex)
            if (l-(d-1)) >= 0
                monomials_domain = compute_monomials(n+1, l-(d-1), precisionringpoly,:invlex)
                len_domain = length(monomials_domain)
                basis = Basis[n-i]
                len_basis = length(basis)
                change_basis,change_basis_inverse = monomial_change_basis_inverse_lifted(f,l,basis,M)
                change_basis = matrix(precisionring,[precisionring(x) for x in Array(change_basis)])
                tmp = zero_matrix(precisionring, len_basis, len)
                for j in 1:len # indexing over monomials 
                    # row vector for monomials[j] with respect to standard monomial basis
                    #row_vec = matrix(precisionring, 1, length(exp_vec), convert_p_to_m([monomials[j]],exp_vec))
                    row_vec = convert_p_to_m([monomials[j]],exp_vec)
                    vec = change_basis_inverse * transpose(row_vec)
                    for k in 1:length(basis)
                        #push!(tmp, precisionring(ZZ(factorial(n-1-i)))*vec[length(exp_vec)-len_basis+k,1])
                        #push!(tmp, precisionring(ZZ(factorial(n-1-i)))*vec[k,1])
                        tmp[k,j] = precisionring(ZZ(factorial(n-1-i)))*vec[k,1]
                    end
                    term = 0
                    for t in 1:len_domain*(n+1)
                        #t_partial = ceil(Int, t/(n+1))
                        t_partial = ceil(Int, t/(n+1)) % (n+1)
                        if t_partial == 0
                            t_partial = n+1
                        end
                        t_domain = t%(len_domain)
                        if t_domain == 0
                            t_domain = len_domain
                        end 
                        term = term + vec[t+len_basis,1]*derivative(monomials_domain[t_domain],partials_index[t_partial])
                    end 
                    monomials[j] = term
                end 
                #T = vcat(matrix(precisionring,1,len,tmp),T)
                T = vcat(tmp, T)
            end 
        end 
    end 
    T = vcat(matrix(precisionring,1,len,[trailing_coefficient(x) for x in monomials]), T)
    return T 
end 

#end 



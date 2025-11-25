
"""
    computeT(f,Basis,M,params)

Computes the "T matrix" as in the notation in Costa's thesis.
"""
function computeT(f, Basis, M, params, cache)
    p = characteristic(parent(f))
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    # if params.vars_reversed 
    #     partials_index = [(n+1)-i for i in 0:n]
    # else
        partials_index = [i for i in 1:n+1]
    # end 

    precisionring, pi = residue_ring(ZZ,p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])
    f_lift = liftCoefficients(precisionring,precisionringpoly,f)

    exp_vec = cache[d*n-n-1]#gen_exp_vec(n+1, d*n-n-1, params.termorder)

    # if cache.vars_reversed
    #     reverse!.(exp_vec)
    # end

    monomials = gen_mon(exp_vec, precisionringpoly)

    # if cache.vars_reversed
    #     reverse!.(exp_vec)
    # end

    #len = binomial(n+(d*n-n-1), n)
    len = length(monomials)

    partials = [ derivative(f_lift, i) for i in 1:n+1 ]
    # if params.vars_reversed == true
    #     partials = reverse(partials)
    # end
    #println(partials)

    T = zero_matrix(precisionring, 0, len)

    for i in 0:n-1
        l = d*(n-i) - n - 1
        basis = Basis[n-i]
        len_basis = length(basis)
        if l >=0 
            exp_vec = cache[l]#gen_exp_vec(n+1,l,params.termorder)
            if (l-(d-1)) >= 0
                monomials_domain = compute_monomials(n+1, l-(d-1), precisionringpoly,params.termorder,cache)#, vars_reversed=cache.vars_reversed)
                len_domain = length(monomials_domain)
                change_basis,change_basis_inverse = monomial_change_basis_inverse_lifted(f,l,basis,M,params,cache)
                #println("change basis = $change_basis")
                #println("change_basis_inverse = $change_basis_inverse")
                #change_basis = matrix(precisionring,[precisionring(x) for x in Array(change_basis)])
                tmp = zero_matrix(precisionring, len_basis, len)
                
                for j in 1:len # indexing over monomials 
                    row_vec = convert_p_to_m([monomials[j]],exp_vec)#,vars_reversed=cache.vars_reversed)
                    vec = change_basis_inverse * transpose(row_vec)
                    
                    for k in 1:len_basis 
                        tmp[k,j] = precisionring(ZZ(factorial(n-1-i)))*vec[k,1]
                    end

                    # standard reduction step applying Griffiths-Dwork relation 
                    term = 0
                    for t in 1:len_domain*(n+1)
                        t_partial = ceil(Int, t/len_domain) 
                        t_domain = t%(len_domain)
                        if t_domain == 0
                            t_domain = len_domain
                        end 
                        term = term + vec[t+len_basis,1] * derivative(monomials_domain[t_domain], partials_index[t_partial])
                    end 
                    monomials[j] = term
                end 
                T = vcat(tmp, T)
            else 
                tmp = zero_matrix(precisionring, len_basis, len)
                for j in 1:len # indexing over monomials
                    for k in 1:length(basis)
                        tmp[k,j] = coeff(monomials[j],liftCoefficients(precisionring,precisionringpoly,basis[k]))
                    end
                end 
                T = vcat(tmp, T)
            end 
        end 
    end 
   # T = vcat(matrix(precisionring,1,len,[trailing_coefficient(x) for x in monomials]), T)
    return T 
end 

#end 



"""
    computeT(f,Basis,M,params)

Computes the "T matrix" as in the notation in Costa's thesis.

fields
------
f - polynomial
Basis - output of compute_monomial_bases
M - algorithm precision
params - the ControlledReductionParamaters
cache - the GradedExpCache used for this controlled reduction
"""
function computeT(f, Basis, M, params, cache)
    p = characteristic(parent(f))
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    partials_index = [i for i in 1:n+1]

    precisionring, pi = residue_ring(ZZ,p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])
    f_lift = liftCoefficients(precisionring,precisionringpoly,f)

    exp_vec = cache[d*n-n-1]#gen_exp_vec(n+1, d*n-n-1, params.termorder)

    monomials = gen_mon(exp_vec, precisionringpoly)

    #len = binomial(n+(d*n-n-1), n)
    len = length(monomials)

    partials = [ derivative(f_lift, i) for i in 1:n+1 ]

    T = zero_matrix(precisionring, 0, len)

    for i in 0:n-1
        l = d*(n-i) - n - 1
        basis = Basis[n-i]
        len_basis = length(basis)
        if l >=0 
            exp_vec = cache[l]#gen_exp_vec(n+1,l,params.termorder)
            if (l-(d-1)) >= 0
                monomials_domain = compute_monomials(n+1, l-(d-1), precisionringpoly,params.termorder,cache)
                len_domain = length(monomials_domain)
                change_basis,change_basis_inverse = monomial_change_basis_inverse_lifted(f,l,basis,M,params,cache)
                tmp = zero_matrix(precisionring, len_basis, len)
                
                for j in 1:len # indexing over monomials 
                    row_vec = convert_p_to_m([monomials[j]],exp_vec)
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

    return T 
end 

"""
    reduce_order_n_to_basis

Takes a polynomial with pole of order n, and converts
it to a vector in terms of the basis of cohomology


e - the original pole order of the element before the frobenius was applied
precision - the series precision for the vector
VS - a matrix space for where the vector will end up
poly_to_reduce - a polynomial with pole that will 

"""
function reduce_order_n_to_basis(n,p,e,precision,VS,poly_to_reduce,T,R,cache,ev)

    N = precision

    if N == 0
        return zero_matrix(R, nrows(T), 1)
    end

    ff = factorial(ZZ(p*(e+N-1)-1))
    val_ff = valuation(ff,p)
    final_val = (n-1) - val_ff  
    ff_invertible = ff / p^val_ff

    inverse_ff = inv(R(ff_invertible))

    temp = VS()
    temp2 = convert_p_to_m([poly_to_reduce[1]],ev,vars_reversed=cache.vars_reversed)
    for i in 1:length(ev)
        temp[i,1] = R(temp2[i])
    end
    temp = T * temp
    for i in 1:length(temp)
        ele = inverse_ff * temp[i]
        if 0 ≤ final_val
            ele *= p^final_val
        else
            lifted = lift(ZZ,ele)
            ele = divexact(ele,p^(-final_val)) 
        end
        temp[i] = ele
    end 
    
    temp

end

"""
    compute_frobenius_matrix(n,d,Reductions,T)

Computes Frobenius Matrix

INPUTS: 
* "n" -- dimension of ambient projective space 
* "d" -- degree of the polynomial f2
* "N" -- integer, series precision
* "Reductions" -- output of computeReductionOfTransformLA
* "T" -- output of computeT
* "Basis" -- ?????
* "params" -- the ControlledReductionParamaters
* "cache" -- the GradedExpCache used for this controlled reduction
"""
function compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, params, cache)
    verbose = params.verbose
    termorder = params.termorder
    (9 < verbose) && println("Terms after controlled reduction: $Reductions")
    R = parent(T[1,1])
    frob_mat_temp = []
    denomArray = []
    ev = cache[d*n-n-1]#gen_exp_vec(n+1,d*n-n-1,termorder)
    VS = matrix_space(R,length(ev),1)
    for i in 1:length(Reductions)
        e = Basis[i][2] # pole order of basis element 
        N = N_m[e]
        poly_to_reduce = Reductions[i][1]
        temp = reduce_order_n_to_basis(n,p,e,N,VS,poly_to_reduce,T,R,cache,ev)
        
        push!(frob_mat_temp, temp)
    end
    frob_mat = hcat(frob_mat_temp...)

    return frob_mat
end

#end 



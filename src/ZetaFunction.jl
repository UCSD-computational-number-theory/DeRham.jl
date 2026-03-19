
struct ZetaFunctionParams
    verbose::Int
    givefrobmat::Bool
    algorithm::Symbol
    termorder::Symbol
    vars_reversed::Bool
    fastevaluation::Bool
    always_use_bigints::Bool
    use_gpu::Bool
    use_threads::Bool
end

default_params() = ZetaFunctionParams(0,false,:default,:invlex,false,true,false,false,false)

"""
    controlled_reduction_cache

Give the minimal PolyExpCache needed for controlled reduction
as in Prop 1.15 of Costa's thesis.

n in this function is the number of variables minus one,
i.e. the weight of the affine Monsky-Washnitzer homology module.

the d we need are

FORWARD:
h*d - n - 1 for h in 1:n (if such values are positive)
d*n - n 
d*n - n - 1
d*n - n - length(S),
d*n - n - length(S) + 1,
d*n - n + d - length(S) + 1, #TODO: is this +1 wrong??
d

REVERSE:
d
n*d - n

fields
------
n - int
d - int
S - vector of ints
params - the ControlledReductionParamaters
"""
function controlled_reduction_cache(n,d,S,params)
    degsforward = [d*n - n,
                   d*n - n - 1,
                   d*n - n - length(S),
                   d*n - n - length(S) + 1,
                   d*n - n + d - length(S),
                   d]
    for h in 1:n
        if 0 ≤ h*d - n - 1
            append!(degsforward,h*d - n - 1)
        end
        if 0 ≤ h*d - n - d
            append!(degsforward,h*d - n - d)
        end
    end

    degsreverse = [d,
                   n*d - n,
                   d*n - n + d - length(S)]

    PolyExpCache(n+1,
                 params.termorder,
                 degsforward,
                 degsreverse,
                 vars_reversed=false)

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

"""
    LPolynomial(FM, q)

Given the Frobenius matrix, computes the corresponding L-polynomial det(1-tq^{-1}FM)

INPUT: 
* "FM" -- Frobenius matrix 
"""

function LPolynomial(FM, n, q, polygon, relative_precision, verbose)
    @assert size(FM, 1) == size(FM, 2) "FM is not a square matrix"

    P, T = polynomial_ring(ZZ, "T")
    lift_to_int(s) = map(x -> lift(ZZ,x),s)

    f = charpoly(P, lift_to_int(FM))
    cp_coeffs = collect(coefficients(f))
    return compute_Lpolynomial(n, q, polygon, relative_precision, cp_coeffs, verbose)
    

    """
    k = degree(f)
    bound = binomial(k,Int(ceil(k/2)))*(q^(n/2))*ceil(k/2)
    R, t = polynomial_ring(ZZ,"t")
    result = R(0)
    for i in 0:k
        if abs(ZZ(coeff(f,i))) > bound
            result = result + (ZZ(coeff(f,i)) - characteristic(P))*t^(k-i)
        else
            result = result + (ZZ(coeff(f,i)))*t^(k-i)
        end
    end
    """

end 

function get_basis_of_cohomology_twoflavors(f,S,params,cache)
    n = nvars(parent(f)) - 1

    basis = compute_monomial_bases(f, params, cache) # basis of cohomology 

    if basis == nothing
        (0 < verbose) && println("Cannont compute monomial basis, this f appears to be non-smooth")
        return false
    end

    Basis = []
    for i in 1:n
        for j in basis[i]
            push!(Basis,[j,i])
        end
    end

    (basis,Basis)
end

function get_basis_of_cohomology(f,S,params,cache)
    get_basis_of_cohomology_twoflavors(f,S,params,cache)[2]
end






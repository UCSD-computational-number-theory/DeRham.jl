
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

"""
   precision_information(f,basis)

returns all the precision information
related to f

basis -- the basis of cohomology, in the format of 
         polynomials with pole
"""
function precision_information(f,basis,verbose=0)
    p = Int64(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    hodge_polygon = hodgepolygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths

    k = sum(hodge_numbers) # dimension of H^n
    (0 < verbose) && println("There are $k basis total elements in H^$n to be reduced")

    r_m = calculate_relative_precision(hodge_polygon, n-1, p) 
    N_m = series_precision(p,n,d,r_m)
    N_m[N_m .== 0] .= 1 #TODO: understand the meaning of zero precision better
    M = algorithm_precision(p,n,d,r_m,N_m)

    (hodge_polygon,r_m,N_m,M)
end

"""
    precision_information(f)

returns the precision information for a polynomial f

"""
function precision_information(f)
    n = nvars(parent(f)) - 1
    PR = parent(f)
    R = coefficient_ring(parent(f))

    params = default_params()

    d = total_degree(f)
    S = collect(0:n) #collect(0:d-1)
    cache = controlled_reduction_cache(n,d,S,params)
    basis = compute_monomial_bases(f, params, cache) # basis of cohomology 
    Basis = []
    for i in 1:n
        for j in basis[i]
            push!(Basis,[j,i])
        end
    end

    precision_information(f,Basis)
end

"""
    print_precision_info(n,d,p)

Calcuates how much precision would be needed to compute
the zeta functino of a polynomial with
n+1 variables, of degree d, at the prime p
"""
function print_precision_info(n,d,p)
    R, vars = polynomial_ring(GF(p),n+1)

    # any hypersurface will do, all we care
    # about is the hodge polygon
    fermat_hypersurface = sum(vars .^ d)

    (hp, rel, ser, alg) =  precision_information(fermat_hypersurface)

    println("Hodge Numbers: $(hp.slopelengths)")
    println("Relative precision: $rel")
    println("Series precision: $ser")
    println("Algorithm precision: $alg")
end 

function pregen_precision_info(n,d,p)
    R, vars = polynomial_ring(GF(p),n+1)

    # any hypersurface will do, all we care
    # about is the hodge polygon
    fermat_hypersurface = sum(vars .^ d)

    (hp, rel, ser, alg) =  precision_information(fermat_hypersurface)

    return alg
end 
"""
    zeta_coefficients(f; verbose=false, givefrobmat=false, algorithm=:costachunks, termorder=:invlex, vars_reversed=true)


Wrapper function that outputs the zeta function of
the projective hypersruface defined by `f`.

INPUTS: 
* "f" -- Oscar polynomial (should be homogeneous)

KEYWORD ARGUMENTS:
verbose -- prints various diagnostic statements based on the level
    verbose levels are by convention a number between 0 and 10.
    0: print nothing
    1: print basic diagnostic info only
    2-9: TBD (To be determined/documented)
    2: more basic diagnostic info
    3: print the output of controlled reduction
    5: print out the u and v information in controlled reduction
    6: print the value of g at each step of reduction
    7: print out one of the R_uv matrices (you might need to modify the print statement to fit your example right now)
    10: print anything that we might consider useful
givefrobmat -- should the funciton also output the appoximated frobenius matrix
algorithm -- the algorithm used for controlled reduction
termorder -- the term ordering that should be used in vector representations
fastevaluation -- should the algorithm use fast evaluation?
>>>if you don't know what this is, ignore it.
vars_reversed -- reverses the order of basis vectors at various places
>>>if you don't know what this is, ignore it.

"""
function zeta_coefficients(f; S=[-1], verbose=0, changef=true, givefrobmat=false, algorithm=:default, termorder=:invlex, vars_reversed=false, fastevaluation=true, always_use_bigints=false, use_gpu=false, use_threads=false, context=nothing)
    PR = parent(f)
    R = coefficient_ring(PR)
    p = Int64(characteristic(PR))
    q = p
    n = nvars(PR) - 1
    d = total_degree(f)
    
    if S == [-1]
        S = collect(0:n)
    end

    if length(S) > d
        error("Length of S must be <= $d")
    end

    if algorithm == :default
        if S == [n]
            algorithm = :varbyvar
        elseif n > 4
            for i in 0:n
                if is_Ssmooth(f,[n-i])
                    algorithm = :varbyvar
                    S = [n]
                    vars = gens(PR)
                    new_vars = copy(vars)
                    for j in 0:n
                        if j == i
                            new_vars[i+1] = vars[n+1]
                        end
                    end
                    new_vars[n+1] = vars[i+1]
                    break
                elseif i == n
                    algorithm = :depthfirst
                end
            end
        else
            algorithm = :depthfirst
        end
    end

    if algorithm==:varbyvar && (S != [n])
        throw("S must be [$n] for varbyvar")
    end

    # vars_reversed = false
    params = ZetaFunctionParams(verbose,givefrobmat,algorithm,termorder,vars_reversed,fastevaluation,always_use_bigints,use_gpu,use_threads)

    cache = controlled_reduction_cache(n,d,S,params)

    (0 < verbose) && println("p = $p")
    (9 < verbose) && println("Working with a degree $d hypersurface in P^$n")

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
    
    #println("Basis of cohomology is $Basis")
    (9 < verbose) && println("Basis of cohomology is $Basis")

    (hodge_polygon, r_m, N_m, M) = precision_information(f,Basis,verbose)

    (1 < verbose) && println("N_m=$N_m")

    (9 < verbose) && println("We work modulo $p^$M, and compute up to the $N_m-th term of the Frobenius power series")

    if always_use_bigints || BigInt(2)^64 < BigInt(p)^M
        residue = BigInt(p)^M
    elseif BigInt(2)^63 < BigInt(p)^M
        # If we use UInt, we get things between 2^63 and 2^64
        residue = UInt(p)^M
    else
        residue = p^M
    end

    # FOR DEBUGGING
    #residue = BigInt(p)^M

    precisionring, = residue_ring(ZZ, residue)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])

    if (0 < verbose)
        println("Computing the relations matrix...")
        @time f_changed, f, pseudo_inverse_mat_new = find_Ssmooth_model(f, M, S, params, changef, cache)
    else
        f_changed, f, pseudo_inverse_mat_new = find_Ssmooth_model(f, M, S, params, changef, cache)
    end

    # For large examples, force a GC here; Julia's GC doesn't seem to automatically cause the
    # garbage collector to run if GPU allocations become too much.
    if 5 <= n
        if (0 < verbose)
            println("collecting garbage...")
            @time GC.gc()
        else
            GC.gc()
        end
    end


    if (f == false)
        (0 < verbose) && println("f is not smooth (or not S-smooth) and we're done. ")
        return false
    end 

    (9 < verbose) && println("pseudo_inverse_mat is $pseudo_inverse_mat_new")

    # recomputes basis if f is different 
    if f_changed 
        (0 < verbose) && println("New model is $f")
        basis = compute_monomial_bases(f, params, cache) # basis of cohomology 
        Basis = []
        for i in 1:n
            for j in basis[i]
                push!(Basis,[j,i])
            end
        end
        (9 < verbose) && println("New basis of cohomology is $Basis")
    end 
    
    for i in 1:length(Basis)
        Basis[i][1] = liftCoefficients(precisionring, precisionringpoly, Basis[i][1])
    end 

    (9 < verbose) && println("Basis of cohomology is $Basis")


    fLift = liftCoefficients(precisionring, precisionringpoly, f)
    if (0 < verbose)
        println("\nApplying Frobenius to basis elements")
        @time FBasis = applyFrobeniusToBasis(Basis,fLift, N_m, p, params)
    else
        FBasis = applyFrobeniusToBasis(Basis,fLift, N_m, p, params)
    end

    l = d * n - n + d - length(S)


    MS = matrix_space(precisionring, nrows(pseudo_inverse_mat_new), ncols(pseudo_inverse_mat_new))
    pseudo_inverse_mat = MS()

    if (0 < verbose)
        println("Computing T matrix...")
        @time T = computeT(f, basis, M, params, cache)
    else
        T = computeT(f, basis, M, params, cache)
    end
    (9 < verbose) && println("T matrix is $T")

    for i in 1:nrows(pseudo_inverse_mat_new)
        for j in 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i,j] = ZZ(pseudo_inverse_mat_new[i,j])
        end
    end

    #TODO: check which algorithm we're using
    (0 < verbose) && println("\nStarting controlled reduction...")
    Reductions = reducetransform(FBasis, N_m, S, fLift, pseudo_inverse_mat, p,  params, cache,context) 
    (2 < verbose) && println(Reductions)
    
    ev = cache[n*d - n - 1] 
    (9 < verbose) && println(convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev))
    FM = compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, params, cache)
    (9 < verbose) && println("The Frobenius matrix is $FM")

    if givefrobmat
        (FM,LPolynomial(FM,n,q,hodge_polygon,r_m, verbose))
    else
        LPolynomial(FM,n,q,hodge_polygon,r_m, verbose)
    end
end

"""
a wrapper to zeta_coefficients

INPUTS: 
* "f" -- Oscar polynomial (should be homogeneous) over the integers
* "p" -- a prime number, integer 

"""
function zeta_coefficients(f, p; S=[-1], verbose=0, changef=true, givefrobmat=false, algorithm=:default, termorder=:invlex, vars_reversed=false, fastevaluation=false, always_use_bigints=false, use_gpu=false, use_threads=false, context=nothing)
    @assert is_prime(p) "p must be prime"
    PR = parent(f)
    PRmodp, hom = change_base_ring(GF(p), PR)

    return zeta_coefficients(hom(f);S=S, verbose=verbose, changef=changef, givefrobmat=givefrobmat, algorithm=algorithm, termorder=termorder, vars_reversed=vars_reversed, fastevaluation=fastevaluation, always_use_bigints=always_use_bigints, use_gpu=use_gpu, use_threads=use_threads, context=context)
end 

function newton_polygon(f; S=[-1], verbose=0, changef=true, algorithm=:default, termorder=:invlex, vars_reversed=false, fastevaluation=true, always_use_bigints=false, use_gpu=false)
    PR = parent(f)
    R = coefficient_ring(PR)
    p = Int64(characteristic(PR))
    zf = zeta_coefficients(f; S=S, verbose=verbose, changef=changef, algorithm=algorithm, termorder=termorder, vars_reversed=vars_reversed, fastevaluation=fastevaluation, always_use_bigints=always_use_bigints, use_gpu=use_gpu)
    
    if zf == false # f isn't smooth
        return false
    end 

    return newton_polygon(p, zf)
end 



"""
    hodgepolygon(basis::Array)

Calculates the hodge polygon of the cohomology module with 
griffiths-dwork basis basis

basis -- an array of "polynomials with pole" as descirbed in PolynomialWithPole.jl
"""
function hodgepolygon(basis::Vector, n)
    #WRONG: n = highestpoleorder(basis)
    hodgenumbers = zeros(Int,n)
    for i in 0:n-1
        h = length(termsoforder(basis,n-i))
        hodgenumbers[i+1] = h
    end

    SlopesPolygon(hodgenumbers)
end

"""
    hodgepolygon(f; basis=nothing, params=default_params())

Calculates the hodge polygon of f

f - the polynomial to get the hodge polygon of
"""
function hodgepolygon(f::RingElem; basis=nothing, params=default_params())
    n = nvars(parent(f)) - 1
    PR = parent(f)
    R = coefficient_ring(parent(f))

    if basis == nothing
        d = total_degree(f)
        S = collect(0:n) #collect(0:d-1)
        cache = controlled_reduction_cache(n,d,S,params)
        basis = compute_monomial_bases(f, params, cache) # basis of cohomology 
    end
    
    Basis = []
    for i in 1:n
        for j in basis[i]
            push!(Basis,[j,i])
        end
    end

    hodgepolygon(Basis,n)
end


"""
    check_smoothness(f)

Using Oscar, check whether f defines a smooth hypersurface.

f is assumed to be homogeneous already. Otherwise Oscar
will throw an error.

note: the name issmooth / is_smooth is already taken by oscar,
and really we're just wrapping that method.
"""
function check_smoothness(f)
    p = characteristic(parent(f))
    nVars = length(gens(parent(f)))

    graded, _ = grade(parent(f))

    R, _ = quo(graded, ideal(graded, [graded(f)]))

    V = proj(R)

    is_smooth(V)
end
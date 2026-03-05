# wrapper function to compute zeta functions for quartic K3 surfaces over F_3s ONLY 
function zeta_function_k3f3(f; S=[-1],verbose=false, precise=true, changef=true, givefrobmat=false, algorithm=:costachunks, termorder=:invlex, vars_reversed=true, fastevaluation=false, always_use_bigints=false, use_gpu=false)
    PR = parent(f)
    R = coefficient_ring(PR)
    p = Int64(3)
    q = p
    n = 3
    d = 4
    @assert total_degree(f) == 4 "f does not define a quartic surface"
    @assert nvars(PR) == 4 "f does not define a hypersurface in P^3"
    @assert characteristic(R) == 3 "f is not defined over FF_3"

    if S == [-1]
        S = collect(0:n)
    end

    params = ZetaFunctionParams(verbose,givefrobmat,algorithm,termorder,vars_reversed,fastevaluation,always_use_bigints,use_gpu)

    cache = controlled_reduction_cache(n,d,S,termorder)


    basis = compute_monomial_bases(f, R, PR, params.termorder) # basis of cohomology 
    Basis = []
    for i in 1:n
        for j in basis[i]
            push!(Basis,[j,i])
        end
    end

    (9 < verbose) && println("Basis of cohomology is $Basis")
    
    hodge_polygon = SlopesPolygon([1, 19, 1])

    if precise 
        (r_m, N_m, M) = ([3, 4, 5], [7, 7, 8], 16)  # guaranteed to give the correct zeta function
    else
        (r_m, N_m, M) = ([1, 2, 2], [4, 4, 3], 6)  # only the first two p-adic digits are guaranteed to be correct 
    end 

    (9 < verbose) && println("We work modulo $p^$M, and compute up to the $N_m-th term of the Frobenius power series")
    (0 < verbose) && println("algorithm precision: $M, series precision: $N_m") 

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
        println("Starting linear algebra problem")
        @time f_changed, f, pseudo_inverse_mat_new = find_Ssmooth_model(f, M, S, params, changef)
    else
        f_changed, f, pseudo_inverse_mat_new = find_Ssmooth_model(f, M, S, params, changef)
    end

    if f == false 
        println("f is not smooth and we're done.")
        return false
    end 

    (9 < verbose) && println("pseudo_inverse_mat is $pseudo_inverse_mat_new")

    # recomputes basis if f is different 
    if f_changed 
        (0 < verbose) && println("New model is $f")
        basis = compute_monomial_bases(f, R, PR, params.termorder) # basis of cohomology 
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
        println("Applying Frobenius to basis...")
        @time FBasis = applyFrobeniusToBasis(Basis,fLift, N_m, p, params)
    else
        FBasis = applyFrobeniusToBasis(Basis,fLift, N_m, p, params)
    end

    l = d * n - n + d - length(S)

    MS = matrix_space(precisionring, nrows(pseudo_inverse_mat_new), ncols(pseudo_inverse_mat_new))
    pseudo_inverse_mat = MS()

    T = computeT(f, basis, M, params)
    (9 < verbose) && println("T matrix is $T")

    for i in 1:nrows(pseudo_inverse_mat_new)
        for j in 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i,j] = ZZ(pseudo_inverse_mat_new[i,j])
        end
    end

    #TODO: check which algorithm we're using
    (2 < verbose) && println("Pseudo inverse matrix:\n$pseudo_inverse_mat")
    Reductions = reducetransform(FBasis, N_m, S, fLift, pseudo_inverse_mat, p,  params, cache) 
    (1 < verbose) && println(Reductions)

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
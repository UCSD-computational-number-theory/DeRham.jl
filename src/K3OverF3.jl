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

    #S = SmallestSubsetSmooth.smallest_subset_s_smooth(fLift,n)

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

    #=
    BasisTLift = []
    for i in basis
        temp = []
        for j in i
            push!(temp, liftCoefficients(precisionring,precisionringpoly,j))
        end
        push!(BasisTLift,temp)
    end

    Basis = []
    
    for i in 1:n
        for j in 1:length(Basis[i])
            #push!(Basis,[j,i])
            Basis[i,j] = liftCoefficients(precisionring, precisionringpoly, Basis[i,j])
        end
    end
    =#
    
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
    #println(FBasis)
    #for e in FBasis
    #    println(length(e))
    #    for t in e
    #        print("   " * "$(length(terms(t[1]))): ")
    #        println(total_degree(t[1]))
    #    end
    #end
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
    #=
    pseudo_inverse_mat = zeros(Int, nrows(pseudo_inverse_mat_new),ncols(pseudo_inverse_mat_new))
    for i in 1:nrows(pseudo_inverse_mat_new)
        for j in 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i,j] = ZZ(pseudo_inverse_mat_new[i,j])
        end
    end
    printMat(pseudo_inverse_mat)
    =#
    #pseudoInverseMat = zeros(PrecisionRing, nrows(pseudoInverseMatTemp), ncols(pseudoInverseMatTemp))

    #PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    #fLift = liftCoefficients(ZZ,PRZZ,f)
    #controlledMatrixZZ = compute_controlled_matrix(fLift, d * n - n + d - length(S), S, ZZ, PRZZ)
    #pseudoInverseMatModP = matrix(ZZ, [lift(ZZ,x) for x in Array(pseudoInverseMatTemp)])
    #pseudo_inverse_mat_new = henselLift(p,M,controlledMatrixZZ, pseudoInverseMatModP)
    
    #for i in 1:nrows(pseudoInverseMat)
    #    for j in 1:ncols(pseudoInverseMat)
    #        pseudoInverseMat[i,j] = PrecisionRing(lift(ZZ, pseudoInverseMatTemp[i,j]))
    #    end
    #end
    #TODO: check which algorithm we're using
    (2 < verbose) && println("Pseudo inverse matrix:\n$pseudo_inverse_mat")
    Reductions = reducetransform(FBasis, N_m, S, fLift, pseudo_inverse_mat, p,  params, cache) 
    (1 < verbose) && println(Reductions)
    #return Reductions
    #if (1 < verbose)
    #    for i in 1:length(Basis)
    #        basis_elt = Basis[i]
    #        after_reduction = Reductions[i]
    #        println("Basis element $basis_elt becomes $after_reduction after controlled reduction")
    #        println()
    #    end 
    #end 
    ev = cache[n*d - n - 1] 
    #ev = gen_exp_vec(n+1,n*d-n-1,termorder)
    (9 < verbose) && println(convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev))
    FM = compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, params, cache)
    # display(FM)
    (9 < verbose) && println("The Frobenius matrix is $FM")

    #reductions_verbose = convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev)

    #(9 < verbose) && println("convert_p_to_m is $reductions_verbose")

    #FM = compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis)

   # (9 < verbose) && println("The Frobenius matrix is $FM")

    if givefrobmat
        (FM,LPolynomial(FM,n,q,hodge_polygon,r_m, verbose))
    else
        LPolynomial(FM,n,q,hodge_polygon,r_m, verbose)
    end
end
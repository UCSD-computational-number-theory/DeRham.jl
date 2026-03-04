
function reducechain_varbyvar(u,g,m,S,f,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))

    J = copy(u)

    modm = Integer(modulus(base_ring(parent(f))))

    gMat = g
    mins = copy(J)
    tempv = copy(J)

    (4 < params.verbose) && println("Starting: J = $J")
    #=
    (5 < params.verbose) && begin
        CUDA.@allowscalar g_poly = vector_to_polynomial(g,n,d*n-n,PR,params.termorder)
        if params.always_use_bigints || params.use_gpu
            println("Starting: g = $((gMat)) = $g_poly")
        else    
            println("Starting: g = $(Int.(gMat)) = $g_poly")
        end
    end
    =#

    l = 1
    for i in eachindex(J)
        if J[n+1-i+1] > 0
            l = n+1-i+1
            break
        end
    end

    highpole = true
    if m == n
        highpole == false
    end
    while highpole && J[l] > 0
        V = varbyvar_chooseV(J,d)
        
        (4 < params.verbose) && print("Chose V = $V; ")
        (6 < params.verbose) && begin
            # the way that chooseV works right now,
            # the following if statement should never hit.
            # for i in 1:length(V)
            #     if params.vars_reversed && V[i] == 0 && J[i] ≠ 0 && (n+1-i) ∈ S
            #         print("Illegal choice of V!")
            #         println("J = $J, S = $S")
            #     end
            # end
        end

        K = 0
        while true
            @. tempv = mins - V
            isLessThanZero = false
            for j in tempv
                if j < 0
                    isLessThanZero = true
                    break
                end
            end
            if isLessThanZero == true
                break
            end
            if m - K == n
                break
            end
            @. mins = tempv
            K = K+1
        end

        (4 < params.verbose) && println("Getting Ruv matrices; ")
        matrices = context.Ruvs[V]
        (4 < params.verbose) && println("Computing A and B; ")
        if params.use_gpu && !(ZZ(2)^25 < modm < ZZ(2)^106) # so not Karatsuba
            eval_to_linear_gpu!(context.B,context.A,context.temp,matrices,mins,V)
        elseif params.use_gpu && d == 3 && (n == 4 || n == 5)# karatsuba
            eval_to_linear_gpu_karatsuba!(context.B,context.A,context.temp,matrices,mins,V)
        else
            eval_to_linear!(context.B,context.A,context.temp,matrices,mins,V)
        end
        (4 < params.verbose) && println("Starting the reduction steps; ")
        gMat = finitediff_prodeval_linear!(context.B,context.A,0,K-1,gMat,context.temp,context.g_temp)
        @. J = J - K*V
        m = m - K
        if m == n
            highpole = false
        end

        (4 < params.verbose) && print("After $(lpad(K,4,' ')) steps,")
        (4 < params.verbose) && println("J = $J")
        #=
        if (5 < params.verbose) 
            CUDA.@allowscalar g = vector_to_polynomial(gMat,n,d*n-n,PR,params.termorder)
            if params.always_use_bigints || params.use_gpu
                println("g = $((gMat)) = $g")
            elseif params.fastevaluation
                println("g = $(Int.(gMat)) = $g")
            else 
                println("g = $(gMat) = $g")
            end
        end
        =#

    end

    return ((J, gMat), m)
end

function reducepoly_varbyvar(pol,S,f,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    result = PR()

    i = pol
    highpoleorder = i[length(i)][2]
    terms = []
    while highpoleorder >= p
        append!(terms,termsoforder(pol,highpoleorder))
        highpoleorder = highpoleorder - p
    end

    vectype = typeof(context.g)
    allcostadata = Vector{Tuple{Tuple{Vector{Int},vectype},Int}}()
    for term in terms
        #g = my_copy(context.g)
        term_costadata = costadata_of_initial_term!(term,my_copy(context.g),n,d,p,S,cache,params)
        append!(allcostadata,[((tweak(term_costadata[1],n*d-n),term_costadata[2]),term[2])])
    end

    notallred = true
    # TODO: prime the jitter here?

    while notallred
        for i in eachindex(allcostadata)
            allcostadata[i] = reducechain_varbyvar(allcostadata[i][1]...,allcostadata[i][2],S,f,p,context,cache,params)
        end
        if params.use_gpu == true
            remove_duplicates_gpu!(allcostadata)
        else
            remove_duplicates!(allcostadata)
        end
        for i in 1:length(allcostadata)
            if allcostadata[i][2] > n
                break
            elseif i == length(allcostadata)
                notallred = false
            end
        end
    end
    result = PR()
    for i in eachindex(allcostadata)
        (reduced_poly,m) = poly_of_end_costadata(allcostadata[i][1],PR,p,d,n,params)
        result += reduced_poly
    end

    vars = gens(PR)
    XS = prod(PR(vars[i+1]) for i in S; init = PR(1))
    [[div(result,XS), n]]
end

function reducetransform_varbyvar(FT,N_m,S,f,pseudoInverseMat,p,cache,params,context)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g_length = binomial(d*n,d*n-n)

    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)
    m = Integer(modulus(base_ring(MS1)))

    #Ruvs = Dict{Vector{Int64}, Vector{typeof(MS1())}}()

    #explookup = Dict{Vector{Int64}, Int64}()
    #ev1 = gen_exp_vec(n+1,n*d-n,params.termorder)
    #for i in 1:length(ev1)
    #    get!(explookup,ev1[i],i)
    #end

    if (3 < params.verbose)
        computeRuv = V -> begin
            println("Computing Ruv for V = $V for the first time.")
            @time computeRuvS(V,S,f,pseudoInverseMat,cache,params)
        end
    else
        computeRuv = V -> begin
            computeRuvS(V,S,f,pseudoInverseMat,cache,params)
        end
    end

    result = similar(FT)

    if context == nothing
        #TODO: right now, it usually seems better to do lazy computations,
        #  since not all of the Ruv are used. However, I know that for some
        #  classes of examples, they are all pretty much always used. For 
        #  such examples, it's better to use an EagerPEP and do threads.
        lazy_Ruv = true#length(S) < d || d < n

        if (0 < params.verbose)
            println("Creating the Ruv PEP object...")
            #CUDA.@time Ruv = select_Ruv_PEP(params,computeRuv,computeRuv_gpu,lazy_Ruv,MS1,cache,d)
            @time Ruv = select_Ruv_PEP(n,d,S,params,computeRuv,lazy_Ruv,MS1,cache)
        else
            Ruv = select_Ruv_PEP(n,d,S,params,computeRuv,lazy_Ruv,MS1,cache)
        end

        if params.use_threads
        
            context_tlv = OhMyThreads.TaskLocalValue{default_context_type(MS1,params)}(
                () -> default_context(MS1,Ruv,params)
            )
            Threads.@threads for i in 1:length(FT) 
                local context = context_tlv[]

                pol = FT[i]
                if (0 < params.verbose)
                    println("Reducing vector $i")
                    @time reduction = reducepoly_varbyvar(pol,S,f,p,context,cache,params)
                else
                    reduction = reducepoly_varbyvar(pol,S,f,p,context,cache,params)
                end
                result[i] = reduction

                #println("cache info: $(cache_info(Ruv.Ucomponent))")
                #i == 5 && error("stopping after vector $i for testing purposes")
                
                #push!(result, reduction)
            end
        else 

            context = default_context(MS1,Ruv,params)
            for i in 1:length(FT) #pol in FT
                

                pol = FT[i]
                if (0 < params.verbose)
                    println("Reducing vector $i")
                    @time reduction = reducepoly_varbyvar(pol,S,f,p,context,cache,params)
                else
                    reduction = reducepoly_varbyvar(pol,S,f,p,context,cache,params)
                end
                result[i] = reduction

                #println("cache info: $(cache_info(Ruv.Ucomponent))")
                #i == 5 && error("stopping after vector $i for testing purposes")
                    
                #push!(result, reduction)
            end
        end
    else
        if (ZZ(2)^25 < m < ZZ(2)^106)
            compute_gpu = V -> KaratsubaMat.(computeRuv(V))
        else 
            compute_gpu = V -> cuMod.(computeRuv(V))
        end 
    
        compute_float = V -> float_entries.(computeRuv(V))

        if 3 < n && d == 3 && S == [n] && params.use_gpu 
            context.Ruvs.backing.compute = compute_float
        elseif params.use_gpu
            context.Ruvs.compute = compute_gpu
        else
            context.Ruvs.compute = computeRuv
        end
        for i in 1:length(FT) #pol in FT
            

            pol = FT[i]
            if (0 < params.verbose)
                println("Reducing vector $i")
                @time reduction = reducepoly_varbyvar(pol,S,f,p,context,cache,params)
            else
                reduction = reducepoly_varbyvar(pol,S,f,p,context,cache,params)
            end
            result[i] = reduction

            #println("cache info: $(cache_info(Ruv.Ucomponent))")
            #i == 5 && error("stopping after vector $i for testing purposes")
                
            #push!(result, reduction)
        end
    end
    

    #(0 < params.verbose && Ruv isa CachePEP) && begin
    #    println("Ruv cache info: $(cache_info(Ruv.Ucomponent))")
    #end
    (0 < params.verbose) && begin
        println("Created $(length(allpoints(context.Ruvs))) of $(length(cache[d])) possible V")
    end
    (1 < params.verbose) && begin
        println("V that were created: \n$(allpoints(context.Ruvs))")
    end

    return result
end
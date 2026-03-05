
# what in here could possibly be causing a lock conflict?
"""
Iteratively execute reduction chunks in the navie strategy until
the vector g is reduced to pole order n
"""
function reducechain_depthfirst(u,g,m,S,f,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    (9 < params.verbose) && println("u = $u") 
    
    J = tweak(u,n*d-n)

    gMat = context.g
    mins = similar(J)
    tempv = similar(J)
    (4 < params.verbose) && println("Starting: J = $J")
    (5 < params.verbose) && begin
        g_poly = vector_to_polynomial(g,n,d*n-n,PR,params.termorder)
        if params.always_use_bigints || params.use_gpu
            println("Starting: g = $((gMat)) = $g_poly")
        else    
            println("Starting: g = $(Int.(gMat)) = $g_poly")
        end
    end

    firsttime = true


    while m > n
        V = rev_chooseV(J,d,S)

        (4 < params.verbose) && print("Chose V = $V; ")
        (6 < params.verbose) && begin
        end
        @. mins = J
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
        matrices = context.Ruvs[V]

        (6 < params.verbose && V == [0,0,0,3] && firsttime) && begin println(matrices[5][:,25]); firsttime=false; error() end
        
        eval_to_linear!(context.B,context.A,context.temp,matrices,mins,V)

        i = 1
        if params.fastevaluation == false
            params.verbose == 11 && println("gMat before is $gMat")
            while i <= K
                gMat = (context.A+context.B*(K-i))*gMat
                i = i+1
                if params.verbose == 11
                    A = context.A
                    B = context.B
                    g = vector_to_polynomial(gMat,n,d*n-n,PR,params.termorder)
                    println("gMat after $i is $gMat = $g")
                end
            end
        else
            gMat = finitediff_prodeval_linear!(context.B,context.A,0,K-1,gMat,context.temp,context.g_temp)
        end
        @. J = J - K*V
        m = m - K
        (4 < params.verbose) && print("After $(lpad(K,4,' ')) steps,")
        (4 < params.verbose) && println("J = $J")
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
        
    end
    return (J, gMat)
end

function reducepoly_depthfirst(pol,S,f,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    result = PR()
    for term in pol
        terms = termsoforder(pol,term[2])
        for t in terms
            (u,_) = costadata_of_initial_term!(t,context.g,n,d,p,S,cache,params)
            reduced = reducechain_depthfirst(u,context.g,t[2],S,f,p,context,cache,params)
            (reduced_poly,m) = poly_of_end_costadata(reduced,PR,p,d,n,params)
            @assert m == n "Controlled reduction outputted a bad pole order"
            result += reduced_poly
        end
    end

    vars = gens(PR)
    XS =  prod(PR(vars[i+1]) for i in S; init = PR(1))
    [[div(result,XS), n]]
end

function reducetransform_depthfirst(FT,N_m,S,f,pseudoInverseMat,p,cache,params,context)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g_length = binomial(d*n,d*n-n)

    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)
    m = Integer(modulus(base_ring(MS1)))

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
        lazy_Ruv = length(S) < d || d < n

        if (0 < params.verbose)
            println("Creating the Ruv PEP object...")
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
                    println("Reducing vector $i in thread $(Threads.threadid())")
                    @time reduction = reducepoly_depthfirst(pol,S,f,p,context,cache,params)
                else
                    reduction = reducepoly_depthfirst(pol,S,f,p,context,cache,params)
                end
                result[i] = reduction

            end
        else 
            context = default_context(MS1,Ruv,params)
            for i in 1:length(FT) #pol in FT
                
                pol = FT[i]
                if (0 < params.verbose)
                    println("Reducing vector $i")
                    @time reduction = reducepoly_depthfirst(pol,S,f,p,context,cache,params)
                else
                    reduction = reducepoly_depthfirst(pol,S,f,p,context,cache,params)
                end
                result[i] = reduction

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
                println("Reducing vector $i in thread $(Threads.threadid())")
                @time reduction = reducepoly_depthfirst(pol,S,f,p,context,cache,params)
            else
                reduction = reducepoly_depthfirst(pol,S,f,p,context,cache,params)
            end
            result[i] = reduction
        end
    end

    (0 < params.verbose) && begin
        println("Created $(length(allpoints(Ruv))) of $(length(cache[d])) possible V")
    end
    (1 < params.verbose) && begin
        println("V that were created: \n$(allpoints(Ruv))")
    end

    return result
end
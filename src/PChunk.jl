
"""
    reducechain_costachunks(u,g,n,d,p,m,S,f,pseudoInverseMat,R,PR)

takes single monomial in frobenius and reduces to pole order n, currently only does one chunk of reduction


if the reduction hits the end, returns u as the "true" value, otherwise returns it in Costa's format
(i.e. entries will be multiplies of p in Costa's format)
"""
function reducechain_pchunk(u,g,m,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,g_temp,params)
    verbose = params.verbose
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    

    
    I = u
   
    (4 < verbose) && println("Expanded I: $I")

    gMat = g
    #(4 < verbose) && println("This is I: $I_edgar")
    J = copy(I)

    #TODO?
    # if params.vars_reversed == false
        V = rev_chooseV(Array{Int}(divexact.(I,p)),d, S)
    # else
    #     V = chooseV(Array{Int}(divexact.(I,p)),d, S)
    # end
    (4 < verbose) && println("LOOK! I=$I, V = $V")


    # if params.vars_reversed == true
    #      gVec = I .- rev_tweak(I,n*d-n)
    # else
         gVec = I .- tweak(I,n*d-n)
    # end
    #ev = gen_exp_vec(n+1,n*d-n,termorder)

    @. I = I - gVec
    if m - n < p
        nend = m - n
    else
        nend = p
    end

    matrices = Ruv[V]#computeRuvS(V,S,f,pseudoInverseMat,Ruvs,cache,params)

    #TODO: the following was changed to use reverse! so 
    #    it doesn't allocate as much, but I realized that
    #    this isn't our bottleneck. If it ever does
    #    become the bottleneck, fix the rest of this method
    #    so it doesn't allocate.
    U = I .- (nend-(d*n-n))*V
    #reverse!(U)

    #reverse!(V)
    B,A = eval_to_linear!(B,A,temp,matrices,U,V)
    #reverse!(V) # put V back to normal

    i = 1

    
    (4 < verbose) && println("Before reduction chunk, I is $I")
    if params.fastevaluation && 1 ≤ nend-(d*n-n)
      gMat = finitediff_prodeval_linear!(B,A,0,nend-(d*n-n)-1,gMat,temp,g_temp)
      i = nend-(d*n-n) + 1
    else
      while i <= (nend-(d*n-n))
        my_mul!(temp,B,nend-(d*n-n)-i)
        my_add!(temp,temp,A)
        gMat = temp*gMat

        #gMat = (A+B*(nend-(d*n-n)-i))*gMat

        #(9 < verbose) && println("After step $i: $(convert.(Int,gMat))")

        i = i+1
        #println(gMat)
      end
    end
    # TODO: test how much of a difference the fast evaluation actually makes
    # i > 1 iff the while loop above is executed at least once 
    if i > 1 # TODO: this will have a problem with fastevaluation
        # UPDATE: I think the problem with fastevaluation is fixed... right?
        @. I = I - (nend-(d*n-n))*V
    end
    (4 < verbose) && println("After steps 1-$i, I is $I")
    i = i-1
    while i <= nend-1
        # if params.vars_reversed == true
        #     y = rev_tweak(J - i*V,d*n-n) .- rev_tweak(J - (i+1)*V,d*n-n)
        # else
            y = tweak(J - i*V,d*n-n) .- tweak(J - (i+1)*V,d*n-n)
        # end
        (4 < verbose) && println("Getting y direction reduction matrix for V = $(y)") 
        
        # there's some sort of parity issue between our code and Costa's
        #A,B = computeRPoly_LAOneVar(y,rev_tweak(J - (i+1)*V,d*n-n) - y,S,n,d,f,pseudoInverseMat,R,PR,termorder)
        
        matrices1 = Ruv[y]#computeRuvS(y,S,f,pseudoInverseMat,Ruvs,cache,params)
        #println(matrices1)

        #if params.vars_reversed == true
        #    #B,A = eval_to_linear!(B,A,temp,matrices1,reverse(rev_tweak(J - (i+1)*V,d*n-n) - y),reverse(y))
        #    B,A = eval_to_linear!(B,A,temp,matrices1,rev_tweak(J - (i+1)*V,d*n-n) - y,y)
        #else
            #B,A = eval_to_linear!(B,A,temp,matrices1,reverse(tweak(J - (i+1)*V,d*n-n) - y),reverse(y))
            B,A = eval_to_linear!(B,A,temp,matrices1,tweak(J - (i+1)*V,d*n-n) - y,y)
        # end
        
        my_add!(temp,A,B)
        my_matvecmul!(g_temp,temp,gMat)
        my_copy!(gMat,g_temp)
        #gMat = temp*gMat

        #gMat = (A+B)*gMat

        (4 < verbose) && println("After step $(i+1): $(gMat))")
        

        i = i+1
        @. I = I - y
        (4 < verbose) && println("After step $(i+1), I is $I")
    end
    
    if nend == p
        newI = J .- p*V
        #@assert undo_rev_tweak(I,p) == newI

        return (newI, gMat)
    else
        return (I,gMat) # gives the "true" u
    end
end

"""
    reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruvs,termorder)

Implements Costa's algorithm for controlled reduction,
sweeping down the terms of the series expansion by the pole order.
"""
function reducepoly_pchunk(pol,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,params)
    #p = Int64(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    #(9 < verbose) && println(pol)
    g_length = binomial(d*n,d*n-n)

    i = pol
    highpoleorder = i[length(i)][2]

    # this is the omega from section 1.5.5 of Costa's thesis.
    ω = [] # this will be an array of costa data

    poleorder = highpoleorder
    while n < poleorder
        (9 < params.verbose) && println("pole order is $poleorder")
        # this is an array of polynomials
        ωₑ = termsoforder(pol,poleorder)

        #(9 < verbose) && println("ωₑ: $ωₑ")

        for term in ωₑ
            #(9 < verbose) && println("term: $term")
            g = zeros(UInt,g_length) 
            term_costadata = costadata_of_initial_term!(term,g,n,d,p,S,cache,params)
            #(9 < verbose) && println("term, in Costa's format: $term_costadata")
            #ω = ω + ωₑ
            incorporate_initial_term!(ω,term_costadata)
        end

        #(9 < verbose) && println("ω: $ω")
        #ω = reducepoly_LA(ω,n,d,p,S,f,pseudoInverseMat,R,PR)
        for i in eachindex(ω)
            #ω[i] = reducechain...
            #(9 < verbose) && println("u is type $(typeof(ω[i][1]))")
            g_temp = similar(ω[i][2])
            ω[i] = reducechain_pchunk(ω[i]...,poleorder,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,g_temp,params)
        end

        poleorder = poleorder - p
    end

    #println("ω: $ω")
    #(9 < verbose) && println(poly_of_end_costadatas(ω,PR,p,d,n,S,termorder))

    #println(gen_exp_vec(n,n*d-n-1,termorder))
           

    return poly_of_end_costadatas(ω,PR,p,d,n,S,params)
end

"""
    reducetransform_costachunks(FT,N_m,S,f,pseudoInverseMat,p,cache,params)

trying to emulate Costa's controlled reduction, changes the order that polynomials are reduced, starts from highest pole order and accumulates the lower order poles as reduction proceeds

N_m - the precision
"""
function reducetransform_pchunk(FT,N_m,S,f,pseudoInverseMat,p,cache,params)

    d = total_degree(f)
    n = nvars(parent(f)) - 1
    MS1 = matrix_space(coefficient_ring(parent(f)), binomial(d*n,d*n-n), binomial(d*n,d*n-n))
    A = MS1()
    B = MS1()
    temp = MS1()
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
    Ruv = LazyPEP{typeof(MS1())}(computeRuv)

    result = []
    i = 1
    for pol in FT
        (0 < params.verbose) && println("Reducing vector $i")
        i += 1
        if (0 < params.verbose)
            @time reduction = reducepoly_pchunk(pol,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,params)
        else
            reduction = reducepoly_pchunk(pol,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,params)
        end

        push!(result, reduction)
    end

    return result
end
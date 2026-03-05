
"""
    reducechain_pchunks(u,g,n,d,p,m,S,f,pseudoInverseMat,R,PR)

takes single monomial in frobenius and reduces to pole order n, currently only does one chunk of reduction


if the reduction hits the end, returns u as the "true" value, otherwise returns it in Costa's format
(i.e. entries will be multiplies of p in Costa's format)

fields
------
u -  vector of ints
g - vector
m - pole order
S - vector of ints
f - polynomial
pseudoInverseMat - output of pseudo_inverse_controlled_lifted
p - prime number
Ruv - output of computeRuvS
cache - the GradedExpCache used for this controlled reduction
A - matrix
B - matrix
temp - preallocated storage matrix
g_temp - preallocated storage vector
params - the ControlledReductionParamaters
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
    J = copy(I)

    #TODO?
    V = rev_chooseV(Array{Int}(divexact.(I,p)),d, S)
    (4 < verbose) && println("LOOK! I=$I, V = $V")

    gVec = I .- tweak(I,n*d-n)

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

    B,A = eval_to_linear!(B,A,temp,matrices,U,V)

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

        i = i+1
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
        y = tweak(J - i*V,d*n-n) .- tweak(J - (i+1)*V,d*n-n)
        (4 < verbose) && println("Getting y direction reduction matrix for V = $(y)") 
        
        matrices1 = Ruv[y]#computeRuvS(y,S,f,pseudoInverseMat,Ruvs,cache,params)

        B,A = eval_to_linear!(B,A,temp,matrices1,tweak(J - (i+1)*V,d*n-n) - y,y)
        
        my_add!(temp,A,B)
        my_matvecmul!(g_temp,temp,gMat)
        my_copy!(gMat,g_temp)

        (4 < verbose) && println("After step $(i+1): $(gMat))")
        

        i = i+1
        @. I = I - y
        (4 < verbose) && println("After step $(i+1), I is $I")
    end
    
    if nend == p
        newI = J .- p*V

        return (newI, gMat)
    else
        return (I,gMat) # gives the "true" u
    end
end

"""
    reducepoly_pchunk(pol,S,f,pseudoInverseMat,p,Ruvs,termorder)

Implements Costa's algorithm for controlled reduction,
sweeping down the terms of the series expansion by the pole order.

fields
------
pol - polynomial
S - vector of ints
f - polynomial
pseudoInverseMat - output of pseudo_inverse_controlled_lifted
p - prime number
Ruv - output of computeRuvS
cache - the GradedExpCache used for this controlled reduction
A - matrix
B - matrix
temp - preallocated storage matrix
params - the ControlledReductionParamaters
"""
function reducepoly_pchunk(pol,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
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

        for term in ωₑ
            g = zeros(UInt,g_length) 
            term_costadata = costadata_of_initial_term!(term,g,n,d,p,S,cache,params)
            incorporate_initial_term!(ω,term_costadata)
        end

        for i in eachindex(ω)
            g_temp = similar(ω[i][2])
            ω[i] = reducechain_pchunk(ω[i]...,poleorder,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,g_temp,params)
        end

        poleorder = poleorder - p
    end    

    return poly_of_end_costadatas(ω,PR,p,d,n,S,params)
end

"""
    reducetransform_pchunk(FT,N_m,S,f,pseudoInverseMat,p,cache,params)

trying to emulate Costa's controlled reduction, changes the order that polynomials are reduced, starts from highest pole order and accumulates the lower order poles as reduction proceeds

fields
------
FT - applyFrobeniusToBasis
N_m - the precision
S - vector of ints
f - polynomial
pseudoInverseMat - output of pseudo_inverse_controlled_lifted
p - prime number
cache - the GradedExpCache used for this controlled reduction
params - the ControlledReductionParamaters
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
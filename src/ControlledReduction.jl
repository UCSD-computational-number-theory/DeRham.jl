
"""
    ControlledReductionContext{MatrixType,VectorType}

A controlled reduction context is used for in-place evaluation
chunks of controlled reduction.
It contains all of the major allocations necessary, so you can
create one context at the beginning (or one per thread)
and then not have to worry about allocating anything else.

MatrixType - the type of A and B, for example zzModMatrix
VectorType - the type of g, for exampel Vector{UInt64}

fields
------
Ruvs - an AbstractPEP which containes the R_uv necessary for the reduction.
  This is intended to be shared across all contexts (i.e. all threads)
  for a single example.
A - the matrix which will store A in Ay + B
B - the matrix which will store B in Ay + B
temp - temporary matrix used for intermediate calculations
g - the vector g which contains the vector which is being reduced
g_temp - temporary vector used for intermediate calculations

Current examples that are in use:

ControlledReductionContext{zzModMatrix,Vector{UInt64}}
ControlledReductionContext{ZZModMatrix,Vector{ZZRingElem}}
ControlledReductionContext{CuModMatrix,CuModVector}

"""
struct ControlledReductionContext{MatrixType,VectorType}
    Ruvs::AbstractPEP{MatrixType}
    A::MatrixType
    B::MatrixType
    temp::MatrixType
    g::VectorType
    g_temp::VectorType
end


#"""
#    reduce_LA(U,V,S,f,pseudoInverseMat,g,PR,termorder)
#
#applies reduction formula from Prop 1.15 in Costa's thesis to 
#basis elements of Homog(dn-d), returns them as polynomials
#will only work with vars_reversed=true
#"""
#function reduce_LA(U,V,S,f,pseudoInverseMat,g,PR,termorder)
#    R = coefficient_ring(PR)
#    Vars = gens(PR)
#    n = nvars(parent(f)) - 1
#    d = total_degree(f)
#    SC = []
#    B = MPolyBuildCtx(PR)
#    push_term!(B, R(1), V)
#    XV = finish(B)
#    for i in 0:n
#        if i in S
#        else
#            push!(SC,i)
#        end
#    end
#    # get gi's using pseudoinverse
#    XS =  prod(PR(Vars[i+1]) for i in S; init = PR(1))
#    gVec = convert_p_to_m([div(XV*(g[1]),XS)],gen_exp_vec(n+1,n*d-n+d-length(S),termorder))
#    MS = matrix_space(parent(gVec[1]), nrows(pseudoInverseMat),1)
#    gJS = MS()
#    gJS = pseudoInverseMat*transpose(gVec)
#    gc = []
#    for i in 1:(n+1)
#        push!(gc, convert_m_to_p(transpose(gJS[Int((i-1)*(length(gJS)/(n+1))+1):Int(i*(length(gJS)/(n+1))),:]),gen_exp_vec(n+1,n*d-n-d+1,termorder),R,PR)[1])
#    end
#    gc = reverse(gc)
#    gcpartials = [ derivative(gc[i], i) for i in 1:(n+1) ]
#    
#    reverse!(gcpartials) # TODO: make this an option, this is the way it is in Costa's code, 
#
#    #return [sum(PR(U[i+1])*XS*gc[i+1] + div(XS,Vars[i+1])*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]
#    return [sum(PR(U[i+1])*div(XS,Vars[i+1])*gc[i+1] + XS*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0))]
#c
#end

"""
    chooseV(I, d, S)
Choose the direction of reduction V following Edgar's method

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f
* "S" -- list 
"""
function chooseV(I, d, S)
    v = zeros(Int,length(I))
    n = length(I)

    sum = 0
    
    #=
    for i in 1:n 
        if I[i] > 0
            v[i] = v[i] + 1
            sum = sum + 1
            if sum >= d 
                break
            end
        end 
    end 
    =#
    
    
    for i in S 
        if I[n-i] > 0
            v[n-i] = v[n-i] + 1
            sum = sum + 1
            if sum >= d 
                break
            end
        end 
    end 
    

    while sum < d
        for i in 0:(n-1) 
            if I[n-i] > v[n-i]
                v[n-i] = v[n-i] + 1
                sum = sum + 1
                break
            end 
        end 
    end 
    return v
end 

function varbyvar_chooseV(I, d)
    v = zeros(Int,length(I))
    n = length(I)
    J = copy(I)
    i = 0
    while i < d
        for j in eachindex(J)
            if J[n-j+1] > 0
                J[n-j+1] = J[n-j+1] - 1
                v[n-j+1] = v[n-j+1] + 1
                i = i + 1
                break
            elseif j == n
                throw(error)
            end
        end
    end
    return v
end


"""
    skim_chooseV(I, d)

choose direction of reduction V by scanning from right to left
Note: this differs from Edgar's method, see chooseV

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f 
"""
function skim_chooseV(I, d)
    V = zeros(Int,length(I))
    i = 0
    #s = 1
    s = length(I)
    foundNonZero = false
    while i < d
        #=
        if s > length(I) && foundNonZero == false
            return V
        elseif s > length(I)
            s = 1
            foundNonZero = false
        end
        if (I - V)[s] > 0
            V[s] = V[s] + 1
            i = i + 1
            foundNonZero = true
        end
        =#
        #FIXME reversed to match Costa's
        if s == 0 && foundNonZero == false
            return V
        elseif s == 0
            s = length(I)
            foundNonZero = false
        end
        if (I - V)[s] > 0
            V[s] = V[s] + 1
            i = i + 1
            foundNonZero = true
        end
        #s = s + 1
        s = s-1
    end
    return V
end

"""
    rev_chooseV(I, d)

choose direction of reduction in the same way as Costa's code

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f 
"""
function rev_chooseV(I, d)
    reverse!(I)

    V = zeros(Int,length(I))
    i = 0
    s = 1
    foundNonZero = false
    while i < d
        if s > length(I) && foundNonZero == false
            return V
        elseif s > length(I)
            s = 1
            foundNonZero = false
        end
        if (I - V)[s] > 0
            V[s] = V[s] + 1
            i = i + 1
            foundNonZero = true
        end
        s = s + 1
    end

    reverse!(V)

    return V
end

"""
    tweak(I,m)

If for a vectors of ints, I, we let |I| = sum(I[i]). This function returns I after removing an integer vector, J, with |J| = m
from it

Removes from the "front" of the vector"

INPUTS
* "I" -- vector of nonnegative integers 
* "m" -- nonnegative integer
"""
function tweak(J,m)
    count = 0
    o = m
    I = copy(J)
    while m > 0
        for i in axes(I,1)
            if (I[i] > 0)
                I[i] = I[i] - 1
                m = m - 1
                break
            end
        end
        if count > length(I)*o
            return I
        end
        count = count+1
    end
    return I
end


"""
    rev_tweak(J,m)

If for a vectors of ints, I, we let |I| = sum(I[i]). This function returns I after removing an integer vector, J, with |J| = m
from it

unlike tweak, this, removes from the "back" of the vector instead of the front.

INPUTS
* "I" -- vector of nonnegative integers 
* "m" -- nonnegative integer
"""
function rev_tweak(J,m)
    count = 0
    o = m
    I = copy(J)
    while m > 0
        for i in reverse(axes(I,1))
            if (I[i] > 0)
                I[i] = I[i] - 1
                m = m - 1
                break
            end
        end
        if count > length(I)*o
            return I
        end
        count = count+1
    end
    return I
end

"""
    undo_rev_tweak(I,p)

gives the unique vector J containing only
multiples of p such that tweak(J, |J-I|) = I

We're assuming that |J-I| < p here.

"""
function undo_rev_tweak(I,p)

    # is this correct for large d and n?
    J = copy(I)

    k = 1
    #accumulator = 0
    while k ≤ length(J)
        while (J[k] % p) != 0
            J[k] = J[k] + 1
            #accumulator += 1
        end
        k = k + 1
    end

    # @assert accumulator == m

    J
end

"""
    reducechain_costachunks(u,g,n,d,p,m,S,f,pseudoInverseMat,R,PR)

takes single monomial in frobenius and reduces to pole order n, currently only does one chunk of reduction


if the reduction hits the end, returns u as the "true" value, otherwise returns it in Costa's format
(i.e. entries will be multiplies of p in Costa's format)
"""
function reducechain_costachunks(u,g,m,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,g_temp,params)
    verbose = params.verbose
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    

    
    I = u
    #TODO?
    if params.vars_reversed == false
        reverse!(I) # parity issue due to Costa's code being reverse from ours
    end
    (4 < verbose) && println("Expanded I: $I")

    gMat = g
    #(4 < verbose) && println("This is I: $I_edgar")
    J = copy(I)

    #TODO?
    if params.vars_reversed == false
        V = rev_chooseV(Array{Int}(divexact.(I,p)),d)
    else
        V = chooseV(Array{Int}(divexact.(I,p)),d, S)
    end
    (4 < verbose) && println("LOOK! I=$I, V = $V")


    if params.vars_reversed == true
         gVec = I .- rev_tweak(I,n*d-n)
    else
         gVec = I .- tweak(I,n*d-n)
    end
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
      gMat = finitediff_prodeval_linear!(B,A,0,nend-(d*n-n)-1,gMat,temp,ui)
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
        if params.vars_reversed == true
            y = rev_tweak(J - i*V,d*n-n) .- rev_tweak(J - (i+1)*V,d*n-n)
        else
            y = tweak(J - i*V,d*n-n) .- tweak(J - (i+1)*V,d*n-n)
        end
        (4 < verbose) && println("Getting y direction reduction matrix for V = $(y)") 
        
        # there's some sort of parity issue between our code and Costa's
        #A,B = computeRPoly_LAOneVar(y,rev_tweak(J - (i+1)*V,d*n-n) - y,S,n,d,f,pseudoInverseMat,R,PR,termorder)
        
        matrices1 = Ruv[y]#computeRuvS(y,S,f,pseudoInverseMat,Ruvs,cache,params)
        #println(matrices1)

        if params.vars_reversed == true
            #B,A = eval_to_linear!(B,A,temp,matrices1,reverse(rev_tweak(J - (i+1)*V,d*n-n) - y),reverse(y))
            B,A = eval_to_linear!(B,A,temp,matrices1,rev_tweak(J - (i+1)*V,d*n-n) - y,y)
        else
            #B,A = eval_to_linear!(B,A,temp,matrices1,reverse(tweak(J - (i+1)*V,d*n-n) - y),reverse(y))
            B,A = eval_to_linear!(B,A,temp,matrices1,tweak(J - (i+1)*V,d*n-n) - y,y)
        end
        
        my_add!(temp,A,B)
        gMat = temp*gMat

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

# what in here could possibly be causing a lock conflict?
"""
Iteratively execute reduction chunks in the navie strategy until
the vector g is reduced to pole order n
"""
function reducechain_naive(u,g,m,S,f,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    (9 < params.verbose) && println("u = $u") 
    

    #if cache.vars_reversed
    #    reverse!(u)
    #end
    
    J = rev_tweak(u,n*d-n)

    #if cache.vars_reversed
    #    reverse!(u)
    #end

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
        V = chooseV(J, d, S)
        (4 < params.verbose) && print("Chose V = $V; ")
        (6 < params.verbose) && begin
            # the way that chooseV works right now,
            # the following if statement should never hit.
            for i in 1:length(V)
                if params.vars_reversed && V[i] == 0 && J[i] ≠ 0 && (n+1-i) ∈ S
                    print("Illegal choice of V!")
                    println("J = $J, S = $S")
                end
            end
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

        #(5 < params.verbose && firsttime) && begin 
        #    for i in 1:length(matrices)
        #        println(matrices[i][:,end])
        #    end
        #end
        (6 < params.verbose && V == [1,1,2] && firsttime) && begin println(matrices); firsttime=false end

        #eval_to_linear!(context.B,context.A,context.temp,matrices,reverse(mins),reverse(V))
        eval_to_linear!(context.B,context.A,context.temp,matrices,mins,V)
        #(5 < params.verbose && firsttime) && begin 
        #    println("---")
        #    println(context.A[:,end])
        #    println(context.B[:,end])
        #end

        i = 1
        if params.fastevaluation == false
            params.verbose == 11 && println("gMat before is $gMat")
            while i <= K
                gMat = (context.A+context.B*(K-i))*gMat
                i = i+1
                if params.verbose == 11
                    A = context.A
                    B = context.B
                    #println("context.A is $A")
                    #println("context.B is $B")
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

function reducechain_varbyvar(u,g,m,S,f,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))

    J = copy(u)

    gMat = g
    mins = copy(J)
    tempv = copy(J)

    (4 < params.verbose) && println("Starting: J = $J")
    (5 < params.verbose) && begin
        g_poly = vector_to_polynomial(g,n,d*n-n,PR,params.termorder)
        if params.always_use_bigints || params.use_gpu
            println("Starting: g = $((gMat)) = $g_poly")
        else    
            println("Starting: g = $(Int.(gMat)) = $g_poly")
        end
    end

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
            for i in 1:length(V)
                if params.vars_reversed && V[i] == 0 && J[i] ≠ 0 && (n+1-i) ∈ S
                    print("Illegal choice of V!")
                    println("J = $J, S = $S")
                end
            end
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

        matrices = context.Ruvs[V]
        eval_to_linear!(context.B,context.A,context.temp,matrices,mins,V)
        gMat = finitediff_prodeval_linear!(context.B,context.A,0,K-1,gMat,context.temp,context.g_temp)
        @. J = J - K*V
        m = m - K
        if m == n
            highpole = false
        end

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
    return ((J, gMat), m)
end

"""
    costadata_of_intial_term
Returns the data used by costa's code given a polynomial term.
More specifically, given a polynomial h, it writes h in the form 
h = x^{rev_tweak(U)}g

Note: this only works with a single term, so it should
only be used at the beginning of reduction
"""
function costadata_of_initial_term!(term,g,n,d,p,S,cache,params)
    termorder = params.termorder
    vars_reversed = params.vars_reversed

    R = base_ring(parent(term[1])) 
    i = term

    Stilda = zeros(Int, n+1)
    for j in S
        Stilda[j+1] = 1
    end

    #ss = zeros(Int,n+1)
    #ss[S .+ 1] .= 1
    #o = ones(Int,length(exponent_vector(i[1],1)))
    if vars_reversed
        U = reverse(exponent_vector(i[1],1) + Stilda)
    else
        U = exponent_vector(i[1],1) + Stilda
    end
    
    gCoeff = coeff(i[1],1)

    if vars_reversed 
        g_exps = U - rev_tweak(U, n*d-n)  
    else
        g_exps = U - tweak(U, n*d-n)  # TODO: make our code work with tweak 
    end 
    ev = cache[n*d-n]#gen_exp_vec(n+1,n*d-n,termorder)
    # this is UInt instead of R to get Oscar to use the fast FLINT method
    #g = zeros(R,length(ev)) 
    my_zero!(g)

    for j in axes(g,1)
        if vars_reversed && !cache.vars_reversed
            if g_exps == reverse(ev[j])
                CUDA.@allowscalar g[j] = lift(gCoeff)
                break
            end 
        elseif vars_reversed && cache.vars_reversed
            if g_exps == ev[j]
                CUDA.@allowscalar g[j] = lift(gCoeff)
                break
            end
        else
            if g_exps == ev[j]
                CUDA.@allowscalar g[j] = lift(gCoeff)
                break
            end
        end 
    end

    #println(U,Int.(g))

    return (U,g)
end

"""
Takes an array of Costa's data tuples, and a new term
to be added. If the new term already exists, it is 
added, if not then it concatenates.

Note: perhaps this would be better served by a custom struct?

Note: right now, we don't have a method that just
takes an array of Costa's data tuples and consolodates it. 
Perhaps this also would be well suited by a custom struct?

^^ such concerns have speed implications, so better to wait
until we're really trying to optimize this.

"""
function incorporate_initial_term!(costadata_arr,costadata)
    ind_already = false

    # if the u is already there, just add the vectors
    for i in eachindex(costadata_arr)
        (u,g) = costadata_arr[i]
        if all(u .== costadata[1])
            costadata_arr[i] = (u, g .+ costadata[2])
            ind_already = true
        end
    end

    # otherwise, push on a new term
    if !ind_already
        #(9 < verbose) && println("incorporation: u has type $(typeof(costadata[1]))")
        push!(costadata_arr,costadata)
    end
end

function remove_duplicates!(costadata_arr)
    i = 1
    while i <= (length(costadata_arr)-1)
        j = i+1
        while j <= length(costadata_arr)
            if all(costadata_arr[i][1][1] .== costadata_arr[j][1][1])
                costadata_arr[i] = ((costadata_arr[i][1][1], costadata_arr[i][1][2] .+ costadata_arr[j][1][2]),costadata_arr[i][2][1])
                deleteat!(costadata_arr,j)
            else
                j = j + 1
            end
        end
        i = i + 1
    end
end


"""
Converts a Costa's data tuple to a polynomial with pole
after reduction.

Thus, the pole order is n.
"""
function poly_of_end_costadata(costadata,PR,p,d,n,params)
    (u,g_vec) = costadata
    vars = gens(PR)

    cpu_g = Array(g_vec)
    g = vector_to_polynomial(cpu_g,n,d*n-n,PR,params.termorder)

    (5 < params.verbose) && begin
        if params.fastevaluation && !params.use_gpu
            println("$(Int.(g_vec)) --> $g")
        else
            println("$(g_vec) --> $g")
        end
    end

    # no need to do rev_tweak since reducechain_costachunks returns the "true" u
    # on the last run
    if params.vars_reversed
        return [prod(vars .^ reverse(u)) * g, n]
    else
        return [prod(vars .^ u) * g, n]
    end
end


"""
Converts an array of Costa's data tuples to an array of polynomials with pole
after reduction.

Thus, the pole order is always n.

"""
function poly_of_end_costadatas(costadatas,PR,p,d,n,S,params)
    res = PR(0)
    vars = gens(PR)
    for costadata in costadatas
        res += poly_of_end_costadata(costadata,PR,p,d,n,params)[1]
    end
    XS =  prod(PR(vars[i+1]) for i in S; init = PR(1))
    [[div(res,XS), n]]
end

"""
    reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruvs,termorder)

Implements Costa's algorithm for controlled reduction,
sweeping down the terms of the series expansion by the pole order.
"""
function reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,params)
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
            ω[i] = reducechain_costachunks(ω[i]...,poleorder,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,g_temp,params)
        end

        poleorder = poleorder - p
    end

    #println("ω: $ω")
    #(9 < verbose) && println(poly_of_end_costadatas(ω,PR,p,d,n,S,termorder))

    #println(gen_exp_vec(n,n*d-n-1,termorder))
           

    return poly_of_end_costadatas(ω,PR,p,d,n,S,params)
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

    allcostadata = []
    for term in terms
        g = copy(context.g)
        term_costadata = costadata_of_initial_term!(term,g,n,d,p,S,cache,params)
        append!(allcostadata,[((rev_tweak(term_costadata[1],n*d-n),term_costadata[2]),term[2])])
    end

    notallred = true
    while notallred
        for i in eachindex(allcostadata)
            allcostadata[i] = reducechain_varbyvar(allcostadata[i][1]...,allcostadata[i][2],S,f,p,context,cache,params)
        end
        remove_duplicates!(allcostadata)
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

function reducepoly_naive(pol,S,f,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    result = PR()
    for term in pol
        #println("Reducing terms of order $(term[2])")
        terms = termsoforder(pol,term[2])
        #println(terms)
        for t in terms
            (u,_) = costadata_of_initial_term!(t,context.g,n,d,p,S,cache,params)
            reduced = reducechain_naive(u,context.g,t[2],S,f,p,context,cache,params)
            (reduced_poly,m) = poly_of_end_costadata(reduced,PR,p,d,n,params)
            @assert m == n "Controlled reduction outputted a bad pole order"
            result += reduced_poly
        end
    end

    vars = gens(PR)
    XS =  prod(PR(vars[i+1]) for i in S; init = PR(1))
    [[div(result,XS), n]]
    #return poly_of_end_costadatas(result,PR,p,d,n,S,params)
end

"""
    reducetransform_costachunks(FT,N_m,S,f,pseudoInverseMat,p,cache,params)

trying to emulate Costa's controlled reduction, changes the order that polynomials are reduced, starts from highest pole order and accumulates the lower order poles as reduction proceeds

N_m - the precision
"""
function reducetransform_costachunks(FT,N_m,S,f,pseudoInverseMat,p,cache,params)

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
            @time reduction = reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,params)
        else
            reduction = reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruv,cache,A,B,temp,params)
        end

        push!(result, reduction)
    end

    return result
end

function reducetransform_varbyvar(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g_length = binomial(d*n,d*n-n)

    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)

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


    #TODO: right now, it usually seems better to do lazy computations,
    #  since not all of the Ruv are used. However, I know that for some
    #  classes of examples, they are all pretty much always used. For 
    #  such examples, it's better to use an EagerPEP and do threads.
    lazy_Ruv = length(S) < d || d < n

    if (0 < params.verbose)
        println("Creating the Ruv PEP object...")
        #CUDA.@time Ruv = select_Ruv_PEP(params,computeRuv,computeRuv_gpu,lazy_Ruv,MS1,cache,d)
        CUDA.@time Ruv = select_Ruv_PEP(n,d,S,params,computeRuv,lazy_Ruv,MS1,cache)
    else
        Ruv = select_Ruv_PEP(n,d,S,params,computeRuv,lazy_Ruv,MS1,cache)
    end

    result = similar(FT)

    #context_tlv = OhMyThreads.TaskLocalValue{default_context_type(MS1,params)}(
    #    () -> default_context(MS1,Ruv,params)
    #)
    #Threads.@threads for i in 1:length(FT) 
    #    context = context_tlv[]

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

    #(0 < params.verbose && Ruv isa CachePEP) && begin
    #    println("Ruv cache info: $(cache_info(Ruv.Ucomponent))")
    #end
    (0 < params.verbose) && begin
        println("Created $(length(allpoints(Ruv))) of $(length(cache[d])) possible V")
    end
    (1 < params.verbose) && begin
        println("V that were created: \n$(allpoints(Ruv))")
    end

    return result
end


function cuMod(A::zzModMatrix)
    m = modulus(base_ring(parent(A)))
    # if this conversion is illegal, 
    #then we're not allowed to make a CuModMatrix
    M = Int(m) 

    float_A = float_entries(A)
    CuModMatrix(float_A,M,elem_type=Float64)
end

function KaratsubaMat(A::zzModMatrix)
    m = modulus(base_ring(parent(A)))
    temp = collect(factor(m))[1]
    d = temp[2]
    p = temp[1]

    N1 = Int(p)^Int(round(d/2))
    N2 = Int(p)^Int(d - round(d/2))

    GPUFiniteFieldMatrices.KaratsubaMatrix(Float64,cuMod(A),N1,N2,N1^2)
end

#TODO: incorporate this idiomatically using adapt.jl
function adapt_to_gpu(Ruvs)
    gpu_Ruvs = Dict{Vector{Int64},Vector{CuModMatrix{Float64}}}()

    for (v,linpoly) in Ruvs

        gpu_Ruvs[v] = Vector(undef,length(linpoly)) 
        for i in 1:length(linpoly)
            gpu_Ruvs[v][i] = cuMod(linpoly[i])
        end
    end

    gpu_Ruvs
end

function default_context(matspace,Ruv,params)
    g_length = number_of_rows(matspace)
    B = base_ring(matspace)
    m = modulus(B)
    #if params.use_gpu == true && (ZZ(2)^22 < m < ZZ(2)^106)
    if params.use_gpu == true && (m < ZZ(2)^106)
        temp = collect(factor(m))[1]
        d = temp[2]
        p = temp[1]

        N1 = Int(p)^Int(round(d/2))
        N2 = Int(p)^Int((d - round(d/2)))

        A = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1^2,true)
        B = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1^2,true)
        temp = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1^2,true)
        GPUFiniteFieldMatrices.initialize_plan!(A)
        GPUFiniteFieldMatrices.initialize_plan!(B)
        GPUFiniteFieldMatrices.initialize_plan!(temp)

        g = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,N1,N2,N1^2,true)
        g_temp = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,N1,N2,N1^2,true)
        GPUFiniteFieldMatrices.initialize_plan!(g)
        GPUFiniteFieldMatrices.initialize_plan!(g_temp)

        return ControlledReductionContext{KaratsubaMatrix{Float64},KaratsubaVector{Float64}}(Ruv,A,B,temp,g,g_temp)
    elseif params.use_gpu == true

        M = Int(m)
        A = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)
        B = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)
        temp = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)

        g = GPUFiniteFieldMatrices.zeros(Float64,g_length,M)
        g_temp = GPUFiniteFieldMatrices.zeros(Float64,g_length,M)

        return ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}(Ruv,A,B,temp,g,g_temp)
    else # use the CPU
        A = matspace()
        B = matspace()
        temp = matspace()

        if params.always_use_bigints || ZZ(2)^64 < ZZ(m)
            # Big modulus
            # NOTE: can't use `similar` for a pointer type
            g = [ZZ(0) for i in 1:g_length]
            g_temp = [ZZ(0) for i in 1:g_length]
        else
            # Small modulus
            g = zeros(UInt,g_length) 
            g_temp = similar(g)
        end

        return ControlledReductionContext(Ruv,A,B,temp,g,g_temp)
    end
end

function default_context_type(matspace,params)
    B = base_ring(matspace)
    m = modulus(B)

    if params.use_gpu == true
        ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}
    elseif params.always_use_bigints || ZZ(2)^64 < ZZ(m)
        ControlledReductionContext{ZZModMatrix,Vector{ZZRingElem}}
    else
        ControlledReductionContext{zzModMatrix,Vector{UInt}}
    end
end

"""
Gives a special PEP for when S = [0] and d=3

This should be all the Ruv ever needed in this situation.
"""
function cubic_S_zero_Vs(n)
    Vs = Vector{Vector{Int}}()

    for i in 1:n+1
        V = zeros(Int,n+1)
        V[i] = 3
        push!(Vs,V)
    end

    for i in 1:n
        V1 = zeros(Int,n+1)
        V1[i] = 1
        V1[i+1] = 2
        push!(Vs,V1)

        V2 = zeros(Int,n+1)
        V2[i] = 2
        V2[i+1] = 1
        push!(Vs,V2)
    end

    Vs
end

"""
This function exists so that we can have the option to time it
in verbose mode.
"""
function select_Ruv_PEP(n,d,S,params,compute,lazy,oscar_matspace,cache)

    m = Integer(modulus(base_ring(oscar_matspace)))

    if (3 < params.verbose)
        compute_gpu = V -> begin println("Moving Ruv for $V to the gpu"); cuMod.(compute(V)) end
    else
        compute_gpu = V -> cuMod.(compute(V))
    end 

    compute_float = V -> float_entries.(compute(V))

    if 3 < n && d == 3 && S == [0] && params.use_gpu #4 < n

        eager_Vs = cubic_S_zero_Vs(n)

        cpu_Ruv = LazyPEP{Matrix{Float64}}(compute_float,eagerVs=eager_Vs,usethreads=false)
        m = M = Int(modulus(base_ring(oscar_matspace)))

        if (3 < params.verbose)
            create_gpu = matrices -> begin
                println("Moving an Ruv to the GPU (first time creation!)")
                CUDA.@time CuModMatrix.(matrices,M,elem_type=Float64)
            end

            convert_gpu = (dest,matrices) -> begin
                println("Moving an Ruv to the GPU (already allocated)")

                CUDA.@time begin
                    for i in 1:length(dest)
                        copyto!(dest[i],matrices[i])
                    end
                end
            end
        else
            create_gpu = matrices -> CuModMatrix.(matrices,M,elem_type=Float64)
            convert_gpu = (dest,matrices) -> begin
                for i in 1:length(dest)
                    copyto!(dest[i],matrices[i])
                end
            end
        end
        s = size(oscar_matspace(),1)

        memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
        # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
        memory = s^2 * 8 * (n+2)

        maxsize = div(memory_cap,memory)

        #println(maxsize)
        #testing
        #maxsize = 13 

        Ruv = CachePEP{Matrix{Float64},CuModMatrix{Float64}}(cpu_Ruv,create_gpu,convert_gpu,maxsize)

        #(1 < params.verbose) && println("Initial cache info: $(cache_info(Ruv.Ucomponent))")
    #elseif params.use_gpu && lazy && (ZZ(2)^22 < m < ZZ(2)^106)
    elseif params.use_gpu && lazy && (m < ZZ(2)^106)
        compute_gpu_karatsuba = V -> KaratsubaMat.(compute(V))
        Ruv = LazyPEP{KaratsubaMatrix{Float64}}(compute_gpu_karatsuba)

    elseif params.use_gpu && lazy
        Ruv = LazyPEP{CuModMatrix{Float64}}(compute_gpu)

    #elseif params.use_gpu && (ZZ(2)^22 < m < ZZ(2)^106)
    elseif params.use_gpu && (m < ZZ(2)^106)
        compute_gpu_karatsuba = V -> KaratsubaMat.(compute(V))
        Ruv = EagerPEP{KaratsubaMatrix{Float64}}(cache[d],compute_gpu_karatsuba,usethreads=false)

    elseif params.use_gpu
        Ruv = EagerPEP{CuModMatrix{Float64}}(cache[d],compute_gpu,usethreads=false)

    elseif lazy
        Ruv = LazyPEP{typeof(oscar_matspace())}(compute)

    else
        Ruv = EagerPEP{typeof(oscar_matspace())}(cache[d],compute,usethreads=false)
    end

    Ruv
end


function reducetransform_naive(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g_length = binomial(d*n,d*n-n)

    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)

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


    #TODO: right now, it usually seems better to do lazy computations,
    #  since not all of the Ruv are used. However, I know that for some
    #  classes of examples, they are all pretty much always used. For 
    #  such examples, it's better to use an EagerPEP and do threads.
    lazy_Ruv = length(S) < d || d < n

    if (0 < params.verbose)
        println("Creating the Ruv PEP object...")
        #CUDA.@time Ruv = select_Ruv_PEP(params,computeRuv,computeRuv_gpu,lazy_Ruv,MS1,cache,d)
        CUDA.@time Ruv = select_Ruv_PEP(n,d,S,params,computeRuv,lazy_Ruv,MS1,cache)
    else
        Ruv = select_Ruv_PEP(n,d,S,params,computeRuv,lazy_Ruv,MS1,cache)
    end

    result = similar(FT)

    #context_tlv = OhMyThreads.TaskLocalValue{default_context_type(MS1,params)}(
    #    () -> default_context(MS1,Ruv,params)
    #)
    #Threads.@threads for i in 1:length(FT) 
    #    context = context_tlv[]

    context = default_context(MS1,Ruv,params)
    for i in 1:length(FT) #pol in FT
    

        pol = FT[i]
        if (0 < params.verbose)
            println("Reducing vector $i")
            @time reduction = reducepoly_naive(pol,S,f,p,context,cache,params)
        else
            reduction = reducepoly_naive(pol,S,f,p,context,cache,params)
        end
        result[i] = reduction

        #println("cache info: $(cache_info(Ruv.Ucomponent))")
        #i == 5 && error("stopping after vector $i for testing purposes")
        
        #push!(result, reduction)
    end

    #(0 < params.verbose && Ruv isa CachePEP) && begin
    #    println("Ruv cache info: $(cache_info(Ruv.Ucomponent))")
    #end
    (0 < params.verbose) && begin
        println("Created $(length(allpoints(Ruv))) of $(length(cache[d])) possible V")
    end
    (1 < params.verbose) && begin
        println("V that were created: \n$(allpoints(Ruv))")
    end

    return result
end

function reducetransform(FT,N_m,S,f,pseudoInverseMat,p,params,cache)
    if params.algorithm == :costachunks
        reducetransform_costachunks(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    elseif params.algorithm == :naive
        reducetransform_naive(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    elseif params.algorithm == :varbyvar
        reducetransform_varbyvar(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    else
        throw(ArgumentError("Unsupported Algorithm: $algorithm"))
    end
end

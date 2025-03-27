#module ControlledReduction

#using Oscar
#using BitIntegers
#using LinearAlgebra
#using Combinatorics

#include("PrecisionEstimate.jl")
#include("CopiedFindMonomialBasis.jl")
#include("Utils.jl")
#include("PolynomialWithPole.jl")
#include("SmallestSubsetSmooth.jl")


#verbose = false

"""

KeyType - the type of u and v, for example Vector{Int64}
MatrixType - the type of A and B, for example zzModMatrix
VectorType - the type of g, for exampel Vector{UInt64}


#const OscarSmallReductionContext = ControlledReductionContext{Vector{Int64},
#                                                         zzModMatrix,
#                                                         Vector{UInt64}}
#const OscarBigReductionContext = ControlledReductionContext{Vector{Int64},
#                                                         ZZModMatrix,
#                                                         Vector{ZZRingElem}}

In the future: make a GPUReductionContext

"""
struct ControlledReductionContext{KeyType,MatrixType,VectorType}
    Ruvs::Dict{KeyType, Vector{MatrixType}} #TODO: this type might be too rigid
    A::MatrixType
    B::MatrixType
    temp::MatrixType
    g::VectorType
    g_temp::VectorType
end

#
#function oscar_default_context(R,g_length,n)
#    MS1 = matrix_space(R, g_length, g_length)
#    Ruvs = Dict{Vector{Int64}, Vector{typeof(MS1())}}()
#    A = MS1()
#    B = MS1()
#    temp = MS1()
#    g = zeros(UInt,g_length) 
#    g_temp = similar(g)
#
#    OscarReductionContext(Ruvs,A,B,temp,g,g_temp)
#end

### The following is derived from Nemo.jl, and needs to be licensed with the GPL 

function my_addmul!(A::zzModMatrix, B::zzModMatrix, C::UInt, D::zzModMatrix)
  @ccall Oscar.Nemo.libflint.nmod_mat_scalar_addmul_ui(A::Ref{zzModMatrix}, 
                                                       B::Ref{zzModMatrix}, 
                                                       D::Ref{zzModMatrix}, 
                                                       C::UInt)::Nothing
  return A
end

function my_mul!(A::zzModMatrix, B::zzModMatrix, c)
    p = characteristic(base_ring(parent(A)))
    ui(i) = 0 ≤ i ? UInt(i) : UInt(i + p)
    c = ui(c)
    @ccall Oscar.Nemo.libflint.nmod_mat_scalar_mul(A::Ref{zzModMatrix},
                                                   B::Ref{zzModMatrix},
                                                   c::UInt)::Nothing
    return A
end

function my_mul!(A::ZZModMatrix, B::ZZModMatrix, c::Int)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_scalar_mul_si(A::Ref{ZZModMatrix}, 
                                               B::Ref{ZZModMatrix}, 
                                               c::Int, 
                                               base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return A 
end

function my_mul!(A::ZZModMatrix, B::ZZModMatrix, c::ZZRingElem)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_scalar_mul_fmpz(A::Ref{ZZModMatrix}, 
                                                 B::Ref{ZZModMatrix},
                                                 C::Ref{ZZRingElem},
                                                 base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return A
end

function my_matvecmul!(z::Vector{UInt},A::zzModMatrix,b::Vector{UInt})
    @ccall Oscar.Nemo.libflint.nmod_mat_mul_nmod_vec(z::Ptr{UInt},
                                                     A::Ref{zzModMatrix},
                                                     b::Ptr{UInt},
                                                     length(b)::Int)::Nothing
    return z
end

function my_matvecmul!(z::Vector{ZZRingElem},A::ZZModMatrix,b::Vector{ZZRingElem})
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_mul_fmpz_vec_ptr(z::Ptr{Ref{ZZRingElem}}, 
                                                  A::Ref{ZZModMatrix}, 
                                                  b::Ptr{Ref{ZZRingElem}}, 
                                                  length(b)::Int, 
                                                  base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return z
end

#function alt_matvecmul!(z::Vector{ZZRingElem},A::ZZModMatrix,b::Vector{ZZRingElem})
#    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_mul_fmpz_vec(z::Ref{ZZRingElem},
#                                                         A::Ref{ZZModMatrix},
#                                                         b::Ref{ZZRingElem},
#                                                         length(b)::Int,
#                                                         base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
#end

### END stuff derived from Nemo.jl

function my_copy!(a,b)
    copy!(a,b)
end

function my_copy!(a::Vector{ZZRingElem},b::Vector{ZZRingElem})
    for i in 1:length(a)
        Oscar.Nemo.set!(a[i],b[i])
    end
end

function my_zero!(a)
    zero!(a)
end


function my_zero!(a::Vector{ZZRingElem})
    for i in 1:length(a)
        Oscar.Nemo.set!(a[i],0)
    end
end

function my_sub!(A,B,C)
    Oscar.Nemo.sub!(A,B,C)
end

function my_sub!(A::ZZModMatrix,B::ZZModMatrix,C::ZZModMatrix)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_sub(A::Ref{ZZModMatrix},
                                                B::Ref{ZZModMatrix},
                                                C::Ref{ZZModMatrix},
                                                base_ring(B).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return A
end


"""
    reduce_LA(U,V,S,f,pseudoInverseMat,g,PR,termorder)

applies reduction formula from Prop 1.15 in Costa's thesis to 
basis elements of Homog(dn-d), returns them as polynomials
will only work with vars_reversed=true
"""
function reduce_LA(U,V,S,f,pseudoInverseMat,g,PR,termorder)
    R = coefficient_ring(PR)
    Vars = gens(PR)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    SC = []
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), V)
    XV = finish(B)
    for i in 0:n
        if i in S
        else
            push!(SC,i)
        end
    end
    # get gi's using pseudoinverse
    XS =  prod(PR(Vars[i+1]) for i in S; init = PR(1))
    gVec = convert_p_to_m([div(XV*(g[1]),XS)],gen_exp_vec(n+1,n*d-n+d-length(S),termorder))
    MS = matrix_space(parent(gVec[1]), nrows(pseudoInverseMat),1)
    gJS = MS()
    gJS = pseudoInverseMat*transpose(gVec)
    gc = []
    for i in 1:(n+1)
        push!(gc, convert_m_to_p(transpose(gJS[Int((i-1)*(length(gJS)/(n+1))+1):Int(i*(length(gJS)/(n+1))),:]),gen_exp_vec(n+1,n*d-n-d+1,termorder),R,PR)[1])
    end
    gc = reverse(gc)
    gcpartials = [ derivative(gc[i], i) for i in 1:(n+1) ]
    
    reverse!(gcpartials) # TODO: make this an option, this is the way it is in Costa's code, 

    #return [sum(PR(U[i+1])*XS*gc[i+1] + div(XS,Vars[i+1])*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]
    return [sum(PR(U[i+1])*div(XS,Vars[i+1])*gc[i+1] + XS*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0))]
c
end

"""
    chooseV(I, d)
Choose the direction of reduction V following Edgar's method

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f 
"""
function chooseV(I, d)
    v = zeros(Int,length(I))
    n = length(I)

    sum = 0
    for i in 1:n 
        if I[i] > 0
            v[i] = v[i] + 1
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
    rev_tweak(I,m)

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
function reducechain_costachunks(u,g,m,S,f,pseudoInverseMat,p,Ruvs,cache,A,B,temp,g_temp,params)
    #p = Int64(characteristic(parent(f)))
    verbose = params.verbose
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    
#    ui(i) = 0 ≤ i ? UInt(i) : UInt(i + characteristic(R))

    
    I = u
    #TODO?
    if params.vars_reversed == false
        reverse!(I) # parity issue due to Costa's code being reverse from ours
    end
    (4 < verbose) && println("Expanded I: $I")

    gMat = g
    #println(gMat)
    #chain = 0
    #I_edgar = [x//7 for x in I]
    #(4 < verbose) && println("This is I: $I_edgar")
    J = copy(I)

    #TODO?
    if params.vars_reversed == false
        V = rev_chooseV(Array{Int}(divexact.(I,p)),d)
    else
        V = chooseV(Array{Int}(divexact.(I,p)),d)
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

    matrices = computeRuvS(V,S,f,pseudoInverseMat,Ruvs,cache,params)

    #TODO: the following was changed to use reverse! so 
    #    it doesn't allocate as much, but I realized that
    #    this isn't our bottleneck. If it ever does
    #    become the bottleneck, fix the rest of this method
    #    so it doesn't allocate.
    U = I .- (nend-(d*n-n))*V
    reverse!(U)

    reverse!(V)
    B,A = computeRPoly_LAOneVar2!(B,A,matrices,U,V,R,temp)
    reverse!(V) # put V back to normal

    i = 1

    
    #(9 < verbose) && println("Before reduction chunk: $(convert.(Int,gMat))")
    (4 < verbose) && println("Before reduction chunk, I is $I")
    if params.fastevaluation && 1 ≤ nend-(d*n-n)
      gMat = finitediff_prodeval_linear!(B,A,0,nend-(d*n-n)-1,gMat,temp,ui)
      i = nend-(d*n-n) + 1
    else
      while i <= (nend-(d*n-n))
        my_mul!(temp,B,nend-(d*n-n)-i)
        add!(temp,temp,A)
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
        
        matrices1 = computeRuvS(y,S,f,pseudoInverseMat,Ruvs,cache,params)
        #println(matrices1)
        #error()

        if params.vars_reversed == true
            B,A = computeRPoly_LAOneVar2!(B,A,matrices1,reverse(rev_tweak(J - (i+1)*V,d*n-n) - y),reverse(y),R,temp)
        else
            B,A = computeRPoly_LAOneVar2!(B,A,matrices1,reverse(tweak(J - (i+1)*V,d*n-n) - y),reverse(y),R,temp)
        end
        
        add!(temp,A,B)
        gMat = temp*gMat

        #gMat = (A+B)*gMat

        (4 < verbose) && println("After step $(i+1): $(gMat))")
        

        i = i+1
        @. I = I - y
        (4 < verbose) && println("After step $(i+1), I is $I")
    end

      #U = I - V
    #=
    K = 0
    mins = I
    while true
        temp = mins - V
        isLessThanZero = false
        for j in temp
            if j < 0
                isLessThanZero = true
                break
            end
        end
        if isLessThanZero == true
            break
        end
        mins = temp
        K = K+1
    end
    =#


    #A,B = computeRPoly_LAOneVar(V,I - Int64((nend-(d*n-n)))*V,S,n,d,f,pseudoInverseMat,R,PR,termorder)
    #=
    matrices = computeRPoly_LAOneVar1(V,S,f,pseudoInverseMat,Ruvs,termorder)
    for i in axes(matrices,1)
        printMat(matrices[i])
    end
    =#
    
    #=
    MK = A + B*K
    MK1 = A + B*(K-1)
    h = MK1*MK*gMat
    A1 = MK - MK1
    j = 2
    while K-j >= 0
        MKj = MK - A1*j
        h = MKj*h
        j = j + 1
    end
            m = m - K
    I = mins
    =# 
    #throw(error)
    
    #error()
    
    if nend == p
        newI = J .- p*V
        #@assert undo_rev_tweak(I,p) == newI

        return (newI, gMat)
    else
        return (I,gMat) # gives the "true" u
    end
    
end

function reducechain_naive(u,g,m,S,f,pseudoInverseMat,p,context,cache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    #ui(i) = 0 ≤ i ? UInt(i) : UInt(i + characteristic(R))
    J = rev_tweak(u,n*d-n)
    gMat = context.g
    mins = similar(J)
    tempv = similar(J)
    (4 < params.verbose) && println("Starting: J = $J")
    (5 < params.verbose) && begin
        if params.always_use_bigints
            println("Starting: g = $((gMat))")
        else
            println("Strting: g = $(Int.(gMat))")
        end
    end

    firsttime = true


    while m > n
        V = chooseV(J,d)
        (4 < params.verbose) && print("Chose V = $V; ")
        (6 < params.verbose) && begin
            # the way that chooseV works right now,
            # the following if statement will never hit.
            for i in 1:length(V)
                if V[i] == 0 && J[i] ≠ 0 && (n+1-i) ∈ S
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
        matrices = computeRuvS(V,S,f,pseudoInverseMat,context.Ruvs,cache,params)

        #(5 < params.verbose && firsttime) && begin 
        #    for i in 1:length(matrices)
        #        println(matrices[i][:,end])
        #    end
        #end
        (6 < params.verbose && V == [1,1,2] && firsttime) && begin println(matrices); firsttime=false end

        computeRPoly_LAOneVar2!(context.B,context.A,matrices,reverse(mins),reverse(V),R,context.temp)
        #(5 < params.verbose && firsttime) && begin 
        #    println("---")
        #    println(context.A[:,end])
        #    println(context.B[:,end])
        #end

        i = 1
        if params.fastevaluation == false
            while i <= K
                gMat = (context.A+context.B*(K-i))*gMat
                i = i+1
            end
        else
            gMat = finitediff_prodeval_linear!(context.B,context.A,0,K-1,gMat,context.temp,context.g_temp)
        end
        @. J = J - K*V
        m = m - K
        (4 < params.verbose) && print("After $(lpad(K,4,' ')) steps,")
        (4 < params.verbose) && println("J = $J")
        if (5 < params.verbose) 
            g = vector_to_polynomial(gMat,n,d*n-n,PR,params.termorder)
            if params.always_use_bigints
                println("g = $((gMat)) = $g")
            else
                println("g = $(Int.(gMat)) = $g")
            end
        end
        
    end
    return (J, gMat)
end

"""
finitediff_prodeval_linear!(a,b,start,stop,g,temp)

Generic function that evaluates 
the function F(x) = a*x + b at 
start, ..., stop and 
computes
F(start)*...*F(stop)*g

Since it's a generic function, it'll work if
a and b and g are numbers, but the intended 
use is for a and b to be matrices and g
to be a vector.

My hope is that since this is generic, we'll
be able to replace arrays with CuArrays for large examples
and it'll still work.

a - linear term coefficient
b - constant term coefficient
start - lowest value to evaluate at (should be an integer)
stop - highest value to evaluate at (should be an integer)
g - the value into which the answer is accumulated.
temp - a temporary pre-allocated matrix that can be used in intermediate computations
ui - a function that converts the output of Julia's mod into an unsigned integer 
    (i.e. if it's negative, add the characteristic)

"""
function finitediff_prodeval_linear!(a,b,start,stop,g,temp,g_temp,ui=nothing)

  #zero!(temp) # temp = F_k

  #Fk = stop*a + b # Fk = F(k), here k = stop
  my_mul!(temp,a,stop)
  add!(temp,temp,b)
  # @. Fk = stop*a + b
  #
  if start == stop
    #Oscar.mul!(g,temp,g)
    #g = temp * g
    #my_matvecmul!(g,temp,g)
    my_matvecmul!(g_temp,temp,g)
    #alt_matvecmul!(g_temp,temp,g)
    #g,g_temp = g_temp,g
    my_copy!(g,g_temp)
    return g
  end


  #g = temp*g
  #gg = copy(g)
  #h = temp * gg
  my_matvecmul!(g_temp,temp,g)
  #alt_matvecmul!(g_temp,temp,g)
  #g,g_temp = g_temp,g
  #copy!(g,g_temp)
  my_copy!(g,g_temp)

  #if h != g
  #    
  #    println(gg)
  #    println(temp)
  #    println("mul! - $g")
  #    println("naive * - $h")
  #    error()
  #end

  for k = stop-1:-1:start
    # right now, Fk = F(k+1)
    
    # Fk = Fk - a
    my_sub!(temp,temp,a)

    # now, Fk = F(k)
    #g = temp * g
    #gg = copy(g)
    #h = temp * gg
    #Oscar.Nemo.mul!(g,temp,g)
    my_matvecmul!(g_temp,temp,g)
    #alt_matvecmul!(g_temp,temp,g)
    #g,g_temp = g_temp,g
    #copy!(g,g_temp)
    my_copy!(g,g_temp)
    #println(g)
    #println(h == g)
    #if h != g
    #    
    #    println(gg)
    #    println("mul! - $g")
    #    println("naive * - $h")
    #    error()
    #end
  end

  # we don't know whether g_temp or g is the "real" g anymore
  #my_copy!(g_temp,g)

  return g
end

"""
Returns the data used by costa's code given a polynomial term.

Note: this only works with a single term, so it should
only be used at the beginning of reduction
"""
function costadata_of_initial_term!(term,g,n,d,p,S,termorder)

    #(9 < verbose) && println("p: $p")

    R = base_ring(parent(term[1])) 
    #mod = modulus(R)
    i = term

    ss = zeros(Int,n+1)
    ss[S .+ 1] .= 1
    #o = ones(Int,length(exponent_vector(i[1],1)))
    II = exponent_vector(i[1],1) + ss
    gCoeff = coeff(i[1],1)

    #V = chooseV(Array{Int}(divexact.(II,p)),d)
    u = II - rev_tweak(II,n*d-n)
    ev = gen_exp_vec(n+1,n*d-n,termorder)
    # this is UInt instead of R to get Oscar to use the fast FLINT method
    #g = zeros(R,length(ev)) 
    #TODO: put a check to make see if we are using BigInts here
    my_zero!(g)

    for j in axes(g,1)
        if u == ev[j]
            g[j] = lift(gCoeff)
            break
        end
    end


    #(9 < verbose) && println("creation: u is type $(typeof(II))")
    return (II,g)
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

"""
Converts a Costa's data tuple to a polynomial with pole
after reduction.

Thus, the pole order is n.
"""
function poly_of_end_costadata(costadata,PR,p,d,n,params)
    (u,g_vec) = costadata
    vars = gens(PR)

    g = vector_to_polynomial(g_vec,n,d*n-n,PR,params.termorder)

    (5 < params.verbose) && println("$(Int.(g_vec)) --> $g")
    # no need to do rev_tweak since reducechain_costachunks returns the "true" u
    # on the last run
    [prod(vars .^ u) * g, n]
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

#"""
#given polynomial, splits into terms and applies reduction to each term
#"""
#function reducepoly_LA(poly,n,d,p,S,f,pseudoInverseMat,R,PR)
#    t = PolynomialWithPole.terms(poly)
#    result = 0
#    for i in t
#        println("reducing term $i")
#        o = ones(Int,length(exponent_vector(i[1],1)))
##        I = exponent_vector(i[1],1) + o
#        I = undo_rev_tweak(exponent_vector(i[1],1),p)
#        gCoeff = coeff(i[1],1)
#        RChain = reducechain_LA(I,gCoeff,n,d,p,i[2],S,f,pseudoInverseMat,R,PR)
#        B = MPolyBuildCtx(PR)
#        println("Result of reduction chunk: $RChain")
#        push_term!(B, R(1), Int.(RChain[2]))
#        XU = finish(B)
#        B = MPolyBuildCtx(PR)
#        push_term!(B, R(1), o)
#        XS = finish(B)
#        ev = gen_exp_vec(n+1,d*n - n,termorder)
#        gReduction = div(XU*convert_m_to_p(transpose(RChain[1]),ev,R,PR)[1],XS)
#        println("Polynomial obtained from reduction: $gReduction")
#        result = result + gReduction
#    end
#    return [result,poly[2] - min(p,poly[2]-n)]
#end
#
#"""
#naive controlled reduction, takes as input the array of frobenius transforms of basis elements and reduces each polynomial
#
#NOTE: what is big N and why isn't it used here?
#"""
#function reducetransform_LA(FT,n,d,p,N,S,f,pseudoInverseMat,R,PR)
#    result = []
#    for i in FT
#        temp = 0
#        for j in axes(i,1)
#            #t = reducepoly_LA([Factorial(R(i[length(i)][2]),R(j[2]))p^(j[2]-n-1)*(j[1]),j[2]],n,d,p,S,f,pseudoInverseMat,R,PR)[1]
#            t = i[j]
#            while t[2] > n
#                t = reducepoly_LA(t,n,d,p,S,f,pseudoInverseMat,R,PR)
#            end
#            temp = temp + t[1]
#        end
#        #push!(result,[temp,n,Factorial(Int1024(i[length(i)][2]),n)])
#        push!(result,[temp,n,Factorial(R(i[length(i)][2]),R(n))])
#    end
#    return result
#end

"""
    reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruvs,termorder)

Implements Costa's algorithm for controlled reduction,
sweeping down the terms of the series expansion by the pole order.
"""
function reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruvs,cache,A,B,temp,params)
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
            term_costadata = costadata_of_initial_term!(term,g,n,d,p,S,params.termorder)
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
            ω[i] = reducechain_costachunks(ω[i]...,poleorder,S,f,pseudoInverseMat,p,Ruvs,cache,A,B,temp,g_temp,params)
        end

        poleorder = poleorder - p
    end

    #println("ω: $ω")
    #(9 < verbose) && println(poly_of_end_costadatas(ω,PR,p,d,n,S,termorder))

    #println(gen_exp_vec(n,n*d-n-1,termorder))
           
    #error()

    return poly_of_end_costadatas(ω,PR,p,d,n,S,params)
end

function reducepoly_naive(pol,S,f,pseudoInverseMat,p,context,cache,params)
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
            (u,_) = costadata_of_initial_term!(t,context.g,n,d,p,S,params.termorder)
            reduced = reducechain_naive(u,context.g,t[2],S,f,pseudoInverseMat,p,context,cache,params)
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
    Ruvs = Dict{Vector{Int64}, Vector{typeof(MS1())}}()
    result = []
    i = 1
    for pol in FT
        (0 < params.verbose) && println("Reducing vector $i")
        i += 1
        if (0 < params.verbose)
            @time reduction = reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruvs,cache,A,B,temp,params)
        else
            reduction = reducepoly_costachunks(pol,S,f,pseudoInverseMat,p,Ruvs,cache,A,B,temp,params)
        end

        push!(result, reduction)
    end

    return result
end

function oscar_default_context(matspace,Ruvs,params)
    g_length = number_of_rows(matspace)
    B = base_ring(matspace)
    m = modulus(B)

    A = matspace()
    B = matspace()
    temp = matspace()

    if params.always_use_bigints || ZZ(2)^64 < ZZ(m)
       # Big modulus
       g = [ZZ(0) for i in 1:g_length]
       g_temp = [ZZ(0) for i in 1:g_length]
    else
        # Small modulus
        g = zeros(UInt,g_length) 
        g_temp = similar(g)
    end

    # NOTE: can't use `similar` for a pointer type
#    g_temp = zeros(eltype(g),length(g))
    ControlledReductionContext(Ruvs,A,B,temp,g,g_temp)
end

function reducetransform_naive(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g_length = binomial(d*n,d*n-n)

    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)

    Ruvs = Dict{Vector{Int64}, Vector{typeof(MS1())}}()

    #explookup = Dict{Vector{Int64}, Int64}()
    #ev1 = gen_exp_vec(n+1,n*d-n,params.termorder)
    #for i in 1:length(ev1)
    #    get!(explookup,ev1[i],i)
    #end

    if (0 < params.verbose)
        println("Calculating the R_uv...")
        @time precomputeRuvs(S,f,pseudoInverseMat,Ruvs,cache,params)
    else
        precomputeRuvs(S,f,pseudoInverseMat,Ruvs,cache,params)
    end
    
    contexts = ControlledReductionContext[]

    for i in 1:length(FT)#Threads.nthreads()
        push!(contexts, oscar_default_context(MS1,Ruvs,params))
    end

    result = similar(FT)

    #TODO: can reduce allocations by changing this for loop
    #  to a nested while inside for. Then only allocate one context
    #  thread, instead of one per reduction vector.
    #Threads.@threads for i in 1:length(FT) #pol in FT
    for i in 1:length(FT) #pol in FT
        context = contexts[i]
        pol = FT[i]
        (0 < params.verbose) && println("Reducing vector $i")
        if (0 < params.verbose)
            @time reduction = reducepoly_naive(pol,S,f,pseudoInverseMat,p,context,cache,params)
        else
            reduction = reducepoly_naive(pol,S,f,pseudoInverseMat,p,context,cache,params)
        end
        result[i] = reduction
        #push!(result, reduction)
    end

    return result
end

function reducetransform(FT,N_m,S,f,pseudoInverseMat,p,params,cache)
    if params.algorithm == :costachunks
        reducetransform_costachunks(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    elseif params.algorithm == :naive
        reducetransform_naive(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    else
        throw(ArgumentError("Unsupported Algorithm: $algorithm"))
    end
end

#"""
#computes reduction matrices
#"""
#function computeRPoly_LAOneVar(V,mins,S,n,d,f,pseudoInverseMat,R,PR,termorder)
#    YRing, y = polynomial_ring(R, "y")
#    PYRing, Vars = polynomial_ring(YRing, ["x$i" for i in 0:n])
#    yV = []
#    for i in axes(V,1)
#        push!(yV, y*V[i])
#    end
#    UVars = mins + yV
#    ev = gen_exp_vec(n+1,n*d-n,termorder)
#    monomials = gen_mon(ev,YRing,PYRing)
#    reductions = []
#    for m in monomials
#        push!(reductions, reduce_LA(UVars,V,S,f,pseudoInverseMat,[m,1],PYRing,termorder)[1])
#    end
#    polyMatrix = Matrix(transpose(convert_p_to_m(reductions,ev)))
#    matSpace = matrix_space(R,nrows(polyMatrix),ncols(polyMatrix))
#    A = matSpace()
#    B = matSpace()
#    for i in 1:nrows(polyMatrix)
#        for j in 1:ncols(polyMatrix)
#            A[i,j] = coeff(polyMatrix[i,j],0)
#            B[i,j] = coeff(polyMatrix[i,j],1)
#        end
#    end
#    return [A,B]
#end

function precomputeRuvs(S,f,pseudoInverseMat,Ruvs,cache,params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1

    evs = gen_exp_vec(n+1,d,params.termorder)
    #Threads.@threads for V in evs
    for V in evs
        computeRuvS(V,S,f,pseudoInverseMat,Ruvs,cache,params)
    end

    
    #(6 < params.verbose) && begin 
    #    println("V = [1,2,0] : $(Ruvs[[1,2,0]])")
    #    println("V = [0,2,1] : $(Ruvs[[0,2,1]])")
    #end
end

function computeRuvOld(V,S,f,pseudoInverseMat,Ruvs,cache,params)
    vars_reversed = params.vars_reversed
    termorder = params.termorder
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    R = coefficient_ring(parent(f))
    MS1 = matrix_space(R, binomial(n*d,n*d-n), binomial(n*d,n*d-n))
    if S != collect(0:n-1)
        throw(ArgumentError("computeRuvOld only works for S={0...n-1}"))
    end
    if haskey(Ruvs, V)
        return get(Ruvs, V, 0)
    else
        (4 < params.verbose) && println("New key: $V")
    end
    ev1 = cache[n*d-n]#gen_exp_vec(n+1,n*d-n,termorder)
    ev2 = cache[n*d-n+d-length(S)]#gen_exp_vec(n+1,n*d-n+d-length(S),termorder)
    ev3 = cache[n*d-n-length(S)+1]#gen_exp_vec(n+1,n*d-n-length(S)+1,termorder)
    explookup = cache[n*d - n,:reverse]
    temp = Vector{Int64}(undef, n+1)
    MS2 = matrix_space(R, length(ev2),1)
    result = Vector{typeof(MS1())}(undef, n+2)
    for i in 1:n+2
        result[i] = MS1(0)
    end
    Stilda = zeros(Int, length(S))
    for i in S
        Stilda[n+1-i] = 1
    end
    for i in 1:length(ev1)
        mon = Vector{Int64}(undef, n+1)
        for m in 1:(n+1)
            mon[m] = ev1[i][m] + V[m] - Stilda[m]
        end
        gVec = MS2()
        for j in 1:length(ev2)
            if ev2[j] == mon
                gVec[j] = R(1)
            else
                gVec[j] = R(0)
            end
        end
        gJS = pseudoInverseMat*gVec
        #println("After LingAlg problem: $gJS")
        for j in 1:(n+1)
            for k in 1:length(ev3)
                for m in 1:(n+1)
                    if m == n+1-j+1
                        temp[m] = ev3[k][m] + Stilda[m] - 1
                    else
                        temp[m] = ev3[k][m] + Stilda[m]
                    end
                end
                #print("ev1[l]: $((ev1[l],typeof(ev1[l])));")
                #print("ev3[k]: $((ev3[k],typeof(ev3[k])));") 
                #println(" $(ev1[l] == ev3[k])")
                l = get(explookup,temp,-1)
                #result[j+1][l,i] = gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
                #result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
                if params.vars_reversed == true
                    result[j+1][l,i] = gJS[Int((j-1)*(length(gJS)/(n+1))+1)+k-1,1]
                    result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((j-1)*(length(gJS)/(n+1))+1)+k-1,1]
                else
                    result[j+1][l,i] = gJS[Int((n+1-j)*(length(gJS)/(n+1))+1)+k-1,1]
                    result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((n+1-j)*(length(gJS)/(n+1))+1)+k-1,1]
                end
                #println("$(result[j+1][l,i]) in $(j+1) && $(result[1][l,i]) in $(1)")
            end
        end
    end
    get!(Ruvs, V, result)
    return result
end

function computeRuvS(V,S,f,pseudoInverseMat,Ruvs,cache,params)
    vars_reversed = params.vars_reversed
    termorder = params.termorder
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    R = coefficient_ring(parent(f))
    #MS1 = matrix_space(R, binomial(n*d,n*d-n), binomial(n*d,n*d-n))
    g_length = binomial(n*d,n*d-n)
    if haskey(Ruvs, V)
        return get(Ruvs, V, 0)
    else
        (4 < params.verbose) && println("New key: $V")
    end
    ev1 = cache[n*d - n]#gen_exp_vec(n+1,n*d-n,termorder)
    ev2 = cache[n*d-n+d-length(S)]#gen_exp_vec(n+1,n*d-n+d-length(S),termorder)
    ev3 = cache[n*d-n-length(S)+1]#gen_exp_vec(n+1,n*d-n-length(S)+1,termorder)
    ev4 = cache[n*d-n-length(S)]#gen_exp_vec(n+1,n*d-n-length(S),termorder)
    explookup = cache[n*d - n,:reverse]
    temp = Vector{Int64}(undef, n+1)
    MS2 = matrix_space(R, length(ev2),1)
    matrixtype = eltype(valtype(Ruvs))
    result = Vector{matrixtype}(undef, n+2)
    for i in 1:n+2
        result[i] = zero_matrix(R,g_length,g_length)
    end
    Stilda = zeros(Int, n+1)
    for i in S
        Stilda[n+1-i] = 1
    end
    distances = Vector{Int64}(undef, n+1)
    for i in 1:(n+1)
        if Stilda[i] == 1
            distances[i] = length(ev3)
        else
            distances[i] = length(ev4)
        end
    end
    distance = 0
    for i in 1:length(ev1)
        mon = Vector{Int64}(undef, n+1)
        for m in 1:(n+1)
            mon[m] = ev1[i][m] + V[m] - Stilda[n+1-m+1]
        end
        gVec = MS2()
        for j in 1:length(ev2)
            if ev2[j] == mon
                gVec[j] = R(1)
            else
                gVec[j] = R(0)
            end
        end
        gJS = pseudoInverseMat*gVec
        #println("After LingAlg problem: $gJS")
        distance = 0
        for j in 1:(n+1)
            if Stilda[j] == 1
                for k in 1:length(ev3)
                    for m in 1:(n+1)
                        if m == n+1-j+1
                            temp[m] = ev3[k][m] + Stilda[n+1-m+1] - 1
                        else
                            temp[m] = ev3[k][m] + Stilda[n+1-m+1]
                        end
                    end
                    #print("ev1[l]: $((ev1[l],typeof(ev1[l])));")
                    #print("ev3[k]: $((ev3[k],typeof(ev3[k])));") 
                    #println(" $(ev1[l] == ev3[k])")
                    l = get(explookup,temp,-1)
                    #result[j+1][l,i] = gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
                    #result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
                    result[j+1][l,i] = result[j+1][l,i] + gJS[distance+k,1]
                    result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[distance+k,1]
                    if V == [0,2,1] && ev1[i] == [0,0,4]
                        println("distance: $distance")
                        println("u_$(j) : $(gJS[distance+k,1]) for term $(ev3[k])")
                        println("constsant : $((ev3[k][n+1-j+1])*gJS[distance+k,1]) for term $(ev3[k])")
                    end
                    #println("$(result[j+1][l,i]) in $(j+1) && $(result[1][l,i]) in $(1)")
                end
            else
                for k in 1:length(ev4)
                    for m in 1:(n+1)
                        temp[m] = ev4[k][m] + Stilda[n+1-m+1]
                    end
                    #print("ev1[l]: $((ev1[l],typeof(ev1[l])));")
                    #print("ev3[k]: $((ev3[k],typeof(ev3[k])));") 
                    #println(" $(ev1[l] == ev3[k])")
                    l = get(explookup,temp,-1)
                    #result[j+1][l,i] = gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
                    #result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
                    result[j+1][l,i] = result[j+1][l,i] + gJS[distance+k,1]
                    result[1][l,i] = result[1][l,i] + (ev4[k][n+1-j+1])*gJS[distance+k,1] + gJS[distance+k,1]
                    #println("$(result[j+1][l,i]) in $(j+1) && $(result[1][l,i]) in $(1)")
                end
            end
            distance = distance + distances[j]
        end
    end
    #TODO: there is some sort of race condition on 
    # our dictionary, and putting this print statement here 
    # fixes it. Later, we'll need to fix this for reals
    (1 < Threads.nthreads()) && println("New key: $V")
    (2 < params.verbose) && begin
        nnzs = count.(!=(0),result)
        total_entries = g_length^2
        percentages = nnzs * 100 ./ total_entries 
        println("Ruv calculated for V = $V with densities")
        for i in 1:length(result)
            println("    Matrix $i: $(nnzs[i]) of $total_entries entreis, $(percentages[i]) %")
        end
    end
    get!(Ruvs, V, result)
    return result
end

#function computeRuvNoSTilda(V,S,f,pseudoInverseMat,Ruvs,explookup,termorder,vars_reversed)
#    n = nvars(parent(f)) - 1
#    d = total_degree(f)
#    R = coefficient_ring(parent(f))
#    MS1 = matrix_space(R, binomial(n*d,n*d-n), binomial(n*d,n*d-n))
#    if haskey(Ruvs, V)
#        return get(Ruvs, V, 0)
#    else
#        println("New key: $V")
#    end
#    ev1 = gen_exp_vec(n+1,n*d-n,termorder)
#    ev2 = gen_exp_vec(n+1,n*d-n+d-length(S),termorder)
#    ev3 = gen_exp_vec(n+1,n*d-n-length(S)+1,termorder)
#    temp = Vector{Int64}(undef, n+1)
#    MS2 = matrix_space(R, length(ev2),1)
#    result = Vector{typeof(MS1())}(undef, n+2)
#    for i in 1:n+2
#        result[i] = MS1(0)
#    end
#    for i in 1:length(ev1)
#        mon = Vector{Int64}(undef, n+1)
#        for m in 1:(n+1)
#            mon[m] = ev1[i][m] + V[m] - 1
#        end
#        gVec = MS2()
#        for j in 1:length(ev2)
#            if ev2[j] == mon
#                gVec[j] = R(1)
#            else
#                gVec[j] = R(0)
#            end
#        end
#        gJS = pseudoInverseMat*gVec
#        #println("After LingAlg problem: $gJS")
#        for j in 1:(n+1)
#            for k in 1:length(ev3)
#                for m in 1:(n+1)
#                    if m == n+1-j+1
#                        temp[m] = ev3[k][m]
#                    else
#                        temp[m] = ev3[k][m] + 1
#                    end
#                end
#                #print("ev1[l]: $((ev1[l],typeof(ev1[l])));")
#                #print("ev3[k]: $((ev3[k],typeof(ev3[k])));") 
#                #println(" $(ev1[l] == ev3[k])")
#                l = get(explookup,temp,0)
#                #result[j+1][l,i] = gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
#                #result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((j-1)*(length(gJS)/(n+1))+1):Int(j*(length(gJS)/(n+1))),:][k]
#                if vars_reversed == true
#                    result[j+1][l,i] = gJS[Int((j-1)*(length(gJS)/(n+1))+1)+k-1,1]
#                    result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((j-1)*(length(gJS)/(n+1))+1)+k-1,1]
#                else
#                    result[j+1][l,i] = gJS[Int((n+1-j)*(length(gJS)/(n+1))+1)+k-1,1]
#                    result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((n+1-j)*(length(gJS)/(n+1))+1)+k-1,1]
#                end
#                #println("$(result[j+1][l,i]) in $(j+1) && $(result[1][l,i]) in $(1)")
#            end
#        end
#    end
#    get!(Ruvs, V, result)
#    return result
#end

#"""
#Computes the Ruv matrix with the u being variables, stores this as n+2 matrices
#"""
#function computeRPoly_LAOneVar1(V,S,f,pseudoInverseMat,Ruvs,cache,termorder)
#    if haskey(Ruvs, V)
#        return get(Ruvs, V, 0)
#    end
#    n = nvars(parent(f)) - 1
#    d = total_degree(f)
#    R = coefficient_ring(parent(f))
#    URing, UVars = polynomial_ring(R, ["u$i" for i in 0:n])
#    PURing, Vars = polynomial_ring(URing, ["x$i" for i in 0:n])
#    #=
#    yV = []
#    for i in axes(V,1)
#        push!(yV, UVars[1]*V[i])
#    end
#    U = UVars[2:(n+2)] + yV
#    =#
#    ev = cache[n*d-n]#gen_exp_vec(n+1,n*d-n,termorder)
#    monomials = gen_mon(ev,URing,PURing)
#    #monomials = gen_mon(ev,R,parent(f))
#    reductions = []
#    for m in monomials
#        push!(reductions, reduce_LA(UVars,V,S,f,pseudoInverseMat,[m,1],PURing,termorder)[1])
#    end
#    polyMatrix = Matrix(transpose(convert_p_to_m(reductions,ev)))
#    matSpace = matrix_space(R,nrows(polyMatrix),ncols(polyMatrix))
#    matrices = []
#    for k in 0:(n+1)
#        tempMat = matSpace()
#        tempExpVec = zeros(Int,n+1)
#        if k >= 1
#            tempExpVec[n+1-k+1] = 1
#        end
#        for i in 1:nrows(polyMatrix)
#            for j in 1:ncols(polyMatrix)
#                tempMat[i,j] = coeff(polyMatrix[i,j], tempExpVec)
#            end
#        end
#        push!(matrices, tempMat)
#    end
#    get!(Ruvs, V, matrices)
#    return matrices
#end

"""
    computeRPoly_LAOneVar2!(matrices, U)

Takes a list of n+2 matrices and ouputs a list two matrices [A,B] corresponding to R_{(x0,...,xn)+yv, v} = Ay + B

INPUTS: 
* "matrices" -- list, output of computeRPoly_LAOneVar1
* "U" -- vector, U =(x0, ..., xn)
* "R" -- ring, base ring of f
"""
function computeRPoly_LAOneVar2!(A,B,matrices, U, V, R,temp)
    zero!(B)
    add!(B,B,matrices[1])
    #B = matrices[1]

    #@assert size(B) == size(matrices[1])
    #for i in 1:size(B,2)
    #    for j in 1:size(B,1)
    #        B[j,i] = matrices[1][j,i]
    #    end
    #end

#    println(typeof(matrices[1][1,1]))

    zero!(A)


    #A = parent(A)()
    
    #R = base_ring(A)
    #nrows = size(A,2)
    #ncols = size(A,1)

    #zz = zero(base_ring(A))
    #for i in 1:nrows
    #    for j in 1:ncols
    #        B[i,j] = matrices[1][i,j]
    #        A[i,j] = zz
    #    end
    #end

    br = base_ring(A)

    #ui(i) = UInt(lift(ZZ,br(i)))
    ui(i) = 0 ≤ i ? UInt(i) : UInt(i + characteristic(br))

    for k in 2:(length(matrices))

        #B = B + matrices[k] * U[k-1]
        #A = A + matrices[k] * V[k-1]
        #println(test)

        my_mul!(temp,matrices[k],U[k-1])
        add!(B,B,temp)
        my_mul!(temp,matrices[k],V[k-1])
        add!(A,A,temp)

        #add!(B,B,matrices[k] * ui(U[k-1]))
        #add!(A,A,matrices[k] * ui(V[k-1]))


        #my_addmul!(B,B,ui(U[k-1]),matrices[k])
        #my_addmul!(A,A,ui(V[k-1]),matrices[k])

        #@assert BB == B
        #@assert AA == A

        #println(B)
        #error()
        #B += matrices[k] * U[k-1]
        #A += matrices[k] * V[k-1]

        

        #ukm1 = R(U[k-1])
        #vkm1 = R(V[k-1])

        #for i in 1:nrows
        #    for j in 1:ncols
        #        B[i,j] = B[i,j] + matrices[k][i,j] * ukm1
        #        A[i,j] = A[i,j] + matrices[k][i,j] * vkm1
        #    end
        #end
    end 

    return (A, B)
end

function printMat(M)
    println()
    for i in axes(M,1)
        for j in axes(M,2)
            print(M[i,j])
            print(" ")
        end
        println()
    end
    println()
end

#----------------------------------
#=
"""
    computeD(N, m)

Returns a list of length N where D_{j,m} = sum_{i=j}^{N-1} (-1)^{i+j}binom{-m}{i}binom{i}{j}

INPUTS: 
* "N" -- integer
* "m" -- integer
"""
function computeD(N, m)
    D = zeros(Int,N)
    for j in 0:(N-1)
        D[j+1] = sum((-1)^(i+j)*binomial(-m,i)*binomial(i,j) for i in j:(N-1))
    end
    return D
end


"""
    applyFrobeniusToMon(n,d,f,N,p,beta,m,R,PR)

Computes the power series expansion of p^{m-n-1}sigma(x^{beta}Omega/f^m) 
using formula (1.10) in Costa's thesis


INPUTS: 
* "n" -- integer, dimension of ambient projective space
* "d" -- integer, degree of the hypersurface f
* "f" -- polynomial, defining homogeneous equation of the hypersurface lifted to characteristic 0
* "N" -- integer, series precision
* "p" -- integer, a prime number that is the characteristic of the base field of the hypersurface
* "beta" -- vector, representing the exponents in the monomial of the basis element
* "m" -- integer, pole order of the basis element 
* "R" -- ring, precision ring 
* "PR" -- ring, polynomial ring with coefficients in R 
"""
function applyFrobeniusToMon(n, d, f, N, p, beta, m, R, PR)
    println("N=$N, m=$m")
    Factorial = factorial(big(p * (N + m - 1) - 1))
    o = ones(Int64, n+1)
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), o)
    X1 = finish(B)
    D = computeD(N,m)
    result = []
    for j in 0:(N-1)
        e = j + m
        factorial_e = R(ZZ(Factorial/factorial(big(p * e - 1))))
        println("e=$e,factorial_e=$factorial_e")
        ev = gen_exp_vec(n+1,d*j)
        fj = f^j
        sum = 0
        for alpha in ev
            B = MPolyBuildCtx(PR)
            push_term!(B, R(1), p * (beta + alpha + o))
            monomial = div(finish(B), X1)
            sum = sum + R(factorial_e * (D[j+1] * (coeff(fj,alpha)^p))) * monomial
            #println(typeof((D[j+1]*(coeff(map_coefficients(lift,fj),alpha)^p))*monomial))
        end
        push!(result, [sum, p*(m+j)])
    end
    return result
end
=#

#=
function applyFrobenius(n,d,f,N,p,poly,R,PR)
    t = getTerms(poly)
    temp = []
    for i in t
        ev = exponent_vector(i[1],1)
        push!(temp, applyFrobeniusToMon(n,d,f,N,p,ev,poly[2],R,PR))
    end
    return temp
end

"""
    applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)

Applies the frobenius to all the elements of Basis

INPUTS: 
* "Basis" -- array of basis elmenets
* "n" -- integer, dimension of ambient projective space
* "d" -- integer, degree of f
* "f" -- polynomial, defining equation of hypersurface (lifted version)
* "N" -- series precision
* "p" -- the prime
* "R" -- basering(parent(f))
* "PR" -- parent(f)
"""
function applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)
    result = []
    for b in Basis
        Fmon = applyFrobeniusToMon(n,d,f,N,p,exponent_vector(b[1],1),b[2],R,PR)
        #println(Fmon)
        push!(result, Fmon)
    end
    return result
end
=#

#=
function reduce(U,V,S,n,d,g,parts,ev,R,PR,Vars)
    SC = []
    gensJS = copy(parts)
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), V)
    XV = finish(B)
    for i in 0:n
        if i in S
            #push!(gensJS,parts[i+1])
            gensJS[i+1] = parts[i+1]
        else
            #push!(gensJS,Vars[i+1]*parts[i+1])
            gensJS[i+1] = Vars[i+1]*parts[i+1]
            #parts[i+1] = Vars[i+1]*parts[i+1]
            push!(SC,i)
        end
    end
    # get gi's using pseudoinverse
    XS =  prod(PR(Vars[i+1]) for i in S; init = PR(1))
    gc, t = reduce_with_quotients(div(XV*g[1],XS),gensJS)
    gcpartials = [ derivative(gc[1], i) for i in 1:(n+1) ]
    return [sum(PR(U[i+1])*XS*gc[i+1] + div(XS,Vars[i+1])*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]

end
=#

#=
function reduceToBasis(U,V,S,n,d,g,parts,ev,R,PR,Vars)
    while g[2] > n
        g = reduce(U,V,S,n,d,g,parts,ev,R,PR,Vars)
    end
    return g
end
=#

#=
function evaluateRUV(RUV,U,ResRing)
    result = copy(RUV)
    for i in axes(RUV,1)
        for j in axes(RUV,2)
            result[i,j] = evaluate(change_base_ring(ResRing,map_coefficients(lift,RUV[i,j])),U)
        end
    end
    return result
end
=#

#=
function reducechain(I,n,d,m,S,parts,R,PR,ResRing)
    chain = 0
    ev = gen_exp_vec(n+1,n*d-n)
    gVec = chooseV(I,n*d - n)
    I = I - gVec
    Vs = []
    RUVs = []
    while m > n
        V = chooseV(I,d)
        U = I - V
        l = findall(x->x==V, Vs)
        if length(l) > 0
            RPoly = RUVs[l[1]]
        else
            RPoly = computeRPoly(V,S,n,d,parts,R,PR)
            push!(Vs,V)
            push!(RUVs,RPoly)
        end
        RNums = evaluateRUV(RPoly,U,ResRing)
        if chain == 0
            chain = RNums
        else
            chain = RNums*chain
        end
        m = m - 1
        I = U
    end
    return [chain, I]
end
=#

#=
function reducechain_LA1(I,gCoeff,l,n,d,m,S,f,pseudoInverseMat,R,PR)
    #chain = 0
    gVec = chooseV(I,n*d - n)
    ev = gen_exp_vec(n+1,n*d-n)
    gMat = zeros(R,length(ev))
    for j in axes(gMat,1)
        if gVec == ev[j]
            gMat[j] = gCoeff
            break
        end
    end
    I = I - gVec
    #Vs = []
    #RUVs = []
    h = 0
    while m > l
        V = chooseV(I,d)
        #U = I - V
        K = 0
        mins = I
        while true
            if m - K <= l
                break
            end
            temp = mins - V
            isLessThanZero = false
            for j in temp
                if j < 0
                    isLessThanZero = true
                    break
                end
            end
            if isLessThanZero == true
                break
            end
            mins = temp
            K = K+1
        end
        #=
        l = findall(x->x==V, Vs)
        if length(l) > 0
            RPoly = RUVs[l[1]]
        else
        =#
        A,B = computeRPoly_LAOneVar(V,mins,S,n,d,f,pseudoInverseMat,R,PR)
        #push!(Vs,V)
        #push!(RUVs,RPoly)
        #RNums = evaluateRUV(RPoly,U,R)
        MK = A + B*K
        MK1 = A + B*(K-1)
        h = MK*gMat
        if K >= 2
            h = MK1*h
            A1 = MK - MK1
            j = 2
            while K-j >= 0
                MKj = MK - A1*j
                h = MKj*h
                j = j + 1
            end
        end
        #=
        if chain == 0
            chain = RNums
        else
            chain = RNums*chain
        end
        =#
        m = m - K
        I = mins
    end
    return [h, I]
end
=#
#-----------------------------
#=
function reducepoly(poly,n,d,S,parts,R,PR,ResRing)
    t = getTerms(poly)
    result = 0
    for i in t
        o = ones(Int,length(exponent_vector(i[1],1)))
        I = exponent_vector(i[1],1) + o
        gVec = chooseV(I,n*d - n)
        ev = gen_exp_vec(n+1,n*d-n)
        gMat = zeros(Int,length(ev))
        for j in axes(gMat,1)
            if gVec == ev[j]
                gMat[j] = coeff(map_coefficients(lift,i[1]),1)
                break
            end
        end
        RChain = reducechain(I,n,d,i[2],S,parts,R,PR,ResRing)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), RChain[2])
        XU = finish(B)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), o)
        XS = finish(B)
        gReduction = div(XU*convert_m_to_p(transpose(RChain[1]*gMat),ev,R,PR)[1],XS)
        result = result + gReduction
    end
    return [result,n]
end
=#

#=
function reducepoly_LA1(poly,l,n,d,S,f,pseudoInverseMat,R,PR)
    t = getTerms(poly)
    result = 0
    for i in t
        o = ones(Int,length(exponent_vector(i[1],1)))
        I = exponent_vector(i[1],1) + o
        gCoeff = coeff(i[1],1)
        RChain = reducechain_LA1(I,gCoeff,l,n,d,i[2],S,f,pseudoInverseMat,R,PR)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), Int.(RChain[2]))
        XU = finish(B)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), o)
        XS = finish(B)
        ev = gen_exp_vec(n+1,n*d-n)
        gReduction = div(XU*convert_m_to_p(transpose(RChain[1]),ev,R,PR)[1],XS)
        result = result + gReduction
    end
    return [result,l]
end
=#


#=
function computeR(u,v,s,n,d,R,PR,vars)
    ev = gen_exp_vec(n+1,d*n-n)
    monomials = gen_mon(ev,R,PR)
    reductions = []
    for m in monomials
        push!(reductions,reduce(u,v,s,n,d,m,ev,R,PR,vars))
    end
    return transpose(convert_p_to_m(reductions,ev))
end
=#
#=
function computeRPoly(V,S,n,d,parts,R,PR)
    URing, UVars = polynomial_ring(R, ["u$i" for i in 0:n])
    PURing, Vars = polynomial_ring(URing, ["x$i" for i in 0:n])
    #this is temporary
    polynomial = Vars[1]^5 + Vars[2]^5 + Vars[3]^5 + Vars[1]*(Vars[2]^3)*Vars[3]
    parts = [ derivative(polynomial, i) for i in 1:(n+1) ]
    #
    ev = gen_exp_vec(n+1,d*n-n)
    monomials = gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, reduce(UVars,V,S,n,d,[m,1],parts,[],URing,PURing,Vars)[1])
    end
    return Matrix(transpose(convert_p_to_m(reductions,ev)))
end
=#

#LA test ------------------------
#=
function computeRPoly_LA(V,S,n,d,f,pseudoInverseMat,R,PR)
    URing, UVars = polynomial_ring(R, ["u$i" for i in 0:n])
    PURing, Vars = polynomial_ring(URing, ["x$i" for i in 0:n])
    ev = gen_exp_vec(n+1,d*n-n)
    monomials = gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, reduce_LA(UVars,V,S,n,d,f,pseudoInverseMat,[m,1],[],URing,PURing,Vars)[1])
    end
    return Matrix(transpose(convert_p_to_m(reductions,ev)))
end
=#
#-----------------------------------
#=
function reducetransform(FT,n,d,p,N,S,parts,R,PR)
    result = 0
    for i in FT
        for j in i
            s = N + j[2] - 1
            ResRing = residue_ring(ZZ,Int1024(p)^s)
            result = reducepoly(j,n,d,S,parts,R,PR,ResRing)[1]
        end
    end
    return [result,n]
end
=#

#=
function computeT(Basis,f,n,d,R,PR)
    ev = gen_exp_vec(n+1,d*n-n-1)
    mons = gen_mon(ev,R,PR)
    T = []
    for m in mons
        temp = []
        for i in 0:(n-1)
            if i == 0
            else
                mterms = terms(m)
                sum = 0
                for t in mterms
                    sum = sum + StandardReduction.stdRed_step(f,t,n-i+1,1)[1]
                end
                m = inv(R(n-i))*sum
            end
            for j in 0:(length(Basis[n-i])-1)
                c = coeff(PR(m),exponent_vector(Basis[n-i][length(Basis[n-i])-j],1))
                push!(temp, c)
                m = m - c*Basis[n-i][length(Basis[n-i])-j]
            end
        end
        push!(T, transpose(temp))
    end
    return transpose(vcat(T...))
end
=#

#=
"""
    liftCoefficients(R, PR, f)

Lifts the coefficeints of f to the ring R.

Works by lifting coefficients to ZZ and then 
converting elements of ZZ to elements of R.
Thus, this method is mostly useful for 
lifting something mod p^m to p^n,
for m < n.

INPUTS: 
* "f" -- the polynomial to be lifted
* "R" -- the ring for the coefficients to end up in
* "PR" -- the polynomial ring (over R) for the result to end up in 
"""
function liftCoefficients(R, PR, f, positiveLift=true)
    t = terms(f)
    sum = 0 
    for i in t
        ev = exponent_vector(i,1)
        c = coeff(i,1)
        B = MPolyBuildCtx(PR)
        charBaseField = characteristic(parent(f))
        if positiveLift && (lift(ZZ,c) > div(charBaseField, 2))
            push_term!(B, R(lift(ZZ,c)-charBaseField), ev)
        else
            push_term!(B, R(lift(ZZ,c)), ev)
        end
        sum = sum + finish(B)
    end
    return sum
end
=#
    
#end

#=
include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("jl")
include("jl")
include("SmallestSubsetSmooth.jl")
include("ZetaFunction.jl")
n = 2
d = 3
p = 7
R = GF(p,1)
PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
x0,x1,x2 = Vars
f = x1^2*x2 - x0^3 - x0*x2^2 - x2^3
Test = ZetaFunction.computeAll(n,d,f,7,p,R,PR,Vars)
=#

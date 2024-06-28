module ControlledReduction

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("Utils.jl")
include("PolynomialWithPole.jl")
#include("SmallestSubsetSmooth.jl")

"""
    reduce_LA(U,V,S,n,d,f,pseudoInverseMat,g,ev,R,PR,Vars)

applies reduction formula from Prop 1.15 in Costa's thesis to 
basis elements of Homog(dn-d), returns them as polynomials
"""
function reduce_LA(U,V,S,n,d,f,pseudoInverseMat,g,ev,R,PR,Vars)
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
    gVec = Utils.convert_p_to_m([div(XV*(g[1]),XS)],Utils.gen_exp_vec(n+1,n*d-n+d-length(S),:invlex))
    gJS = pseudoInverseMat*transpose(gVec)
    gc = []
    for i in 1:(n+1)
        push!(gc, Utils.convert_m_to_p(transpose(gJS[Int((i-1)*(length(gJS)/(n+1))+1):Int(i*(length(gJS)/(n+1)))]),Utils.gen_exp_vec(n+1,n*d-n-length(S)+1,:invlex),R,PR)[1])
    end
    gc = reverse(gc)
    gcpartials = [ derivative(gc[i], i) for i in 1:(n+1) ]
    
    gcpartials = reverse(gcpartials) # TODO: make this an option, this is the way it is in Costa's code, 

    #return [sum(PR(U[i+1])*XS*gc[i+1] + div(XS,Vars[i+1])*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]
    return [sum(PR(U[i+1])*div(XS,Vars[i+1])*gc[i+1] + XS*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]

end

"""
    chooseV(I, d)

choose direction of reduction in the same way as Costa's code

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f 
"""
function chooseV(I, d)
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
    I = reverse(I)

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

    V = reverse(V)

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
    reducechain_LA(u,g,n,d,p,m,S,f,pseudoInverseMat,R,PR)

takes single monomial in frobenius and reduces to pole order n, currently only does one chunk of reduction


if the reduction hits the end, returns u as the "true" value, otherwise returns it in Costa's format
(i.e. entries will be multiplies of p in Costa's format)
"""
function reducechain_LA(u,g,n,d,p,m,S,f,pseudoInverseMat,R,PR)
    
    #println(pseudoInverseMat)

    # @assert pseudoInverseMat ==
    pseudoInverseMat = 
    [155 0 0 0 0 11 0 0 0 221 0 0 310 0 22;
    0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
    225 0 0 0 0 155 0 0 0 11 0 0 221 0 310;
    0 0 0 0 0 0 0 0 0 0 0 0 0 114 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 114;
    0 0 0 1 0 0 172 0 0 0 0 0 0 0 0;
    59 0 0 0 1 94 0 172 0 166 0 0 61 0 188;
    0 0 0 0 0 0 0 0 172 0 0 0 0 286 0;
    0 57 0 173 0 0 0 0 0 0 172 0 0 114 0;
    83 0 57 0 173 59 0 0 0 94 0 172 166 0 61;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    114 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 114 0 0 0 0 0 0 0 0 0 0 0 0 0;
    166 0 114 0 0 118 0 0 0 188 0 0 332 0 236;
    11 0 0 0 0 221 0 0 0 310 0 0 22 0 214;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;] #"seems to have the wrong pseudo-inverse"
    #I = [28,7,28]
    #gCoeff = R(2)
    
    I = u

    #TODO?
    #I = reverse(I) # parity issue due to Costa's code being reverse from ours

    gMat = g
    #chain = 0
    println("This is I: $I")
    J = copy(I)

    #TODO?
    #V = rev_chooseV(Array{Int}(divexact.(I,p)),d)
    V = chooseV(Array{Int}(divexact.(I,p)),d)



    gVec = I - rev_tweak(I,n*d-n)
    ev = Utils.gen_exp_vec(n+1,n*d-n,:invlex)
    #gMat = zeros(R,length(ev))
    #for j in axes(gMat,1)
    #    if gVec == ev[j]
    #        gMat[j] = gCoeff
    #        break
    #    end
    #end
    I = I - gVec
    if m - n < p
        nend = m - n
    else
        nend = p
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

    println("Getting reduction matrix for V = $V")

    A,B = computeRPoly_LAOneVar(V,I - (nend-(d*n-n))*V,S,n,d,f,pseudoInverseMat,R,PR)
    matrices = computeRPoly_LAOneVar1(V,S,n,d,f,pseudoInverseMat,R,PR)

    #for i in axes(matrices,1)
    #    printMat(matrices[i])
    #end
    
    if V == [1,1,1]
        #println("Using precomputed R_u,[1,1,1]")
    matrices_precomputed = [
                R[[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[169 0 0 0 1 249 0 172 0 177 0 0 282 0 155]
[0 114 0 0 0 0 0 0 1 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[94 0 57 0 173 280 0 0 0 61 0 172 188 0 160]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
],
                R[[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[155 0 0 0 0 11 0 0 0 221 0 0 310 0 22]
[0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[225 0 0 0 0 155 0 0 0 11 0 0 221 0 310]
[0 0 0 0 0 0 0 0 0 0 0 0 0 114 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 114]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
],
                R[[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 1 0 0 172 0 0 0 0 0 0 0 0]
[59 0 0 0 1 94 0 172 0 166 0 0 61 0 188]
[0 0 0 0 0 0 0 0 172 0 0 0 0 286 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 57 0 173 0 0 0 0 0 0 172 0 0 114 0]
[83 0 57 0 173 59 0 0 0 94 0 172 166 0 61]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
],
                R[[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[114 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 114 0 0 0 0 0 0 0 0 0 0 0 0 0]
[166 0 114 0 0 118 0 0 0 188 0 0 332 0 236]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[11 0 0 0 0 221 0 0 0 310 0 0 22 0 214]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
]
]
      #matrices = map(M -> R.(M),matrices)
      for i in 1:4
          @assert matrices[i] == matrices_precomputed[i] "the $i-th R_u,v matrix seems to have been computed wrong"
      end
    end 

    A1,B1 = computeRPoly_LAOneVar2(matrices,I - (nend-(d*n-n))*V,R)
    i = 1
    
    println("Before reduction chunk: $gMat")
    println("Before reduction chunk, I is $I")
    fastevaluation = false
    if fastevaluation
      gMat = finitediff_prodval_linear(B,A,nend-(dn-n),nend,gMat)
    else
      while i <= (nend-(d*n-n))
        gMat = (A+B*(nend-(d*n-n)-i))*gMat

        println("After step $i: $gMat")

        i = i+1
      end
    end
    # TODO: test how much of a difference the fast evaluation actually makes

    I = I - (nend-(d*n-n))*V
    println("After steps 1-$i, I is $I")
    i = i-1
    while i <= nend-1
        y = rev_tweak(J - i*V,d*n-n) - rev_tweak(J - (i+1)*V,d*n-n)
        println("Getting y direction reduction matrix for V = $(y)") 
        # there's some sort of parity issue between our code and edgar's
        A,B = computeRPoly_LAOneVar(y,rev_tweak(J - (i+1)*V,d*n-n) - y,S,n,d,f,pseudoInverseMat,R,PR)
        gMat = (A+B)*gMat
        println("After step $(i+1): $gMat")

        i = i+1
        I = I - y
        println("After step $(i+1), I is $I")
    end
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
    
    
    if nend == p
        newI = J - p*V
        @assert undo_rev_tweak(I,p) == newI

        return (newI, gMat)
    else
        return (I,gMat) # gives the "true" u
    end

end

"""
finitediff_prodeval_linear(a,b,start,stop,g)

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

"""
function finitediff_prodeval_linear(a,b,start,stop,g)
  if start == stop
    return (a .* stop .+ b) * g
  end

  Fk = a .* stop .+ b # Fk = F(k), here k = stop

  g = Fstop*g

  for k = stop-1:-1:start
    # right now, Fk = F(k+1)
    Fk = Fk .- a
    # now, Fk = F(k)
    g = Fk * g
  end

  g
end

"""
Returns the data used by costa's code given a polynomial term.

Note: this only works with a single term, so it should
only be used at the beginning of reduction
"""
function costadata_of_initial_term(term,n,d,p)

    R = base_ring(parent(term[1])) 
    i = term

    o = ones(Int,length(exponent_vector(i[1],1)))
    II = exponent_vector(i[1],1) + o # doh't use I since it's overloaded with the id matrix
    gCoeff = coeff(i[1],1)
    

    V = chooseV(Array{Int}(divexact.(II,p)),d)
    u = II - rev_tweak(II,n*d-n)
    ev = Utils.gen_exp_vec(n+1,n*d-n,:invlex)
    g = zeros(R,length(ev))
    for j in axes(g,1)
        if u == ev[j]
            g[j] = gCoeff
            break
        end
    end


    println("creation: u is type $(typeof(II))")
    (II,g)
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
        println("incorporation: u has type $(typeof(costadata[1]))")
        push!(costadata_arr,costadata)
    end
end

"""
Converts a Costa's data tuple to a polynomial with pole
after reduction.

Thus, the pole order is n.
"""
function poly_of_end_costadata(costadata,PR,p,d,n)
    (u,g_vec) = costadata
    vars = gens(PR)

    g = Utils.vector_to_polynomial(g_vec,n,d*n-n,PR,:invlex)

    # no need to do rev_tweak since reducechain_LA returns the "true" u
    # on the last run
    [prod(vars .^ u) * g, n]
end


"""
Converts an array of Costa's data tuples to an array of polynomials with pole
after reduction.

Thus, the pole order is always n.

"""
function poly_of_end_costadatas(costadatas,PR,p,d,n)
    res = PR(0)
    for costadata in costadatas
        res += poly_of_end_costadata(costadata,PR,p,d,n)[1]
    end
   
    [[res, n]]
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
#        ev = Utils.gen_exp_vec(n+1,d*n - n,:invlex)
#        gReduction = div(XU*Utils.convert_m_to_p(transpose(RChain[1]),ev,R,PR)[1],XS)
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
#        push!(result,[temp,n,Utils.Factorial(R(i[length(i)][2]),R(n))])
#    end
#    return result
#end

"""
Implements Costa's algorithm for controlled reduction,
sweeping down the terms of the series expansion by the pole order.
"""
function reducepoly_LA_descending(pol,n,d,p,S,f,pseudoInverseMat,R,PR)
    println(pol)

    i = pol
    highpoleorder = i[length(i)][2]

    # this is the omega from section 1.5.5 of Costa's thesis.
    ω = [] # this will be an array of costa data

    poleorder = highpoleorder
    while n < poleorder
        # this is an array of polynomials
        ωₑ = PolynomialWithPole.termsoforder(pol,poleorder)

        println("ωₑ: $ωₑ")

        for term in ωₑ
            println("term: $term")
            term_costadata = costadata_of_initial_term(term,n,d,p)
            println("term, in Costa's format: $term_costadata")
            #ω = ω + ωₑ
            incorporate_initial_term!(ω,term_costadata)
        end

        println("ω: $ω")
        #ω = reducepoly_LA(ω,n,d,p,S,f,pseudoInverseMat,R,PR)
        for i in eachindex(ω)
            #ω[i] = reducechain...
            println("u is type $(typeof(ω[i][1]))")
            ω[i] = reducechain_LA(ω[i]...,n,d,p,poleorder,S,f,pseudoInverseMat,R,PR)
        end

        poleorder = poleorder - p
    end

    poly_of_end_costadatas(ω,PR,p,d,n)
    println(poly_of_end_costadatas(ω,PR,p,d,n))
end

"""
trying to emulate Costa's controlled reduction, changes the order that polynomials are reduced, starts from highest pole order and accumulates the lower order poles as reduction proceeds

TODO: what exactly is big N?? Why isn't is used?
"""
function reducetransform_LA_descending(FT,n,d,p,N,S,f,pseudoInverseMat,R,PR)
    result = []
    for pol in FT
        reduction = reducepoly_LA_descending(pol,n,d,p,S,f,pseudoInverseMat,R,PR)
        push!(result, reduction)
    end
    return result
end
        
"""
computes reduction matrices
"""
function computeRPoly_LAOneVar(V,mins,S,n,d,f,pseudoInverseMat,R,PR)
    YRing, y = polynomial_ring(R, "y")
    PYRing, Vars = polynomial_ring(YRing, ["x$i" for i in 0:n])
    yV = []
    for i in axes(V,1)
        push!(yV, y*V[i])
    end
    UVars = mins + yV
    ev = Utils.gen_exp_vec(n+1,n*d-n,:invlex)
    monomials = Utils.gen_mon(ev,YRing,PYRing)
    reductions = []
    for m in monomials
        push!(reductions, reduce_LA(UVars,V,S,n,d,f,pseudoInverseMat,[m,1],[],YRing,PYRing,Vars)[1])
    end
    polyMatrix = Matrix(transpose(Utils.convert_p_to_m(reductions,ev)))
    matSpace = matrix_space(R,nrows(polyMatrix),ncols(polyMatrix))
    A = matSpace()
    B = matSpace()
    for i in 1:nrows(polyMatrix)
        for j in 1:ncols(polyMatrix)
            A[i,j] = coeff(polyMatrix[i,j],0)
            B[i,j] = coeff(polyMatrix[i,j],1)
        end
    end
    return [A,B]
end

"""
Computes the Ruv matrix with the u being variables, stores this as n+2 matrices
"""
function computeRPoly_LAOneVar1(V,S,n,d,f,pseudoInverseMat,R,PR)
    URing, UVars = polynomial_ring(R, ["u$i" for i in 0:n])
    PURing, Vars = polynomial_ring(URing, ["x$i" for i in 0:n])
    #=
    yV = []
    for i in axes(V,1)
        push!(yV, UVars[1]*V[i])
    end
    U = UVars[2:(n+2)] + yV
    =#
    ev = Utils.gen_exp_vec(n+1,n*d-n,:invlex)
    monomials = Utils.gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, reduce_LA(UVars,V,S,n,d,f,pseudoInverseMat,[m,1],[],URing,PURing,Vars)[1])
    end
    polyMatrix = Matrix(transpose(Utils.convert_p_to_m(reductions,ev)))
    matSpace = matrix_space(R,nrows(polyMatrix),ncols(polyMatrix))
    matrices = []
    for k in 0:(n+1)
        tempMat = matSpace()
        tempExpVec = zeros(Int,n+1)
        if k >= 1
            tempExpVec[n+1-k+1] = 1
        end
        for i in 1:nrows(polyMatrix)
            for j in 1:ncols(polyMatrix)
                tempMat[i,j] = coeff(polyMatrix[i,j], tempExpVec)
            end
        end
        push!(matrices, tempMat)
    end
    return matrices
end

"""
    computeRPoly_LAOneVar2(matrices, U)

Takes a list of n+2 matrices and ouputs a list two matrices [A,B] corresponding to R_{(x0,...,xn)+yv, v} = Ay + B

INPUTS: 
* "matrices" -- list, output of computeRPoly_LAOneVar1
* "U" -- vector, U =(x0, ..., xn)
* "R" -- ring, base ring of f
"""
function computeRPoly_LAOneVar2(matrices, U, R)
    B = matrices[1]
    matSpace = matrix_space(R, nrows(B), ncols(B))
    A = matSpace()
    for k in 2:(length(matrices))
        A += matrices[k] * U[k-1]
    end 

    return [A, B]
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
        ev = Utils.gen_exp_vec(n+1,d*j)
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
    ev = Utils.gen_exp_vec(n+1,n*d-n)
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
    ev = Utils.gen_exp_vec(n+1,n*d-n)
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
        ev = Utils.gen_exp_vec(n+1,n*d-n)
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
        gReduction = div(XU*Utils.convert_m_to_p(transpose(RChain[1]*gMat),ev,R,PR)[1],XS)
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
        ev = Utils.gen_exp_vec(n+1,n*d-n)
        gReduction = div(XU*Utils.convert_m_to_p(transpose(RChain[1]),ev,R,PR)[1],XS)
        result = result + gReduction
    end
    return [result,l]
end
=#


#=
function computeR(u,v,s,n,d,R,PR,vars)
    ev = Utils.gen_exp_vec(n+1,d*n-n)
    monomials = Utils.gen_mon(ev,R,PR)
    reductions = []
    for m in monomials
        push!(reductions,reduce(u,v,s,n,d,m,ev,R,PR,vars))
    end
    return transpose(Utils.convert_p_to_m(reductions,ev))
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
    ev = Utils.gen_exp_vec(n+1,d*n-n)
    monomials = Utils.gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, reduce(UVars,V,S,n,d,[m,1],parts,[],URing,PURing,Vars)[1])
    end
    return Matrix(transpose(Utils.convert_p_to_m(reductions,ev)))
end
=#

#LA test ------------------------
#=
function computeRPoly_LA(V,S,n,d,f,pseudoInverseMat,R,PR)
    URing, UVars = polynomial_ring(R, ["u$i" for i in 0:n])
    PURing, Vars = polynomial_ring(URing, ["x$i" for i in 0:n])
    ev = Utils.gen_exp_vec(n+1,d*n-n)
    monomials = Utils.gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, reduce_LA(UVars,V,S,n,d,f,pseudoInverseMat,[m,1],[],URing,PURing,Vars)[1])
    end
    return Matrix(transpose(Utils.convert_p_to_m(reductions,ev)))
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
    ev = Utils.gen_exp_vec(n+1,d*n-n-1)
    mons = Utils.gen_mon(ev,R,PR)
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
    
end

#=
include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("Utils.jl")
include("Utils.jl")
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

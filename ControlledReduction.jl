module ControlledReduction

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
#include("FindMonomialBasis.jl")
include("AutomatedScript.jl")
include("Utils.jl")
#include("SmallestSubsetSmooth.jl")
include("StandardReduction.jl")

function computeReduction(U,V,S,n,d,g,parts,ev,R,PR,Vars)
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
    # get gi's using psuedoinverse
    XS =  prod(PR(Vars[i+1]) for i in S; init = PR(1))
    gc, t = reduce_with_quotients(div(XV*g[1],XS),gensJS)
    gcpartials = [ derivative(gc[1], i) for i in 1:(n+1) ]
    return [sum(PR(U[i+1])*XS*gc[i+1] + div(XS,Vars[i+1])*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]

end

#LA test --------------------------------
function computeReductionLA(U,V,S,n,d,f,psuedoInverseMat,g,ev,R,PR,Vars)
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
    # get gi's using psuedoinverse
    XS =  prod(PR(Vars[i+1]) for i in S; init = PR(1))
    gVec = AutomatedScript.convert_p_to_m([div(XV*(g[1]),XS)],AutomatedScript.gen_exp_vec(n+1,n*d-n+d-length(S)))
    gJS = psuedoInverseMat*transpose(gVec)
    gc = []
    for i in 1:(n+1)
        push!(gc, AutomatedScript.convert_m_to_p(transpose(gJS[Int((i-1)*(length(gJS)/(n+1))+1):Int(i*(length(gJS)/(n+1)))]),AutomatedScript.gen_exp_vec(n+1,n*d-n-length(S)),R,PR)[1])
    end
    gcpartials = [ derivative(gc[i], i) for i in 1:(n+1) ]
    return [sum(PR(U[i+1])*XS*gc[i+1] + div(XS,Vars[i+1])*gcpartials[i+1] for i in S; init = PR(0)) + XS*sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]

end
#----------------------------------------

function chooseV(I,d)
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
    return V
end

function reduceToBasis(U,V,S,n,d,g,parts,ev,R,PR,Vars)
    while g[2] > n
        g = computeReduction(U,V,S,n,d,g,parts,ev,R,PR,Vars)
    end
    return g
end

function evaluateRUV(RUV,U,ResRing)
    result = copy(RUV)
    for i in axes(RUV,1)
        for j in axes(RUV,2)
            result[i,j] = evaluate(change_base_ring(ResRing,map_coefficients(lift,RUV[i,j])),U)
        end
    end
    return result
end

function computeReductionChain(I,n,d,m,S,parts,R,PR,ResRing)
    chain = 0
    ev = AutomatedScript.gen_exp_vec(n+1,n*d-n)
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

#LA test --------------------
function computeReductionChainLA(I,gCoeff,n,d,m,S,f,psuedoInverseMat,R,PR)
    #chain = 0
    gVec = chooseV(I,n*d - n)
    ev = AutomatedScript.gen_exp_vec(n+1,n*d-n)
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
    while m > n
        V = chooseV(I,d)
        #U = I - V
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
        #=
        l = findall(x->x==V, Vs)
        if length(l) > 0
            RPoly = RUVs[l[1]]
        else
        =#
        A,B = computeRPolyLAOneVar(V,mins,S,n,d,f,psuedoInverseMat,R,PR)
        #push!(Vs,V)
        #push!(RUVs,RPoly)
        #RNums = evaluateRUV(RPoly,U,R)
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

function computeReductionChainLA1(I,gCoeff,l,n,d,m,S,f,psuedoInverseMat,R,PR)
    #chain = 0
    gVec = chooseV(I,n*d - n)
    ev = AutomatedScript.gen_exp_vec(n+1,n*d-n)
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
        A,B = computeRPolyLAOneVar(V,mins,S,n,d,f,psuedoInverseMat,R,PR)
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
#-----------------------------

function computeReductionOfPoly(poly,n,d,S,parts,R,PR,ResRing)
    t = getTerms(poly)
    result = 0
    for i in t
        o = ones(Int,length(exponent_vector(i[1],1)))
        I = exponent_vector(i[1],1) + o
        gVec = chooseV(I,n*d - n)
        ev = AutomatedScript.gen_exp_vec(n+1,n*d-n)
        gMat = zeros(Int,length(ev))
        for j in axes(gMat,1)
            if gVec == ev[j]
                gMat[j] = coeff(map_coefficients(lift,i[1]),1)
                break
            end
        end
        RChain = computeReductionChain(I,n,d,i[2],S,parts,R,PR,ResRing)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), RChain[2])
        XU = finish(B)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), o)
        XS = finish(B)
        gReduction = div(XU*AutomatedScript.convert_m_to_p(transpose(RChain[1]*gMat),ev,R,PR)[1],XS)
        result = result + gReduction
    end
    return [result,n]
end

#LA test ------------------------
function computeReductionOfPolyLA(poly,n,d,S,f,psuedoInverseMat,R,PR)
    t = getTerms(poly)
    result = 0
    for i in t
        o = ones(Int,length(exponent_vector(i[1],1)))
        I = exponent_vector(i[1],1) + o
        gCoeff = coeff(i[1],1)
        RChain = computeReductionChainLA(I,gCoeff,n,d,i[2],S,f,psuedoInverseMat,R,PR)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), Int.(RChain[2]))
        XU = finish(B)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), o)
        XS = finish(B)
        ev = AutomatedScript.gen_exp_vec(n+1,d*n - n)
        gReduction = div(XU*AutomatedScript.convert_m_to_p(transpose(RChain[1]),ev,R,PR)[1],XS)
        result = result + gReduction
    end
    return [result,n]
end

function computeReductionOfPolyLA1(poly,l,n,d,S,f,psuedoInverseMat,R,PR)
    t = getTerms(poly)
    result = 0
    for i in t
        o = ones(Int,length(exponent_vector(i[1],1)))
        I = exponent_vector(i[1],1) + o
        gCoeff = coeff(i[1],1)
        RChain = computeReductionChainLA1(I,gCoeff,l,n,d,i[2],S,f,psuedoInverseMat,R,PR)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), Int.(RChain[2]))
        XU = finish(B)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(1), o)
        XS = finish(B)
        ev = AutomatedScript.gen_exp_vec(n+1,n*d-n)
        gReduction = div(XU*AutomatedScript.convert_m_to_p(transpose(RChain[1]),ev,R,PR)[1],XS)
        result = result + gReduction
    end
    return [result,l]
end
#-----------------------------------

function computeReductionOfTransform(FT,n,d,p,N,S,parts,R,PR)
    result = 0
    for i in FT
        for j in i
            s = N + j[2] - 1
            ResRing = residue_ring(ZZ,Int1024(p)^s)
            result = computeReductionOfPoly(j,n,d,S,parts,R,PR,ResRing)[1]
        end
    end
    return [result,n]
end

#LA test ---------------------------
function computeReductionOfTransformLA(FT,n,d,p,N,S,f,psuedoInverseMat,R,PR)
    result = []
    for i in FT
        temp = 0
        for j in i
            t = computeReductionOfPolyLA([Factorial(R(i[length(i)][2]),R(j[2]))p^(j[2]-n-1)*(j[1]),j[2]],n,d,S,f,psuedoInverseMat,R,PR)[1]
            temp = temp + t
        end
        #push!(result,[temp,n,Factorial(Int1024(i[length(i)][2]),n)])
        push!(result,[temp,n,Factorial(R(i[length(i)][2]),R(n))])
    end
    return result
end

function computeReductionOfTransformLA1(FT,n,d,p,N,S,f,psuedoInverseMat,R,PR)
    result = []
    for i in FT
        omega = [PR(0),i[length(i)][2]]
        m = Int(i[1][2]/p)
        for j in 1:(m+N-1)
            e = m+N-j
            omegae = Factorial(R(p*(m+N-1)-1),R(p*e-1))*i[e][1]
            omega[1] = omega[1] + omegae
            l = max(p*(e-1),n)
            omega = computeReductionOfPolyLA1([Factorial(p*(e-1),(l-1))*omega[1],omega[2]],l,n,d,S,f,psuedoInverseMat,R,PR)
        end
        push!(result,[omega[1],n,Factorial(R(p*(m+N-1)-1),R(1))])
    end
    return result
end
#-------------------------------------
        



#=
function computeR(u,v,s,n,d,R,PR,vars)
    ev = AutomatedScript.gen_exp_vec(n+1,d*n-n)
    monomials = AutomatedScript.gen_mon(ev,R,PR)
    reductions = []
    for m in monomials
        push!(reductions,computeReduction(u,v,s,n,d,m,ev,R,PR,vars))
    end
    return transpose(AutomatedScript.convert_p_to_m(reductions,ev))
end
=#

function computeRPoly(V,S,n,d,parts,R,PR)
    URing, UVars = polynomial_ring(R, ["u$i" for i in 0:n])
    PURing, Vars = polynomial_ring(URing, ["x$i" for i in 0:n])
    #this is temporary
    polynomial = Vars[1]^5 + Vars[2]^5 + Vars[3]^5 + Vars[1]*(Vars[2]^3)*Vars[3]
    parts = [ derivative(polynomial, i) for i in 1:(n+1) ]
    #
    ev = AutomatedScript.gen_exp_vec(n+1,d*n-n)
    monomials = AutomatedScript.gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, computeReduction(UVars,V,S,n,d,[m,1],parts,[],URing,PURing,Vars)[1])
    end
    return Matrix(transpose(AutomatedScript.convert_p_to_m(reductions,ev)))
end

#LA test ------------------------
function computeRPolyLA(V,S,n,d,f,psuedoInverseMat,R,PR)
    URing, UVars = polynomial_ring(R, ["u$i" for i in 0:n])
    PURing, Vars = polynomial_ring(URing, ["x$i" for i in 0:n])
    ev = AutomatedScript.gen_exp_vec(n+1,d*n-n)
    monomials = AutomatedScript.gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, computeReductionLA(UVars,V,S,n,d,f,psuedoInverseMat,[m,1],[],URing,PURing,Vars)[1])
    end
    return Matrix(transpose(AutomatedScript.convert_p_to_m(reductions,ev)))
end

function computeRPolyLAOneVar(V,mins,S,n,d,f,psuedoInverseMat,R,PR)
    YRing, y = polynomial_ring(R, "y")
    PYRing, Vars = polynomial_ring(YRing, ["x$i" for i in 0:n])
    yV = []
    for i in axes(V,1)
        push!(yV, y*V[i])
    end
    UVars = mins + yV
    ev = AutomatedScript.gen_exp_vec(n+1,n*d-n)
    monomials = AutomatedScript.gen_mon(ev,YRing,PYRing)
    reductions = []
    for m in monomials
        push!(reductions, computeReductionLA(UVars,V,S,n,d,f,psuedoInverseMat,[m,1],[],YRing,PYRing,Vars)[1])
    end
    polyMatrix = Matrix(transpose(AutomatedScript.convert_p_to_m(reductions,ev)))
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
#----------------------------------

"""
    computeD(N, m)

Returns a list of length N where D_{j,m} = \sum_{i=j}^{N-1} (-1)^{i+j}\binom{-m, i}\binom{i, j}

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

Computes the power series expansion of p^{m-n-1}\sigma(x^{\beta}\Omega/f^m) 
using formula (1.10) in Costa's thesis


INPUTS: 
* "n" -- integer, dimension of ambient projective space
* "d" -- integer, degree of the hypersurface f
* "f" -- polynomial, defining homogeneous equation of the hypersurface lifted to characteristic 0
* "N" -- integer, series precision
* "p" -- integer, a prime number that is the characteristic of the base field of the hypersurface
* "beta" -- vector, representing the exponents in the monomial of the basis element
* "m" -- integer, pole order of the basis element 
* "R" -- ring, basering(parent(f))
* "PR" -- ring, parent(f)
"""
function applyFrobeniusToMon(n, d, f, N, p, beta, m, R, PR)
    s = N + m -1 # never used? 
    o = ones(Int64, n+1)
    B = MPolyBuildCtx(PR)
    push_term!(B, R(1), o)
    X1 = finish(B)
    D = computeD(N,m)
    result = []
    for j in 0:(N-1)
        ev = AutomatedScript.gen_exp_vec(n+1,d*j)
        fj = f^j
        sum = 0
        for alpha in ev
            B = MPolyBuildCtx(PR)
            push_term!(B, R(1), p*(beta + alpha + o))
            monomial = div(finish(B),X1)
            sum = sum + p^(m-1)*((D[j+1]*(coeff(fj,alpha)^p))%(p^s))*monomial
            #println(typeof((D[j+1]*(coeff(map_coefficients(lift,fj),alpha)^p))*monomial))
        end
        push!(result, [sum, p*(m+j)])
    end
    return result
end

function getTerms(poly)
    t = terms(poly[1])
    result = []
    for i in t
        push!(result,[i,poly[2]])
    end
    return result
end

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

Basis - array of basis elmenets
n - number of variables minus 1
d - degree
f - polynomial which is the denominator of poles (lifted version)
N - series precision
p - the prime
R - basering(parent(f))
PR - parent(f)
"""
function applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)
    result = []
    for b in Basis
        Fmon = applyFrobeniusToMon(n,d,f,N,p,exponent_vector(b[1],1),b[2],R,PR)
        println(Fmon)
        push!(result, Fmon)
    end
    return result
end

function Factorial(x,y)
    fact = 1
    if x == 0
        return 1
    end
    while x != y
        fact = fact*x
        x = x - 1
    end
    return fact
end

function computeT(Basis,f,n,d,R,PR)
    ev = AutomatedScript.gen_exp_vec(n+1,d*n-n-1)
    mons = AutomatedScript.gen_mon(ev,R,PR)
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

"""
    liftCoefficients(R,PR,f)

Lifts the coefficeints of f to the ring R.

Works by lifting coefficients to ZZ and then 
converting elements of ZZ to elements of R.
Thus, this method is mostly useful for 
lifting something mod p^m to p^n,
for m < n.

f - the polynomial to be lifted
R - the ring for the coefficients to end up in
PR - the polynomial ring (over R) for the result
to end up in 
"""
function liftCoefficients(R,PR,f)
    t = terms(f)
    sum = 0 
    for i in t
        ev = exponent_vector(i,1)
        c = coeff(i,1)
        B = MPolyBuildCtx(PR)
        push_term!(B, R(lift(ZZ,c)), ev)
        sum = sum + finish(B)
    end
    return sum
end
    
end

#=
include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("AutomatedScript.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")
n = 2
d = 3
p = 7
R = GF(p,1)
PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
x,y,z = Vars
f = y^2*z - x^3 - x*z^2 - z^3
Test = ControlledReduction.computeAll(n,d,f,7,p,R,PR,Vars)
=#

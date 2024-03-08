module ControlledReduction

using Oscar
using BitIntegers
using LinearAlgebra
using Combinatorics

include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("AutomatedScript.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")
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
        push!(gc, AutomatedScript.convert_m_to_p(transpose(gJS[Int((i-1)*(length(gJS)/(n+1))+1):Int(i*(length(gJS)/(n+1)))]),AutomatedScript.gen_exp_vec(n+1,n*d-n),R,PR)[1])
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
        println(m)
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
        println(K)
        #=
        l = findall(x->x==V, Vs)
        if length(l) > 0
            RPoly = RUVs[l[1]]
        else
        =#
        RPoly = computeRPolyLAOneVar(V,mins,S,n,d,f,psuedoInverseMat,R,PR)
        #push!(Vs,V)
        #push!(RUVs,RPoly)
        #RNums = evaluateRUV(RPoly,U,R)
        MK = evaluateRUV(RPoly,K,R)
        MK1 = evaluateRUV(RPoly,K-1,R)
        h = MK1*MK*gMat
        A1 = MK - MK1
        j = 2
        while K-j >= 0
            MKj = MK - A1*j
            h = MKj*h
            println("multipled part of reduction chain")
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
    println("starting reduction of poly")
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
        ev = AutomatedScript.gen_exp_vec(n+1,d*n - n - 1)
        gReduction = div(XU*AutomatedScript.convert_m_to_p(transpose(RChain[1]),ev,R,PR)[1],XS)
        result = result + gReduction
    end
    return [result,n]
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
            temp = temp + computeReductionOfPolyLA([p^(j[2]-n-1)*(j[1]),j[2]],n,d,S,f,psuedoInverseMat,R,PR)[1]
            #R(inv(R(Factorial(j[2])))*Factorial(n))*
            println("computed reduction of part of basis element")
        end
        push!(result,[temp,n])
        println("computed reudction of basis element")
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
    ev = AutomatedScript.gen_exp_vec(n+1,d*n-n)
    monomials = AutomatedScript.gen_mon(ev,YRing,PYRing)
    reductions = []
    for m in monomials
        push!(reductions, computeReductionLA(UVars,V,S,n,d,f,psuedoInverseMat,[m,1],[],YRing,PYRing,Vars)[1])
    end
    return Matrix(transpose(AutomatedScript.convert_p_to_m(reductions,ev)))
end
#----------------------------------

function computeD(N,m)
    D = zeros(Int,N)
    for j in 0:(N-1)
        D[j+1] = sum((-1)^(i+j)*binomial(-m,i)*binomial(i,j) for i in j:(N-1))
    end
    return D
end

function applyFrobeniusToMon(n,d,f,N,p,beta,m,R,PR)
    s = N + m - 1
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
            sum = sum + p^(m-1)*(D[j+1]*(coeff(map_coefficients(lift,fj),alpha)^p))*monomial
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

function applyFrobeniusToBasis(Basis,n,d,f,N,p,R,PR)
    result = []
    for b in Basis
        push!(result, applyFrobeniusToMon(n,d,f,N,p,exponent_vector(b[1],1),b[2],R,PR))
        println("Applied Frob to Basis element")
    end
    return result
end

function Factorial(x)
    fact = 1
    while x != 0
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
                    sum = sum + StandardReduction.stdRed_step(f,exponent_vector(t,1),n-i+1)[1]
                end
                m = sum
            end
            for j in 0:(length(Basis[n-i])-1)
                println(n-i)
                println(length(Basis[n-i])-j)
                c = coeff(m,exponent_vector(Basis[n-i][length(Basis[n-i])-j],1))
                push!(temp, c)
                m = m - c*Basis[n-i][length(Basis[n-i])-j]
            end
        end
        push!(T, transpose(temp))
    end
    return vcat(T...)
end

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


function computeFrobeniusMatrix(n,d,f,precision,p,R,PR,vars)
    Nm = PrecisionEstimate.compute_precisions_each(p,precision,n)
    #N = max(Nm...)
    N = 6
    s = N + n - 1
    M = Int(precision + floor((p*s-1)/(p-1) + 1))
    PrecisionRing = PadicField(p,M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR)
    fLift = liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    T = computeT(BasisTLift,fLift,n,d,PrecisionRing,PrecisionRingPoly)
    println(BasisT)
    println(T)
    #S = SmallestSubsetSmooth.smallest_subset_s_smooth(fLift,n)
    Basis = []
    for i in 1:n
        for j in BasisT[i]
            push!(Basis,[liftCoefficients(PrecisionRing,PrecisionRingPoly,j),i])
        end
    end
    println(Basis)
    FBasis = applyFrobeniusToBasis(Basis,n,d,f,N,p,PrecisionRing,PrecisionRingPoly)
    println(FBasis[1])
    psuedoInverseMatTemp = CopiedFindMonomialBasis.psuedo_inverse_controlled(f,[],R,PR)
    psuedoInverseMat = zeros(PrecisionRing,nrows(psuedoInverseMatTemp),ncols(psuedoInverseMatTemp))
    for i in 1:nrows(psuedoInverseMat)
        for j in 1:ncols(psuedoInverseMat)
            psuedoInverseMat[i,j] = PrecisionRing(lift(psuedoInverseMatTemp[i,j]))
        end
    end
    denoms = []
    for i in FBasis
        for j in i
            push!(denoms,j[2])
        end
    end
    #Reductions = computeReductionOfTransformLA(FBasis,n,d,p,N,[],fLift,psuedoInverseMat,PrecisionRing,PrecisionRingPoly)
    return Reductions
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
n = 3
d = 4
p = 7
Fp = GF(p,1)

R = Fp
PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
x,y,z,w = Vars
f = x^4 + y^4 + z^4 + w^4
Test = ControlledReduction.computeFrobeniusMatrix(n,d,f,7,p,R,PR,Vars)
=#

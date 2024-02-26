module ControlledReduction

using Oscar
using BitIntegers
using LinearAlgebra
include("AutomatedScript.jl")

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
function computeReductionLA(U,V,S,n,d,f,g,ev,R,PR,Vars)
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
    partsMat = LinearAlgebra.pinv(CopiedFindMonomialBasis.find_monomial_basis(f, 2, true)[2])
    gVec = AutomatedScript.convert_p_to_m(g,AutomatedScript.gen_exp_vec(n,n*d-n))
    gJS = partsMat*gVec
    gc = 
    gcpartials = [ derivative(gc[1], i) for i in 1:(n+1) ]
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
            println(V)
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
function computeReductionChainLA(I,n,d,m,S,parts,R,PR,ResRing)
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
            RPoly = computeRPolyLA(V,S,n,d,parts,R,PR)
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
function computeReductionOfPolyLA(poly,n,d,S,parts,R,PR,ResRing)
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
        RChain = computeReductionChainLA(I,n,d,i[2],S,parts,R,PR,ResRing)
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
function computeReductionOfTransformLA(FT,n,d,p,N,S,parts,R,PR)
    result = 0
    for i in FT
        for j in i
            s = N + j[2] - 1
            ResRing = residue_ring(ZZ,Int1024(p)^s)
            result = computeReductionOfPolyLA(j,n,d,S,parts,R,PR,ResRing)[1]
        end
    end
    return [result,n]
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
    URing, UVars = PolynomialRing(R, ["u$i" for i in 0:n])
    PURing, Vars = PolynomialRing(URing, ["x$i" for i in 0:n])
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
function computeRPolyLA(V,S,n,d,parts,R,PR)
    URing, UVars = PolynomialRing(R, ["u$i" for i in 0:n])
    PURing, Vars = PolynomialRing(URing, ["x$i" for i in 0:n])
    #this is temporary
    polynomial = Vars[1]^5 + Vars[2]^5 + Vars[3]^5 + Vars[1]*(Vars[2]^3)*Vars[3]
    parts = [ derivative(polynomial, i) for i in 1:(n+1) ]
    #
    ev = AutomatedScript.gen_exp_vec(n+1,d*n-n)
    monomials = AutomatedScript.gen_mon(ev,URing,PURing)
    reductions = []
    for m in monomials
        push!(reductions, computeReductionLA(UVars,V,S,n,d,[m,1],parts,[],URing,PURing,Vars)[1])
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
            #println(typeof(coeff(fj,alpha)^p))
            sum = sum + p^(m-1)*(D[j+1]*(coeff(map_coefficients(lift,fj),alpha)^p) % p^s)*monomial
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
        push!(result, applyFrobenius(n,d,f,N,p,b,R,PR))
    end
    return result
end
    

end
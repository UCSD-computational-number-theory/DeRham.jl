module ControlledReduction

using Oscar
include("AutomatedScript.jl")

function computeReduction(U,V,S,n,d,g,parts,ev,R,PR,Vars)
    SC = []
    gensJS = []
    for i in 0:n
        if i in S
            push!(gensJS,parts[i])
        else
            push!(gensJS,Vars[i+1]*parts[i+1])
            push!(SC,i)
        end
    end
    # get gi's using psuedoinverse
    gc, t = reduce_with_quotients(g[1],parts)
    gcpartials = [ derivative(gc[1], i) for i in 1:(n+1) ]
    XS =  prod(PR(Vars[i+1]) for i in S; init = PR(1))
    return [sum((PR(U[i+1])*div(XS,Vars[i+1])*gc[i+1] + XS*gcpartials[i+1]) for i in S; init = PR(0)) + sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*Vars[i+1]*gcpartials[i+1]) for i in SC; init = PR(0)), g[2]-1]

end

function computeR(u,v,s,n,d,R,PR,vars)
    ev = AutomatedScript.gen_exp_vec(n,d*n-n)
    monomials = AutomatedScript.gen_mon(ev,R,PR)
    reductions = []
    for m in monomials
        push!(reductions,computeReduction(u,v,s,n,d,m,ev,R,PR,vars))
    end
    return transpose(AutomatedScript.convert_p_to_m(reductions,ev))
end

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
            sum = sum + R(p^(m-1))*(R(D[j+1]*(coeff(fj,alpha)^p) % p^s))*monomial
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
        push!(temp, applyFrobeniusToMon(n,d,map_coefficients(lift,f),N,p,ev,poly[2],R,PR))
    end
    return temp
end



end
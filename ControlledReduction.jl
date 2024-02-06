module ControlledReduction

using Oscar
include("AutomatedScript.jl")

function computeReduction(U,V,S,n,d,g,ev,R,PR,vars)
    SC = []
    for i in 0:n
        if i in S
        else
            push!(SC,i)
        end
    end
    # get gi's using psuedoinverse
    gc = []
    gcpartials = [ derivative(gc[i+1], i) for i in 1:(n+1) ]
    XS =  prod(PR(vars[i+1]) for i in S)
    return [sum((PR(U[i+1])*div(XS,vars[i+1])*gc[i+1] + XS*gcpartials[i+1]) for i in S) + sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*var[i+1]*gcpartials[i+1]) for i in SC), g[2]-1]

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
    D = zeros(N)
    for j in 0:(N-1)
        D[j+1] = sum((-1)^(i+j)*binomial(-m,i)*binomial(i,j) for i in j:(N-1))
    end
    return D
end

function applyFrobenius(n,d,f,N,p,beta,m,R,PR)
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
            sum = sum + R(p^(m-1))*(R(D[j+1]*(coeff(fj,alpha)^p) % p^s))*monomial
        end
        push!(result, [sum, p*(m+j)])
    end
    return result
end

end
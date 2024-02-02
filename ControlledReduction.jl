module ControlledReduction

using Oscar

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
    return sum((PR(U[i+1])*div(XS,vars[i+1])*gc[i+1] + XS*gcpartials[i+1]) for i in S) + sum((PR(U[i+1]+1)*XS*gc[i+1] + XS*var[i+1]*gcpartials[i+1]) for i in SC)

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

end
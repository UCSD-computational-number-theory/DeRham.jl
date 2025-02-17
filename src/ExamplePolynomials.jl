# TODO: PR this into Oscar, hopefully

function artin_shreier_curve(num,denom)
  p = characteristic(parent(num))
  R, (x,y) = polynomial_ring(GF(p),["x","y"]) 

  lhs = evaluate(num,x)
  prod_factor = evaluate(denom,x)

  #TODO: make this into a scheme in Oscar's language
  prod_factor * (y^p - y) - lhs
end

function projective_artin_shreier_curve(num,denom)
  nonhomog = artin_shreier_curve(num,denom)
  P = parent(nonhomog)
  H = homogenizer(P, "z")
  H(nonhomog)
end

#function artih_shreier_curve(p,f)
#  TODO: 
#
#end

function artin_shreier_witt_curve(p,num,denom)
  #TODO: make heuristic algorithm for delta_1 in MMPSingularities.jl

  error("Artin-Shreier-Witt curves not implemented yet")
end

"""
Gives the fermat hypersurface of degree d and dimension n-2 
in characterstic p.
"""
function fermat_hypersurface(n,d,p)
    F = GF(p)
    R, _ = polynomial_ring(F,n)
    

    C = MPolyBuildCtx(R)
    exp_vec = zeros(Int,n)
    for i in 1:n
        zero!(exp_vec)
        exp_vec[i] = d
        push_term!(C,F(1),exp_vec)
    end
    
    finish(C)
end

"""
Gives the hypersurface of degree d and dimension n-2
which is defined by the polynomial which has every 
term with coefficient 1.
"""
function full_hypersurface(n,d,p)
    F = GF(p)
    R, _ = polynomial_ring(F,n)

    evs = gen_exp_vec(n,d)

    C = MPolyBuildCtx(R)

    for ev in evs
        push_term!(C,F(1),ev)
    end
    
    finish(C)
end

"""
Gives a hypersurface of degree d and dimension n-2
which has uniform random coefficients in characteristic p
"""
function random_hypersurface(n,d,p)
    F = GF(p)
    R, _ = polynomial_ring(F,n)

    evs = gen_exp_vec(n,d)
    cs = rand(0:p-1,length(evs))

    C = MPolyBuildCtx(R)

    for i in eachindex(evs)
        push_term!(C,F(cs[i]),evs[i])
    end
    
    finish(C)
end

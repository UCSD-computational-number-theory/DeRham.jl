module StandardReduction

using Oscar

function stdreduction(polynomial, num_vars, f_exp)
       partials = []
       for i in axes(1,num_vars)
           push!(partials, derivatives(polynomial, x[i])
       end
       G = groebner_basis(polynomial)
       result = 0
       for i in axes(1, num_vars)
	   result = result + derivative(G[I], i)
       end
       
       return (f_exp - 1)
end
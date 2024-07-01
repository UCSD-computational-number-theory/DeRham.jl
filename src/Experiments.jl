module Experiments

using Oscar


#Very rough code to generate monomials, will clean up later


#Generates an array of vectors corresponding to the degrees of the variables in the monomial
function gen_exp_vec(n,d)
    result = []
    if n == 1
        return [[d]]
    end
    if d == 1
        for i in 1:n
            s = zeros(Int64,n)
            s[i] = 1
            push!(result,s)
        end
        return result
    end
    for i in 0:d
        y = gen_exp_vec(n-1,d-i)
        for j in axes(y,1)
            append!(y[j],i)
        end
        append!(result,y)
    end
    return result
end

#generates the monomials
function gen_mon(exp_vec, R)
    result = []
    for i in axes(exp_vec,1)
        B = MPolyBuildCtx(R)
        push_term!(B, R(1), exp_vec[i])
        monomial = finish(B)
        push!(result,monomial)
    end
    return result
end

function compute_relations(monomials, partials)
    result = []
    for i in axes(monomials,1)
        for j in axes(partials,1)
            push!(result, monomials[i]*partials[j])
        end
    end
    return result
end

function convert_p_to_m(polys, expvec)
    result = []
    for i in axes(polys,1)
        temp = []
        for j in axes(expvec,1)
            append!(temp, [coeff(polys[i], expvec[j])])
        end
        push!(result, transpose(temp))
    end
    return vcat(result...)
end

function convert_m_to_p(mat, expvec, R, PR)
    result = []
    for i in axes(mat,1)
        B = MPolyBuildCtx(PR)
        for j in axes(expvec,1)
            push_term!(B, mat[i,j], expvec[j])
        end
        push!(result,finish(B))
    end
    return result
end

end
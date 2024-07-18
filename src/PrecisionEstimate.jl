# module PrecisionEstimate

#p = 41
#r = 7 # p-adic precision
#n = 2

"""
    frobenius_precision()

Determines the necessary p-adic precision r for the Frobenius matrix  
"""
function frobenius_precision(k, q)
    if k == 2 && q == 7
        r = 2
    elseif k == 21 && q == 7
        r = 10
    end 
    return r
end 

"""
    relative_precision()

Determines the relative p-adic precision [r_m] for the basis of cohomology 
"""
function relative_precision(k, q)
    if k == 2 && q == 7
        r_m = [2,3]
    elseif k == 21 && q == 7
        r_m = [2,3,3]
    end
    return r_m
end 

function series_precision(r_m, p, n)
    if r_m == [2,3] && p == 7 && n == 2
        N_m = [2,2]
    elseif r_m == [2,3,3] && p == 7 && n == 3
        N_m = [4,4,3]
    end 
    return N_m
end

"""
    algorithm_precision()

Determines M such that the algorithm works in Z/p^M Z 
"""
function algorithm_precision(r_m, N_m, p)
    if r_m == [2,3] && N_m == [2,2] && p == 7
        M = 3
    elseif r_m == [2,3,3] && N_m == [4,4,3] && p == 7
        M = 6
    end
    return M
end 




#=
function compute_N(p, r, m, n)
    num = n * log(n + r)
    denom = (n + r) * log(p) - n
    return ceil(Int, -m + (n + r) * (1 + num / denom))
end

function compute_precisions_each(p, r, n)
    result = []
    for m = 1:n
        push!(result,compute_N(p, r, m, n))
    end
    return result
end

#compute_precisions_each(p, r, n)
=#
#end

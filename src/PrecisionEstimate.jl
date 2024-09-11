# module PrecisionEstimate

#p = 41
#r = 7 # p-adic precision
#n = 2

#TODO: this file needs to be dual licensed under the GPL because I copied hardcoded values from controlledreduction


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
        r_m = [1,2]
    elseif k == 21 && q == 7
        r_m = [2,3,3]
    elseif k == 12 && q == 7
        r_m = [3,4]
    end
    return r_m
end 

#function series_precision(r_m, p, n)
function series_precision(p, n, d)
    #if r_m == [2,3] && p == 7 && n == 2
    #    N_m = [2,2]
    #elseif r_m == [2,3,3] && p == 7 && n == 3
    #    N_m = [4,4,3]
    #end 
    #return N_m
    if ( n == 2 && d == 3 ) 
        if ( p < 5) 
            N = [2 3]
        elseif ( 5 <= p && p < 17 ) 
            N = [2,2]
        elseif ( 17 <= p  ) 
            N = [0,1]
        end
    elseif ( n == 2 && d == 4 ) 
        if ( p < 5) 
            N = [4,4]
        elseif ( 5 <= p && p < 17 ) 
            N = [3,3]
        elseif ( 17 <= p  ) 
            N = [2,2]
        end
    elseif ( n == 2 && d == 5 ) 
        if ( p < 5) 
            N = [6,6]
        elseif ( 5 <= p && p < 7 ) 
            N = [4,5]
        elseif ( 7 <= p  ) 
            N = [4,4]
        end
    elseif ( n == 2 && d == 6 ) 
        if ( p < 5) 
            N = [8,9]
        elseif ( 5 <= p && p < 7 ) 
            N = [7,7]
        elseif ( 7 <= p && p < 11 ) 
            N = [6,7]
        elseif ( 11 <= p  ) 
            N = [6,6]
        end
    elseif ( n == 2 && d == 7 ) 
        if ( p < 5) 
            N = [11,11]
        elseif ( 5 <= p && p < 7 ) 
            N = [10,10]
        elseif ( 7 <= p && p < 11 ) 
            N = [10,10]
        elseif ( 11 <= p && p < 17 ) 
            N = [9,9]
        elseif ( 17 <= p  ) 
            N = [8,8]
        end
    elseif ( n == 2 && d == 8 ) 
        if ( p < 5) 
            N = [14,14]
        elseif ( 5 <= p && p < 7 ) 
            N = [13,13]
        elseif ( 7 <= p && p < 13 ) 
            N = [13,13]
        elseif ( 13 <= p && p < 17 ) 
            N = [12,13]
        elseif ( 17 <= p  ) 
            N = [11,11]
        end
    elseif ( n == 2 && d == 9 ) 
        if ( p < 5) 
            N = [18,18]
        elseif ( 5 <= p && p < 7 ) 
            N = [16,16]
        elseif ( 7 <= p && p < 11 ) 
            N = [16,16]
        elseif ( 11 <= p && p < 17 ) 
            N = [16,16]
        elseif ( 17 <= p  ) 
            N = [15,15]
        end
    elseif ( n == 3 && d == 4 ) 
        if ( p < 5) 
            N = [7,7,8]
        elseif ( 5 <= p && p < 7 ) 
            N = [4,5,5]
        elseif ( 7 <= p && p < 23 ) 
            N = [4,4,3]
        elseif ( 23 <= p && p < 43 ) 
            N = [3,3,3]
        elseif ( 43 <= p  ) 
            N = [3,3,2]
        end
    elseif ( n == 3 && d == 5 ) 
        if ( p < 5) 
            N = [11,11,10]
        elseif ( 5 <= p && p < 7 ) 
            N = [8,8,9]
        elseif ( 7 <= p && p < 11 ) 
            N = [8,8,7]
        elseif ( 11 <= p && p < 23 ) 
            N = [7,7,6]
        elseif ( 23 <= p && p < 29 ) 
            N = [6,6,6]
        elseif ( 29 <= p  ) 
            N = [6,6,5]
        end
    elseif ( n == 3 && d == 6 ) 
        if ( p < 5) 
            N = [17,18,17]
        elseif ( 5 <= p && p < 7 ) 
            N = [15,15,14]
        elseif ( 7 <= p && p < 11 ) 
            N = [15,15,14]
        elseif ( 11 <= p && p < 17 ) 
            N = [14,14,13]
        elseif ( 17 <= p && p < 23 ) 
            N = [13,13,12]
        elseif ( 23 <= p  ) 
            N = [12,12,11]
        end
    end
    N
end

"""
    algorithm_precision()

Determines M such that the algorithm works in Z/p^M Z 
"""
#function algorithm_precision(r_m, N_m, p)
function algorithm_precision(p,n,d,r_m,N_m)
    #if r_m == [2,3] && N_m == [2,2] && p == 7
    #    M = 3
    #elseif r_m == [2,3,3] && N_m == [4,4,3] && p == 7
    #    M = 6
    #end
    #return M
    s_m = [i+x-1 for (i,x) in enumerate(N_m)]
    s_m_valuation = [valuation(ZZ(factorial(big(p*s-1))), p) for s in s_m]

    precision = maximum([r_m[m] + s_m_valuation[m] - m + 1 for m in 1:length(r_m)])
    println(precision)

    if ( n == 2 && d == 3 ) 
        if ( p < 5) 
            @assert precision == 5;
        elseif ( 5 <= p && p < 17 ) 
            @assert precision == 3;
        elseif ( 17 <= p  ) 
            @assert precision == 1;
        end
    elseif ( n == 2 && d == 4 ) 
        if ( p < 5) 
            @assert precision == 7;
        elseif ( 5 <= p && p < 17 ) 
            @assert precision == 5;
        elseif ( 17 <= p  ) 
            @assert precision == 3;
        end
    elseif ( n == 2 && d == 5 ) 
        if ( p < 5) 
            @assert precision == 12;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 9;
        elseif ( 7 <= p  ) 
            @assert precision == 7;
        end
    elseif ( n == 2 && d == 6 ) 
        if ( p < 5) 
            @assert precision == 19;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 13;
        elseif ( 7 <= p && p < 11 ) 
            @assert precision == 13;
        elseif ( 11 <= p  ) 
            @assert precision == 11;
        end
    elseif ( n == 2 && d == 7 ) 
        if ( p < 5) 
            @assert precision == 23;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 20;
        elseif ( 7 <= p && p < 11 ) 
            @assert precision == 19;
        elseif ( 11 <= p && p < 17 ) 
            @assert precision == 17;
        elseif ( 17 <= p  ) 
            @assert precision == 15;
        end
    elseif ( n == 2 && d == 8 ) 
        if ( p < 5) 
            @assert precision == 30;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 26;
        elseif ( 7 <= p && p < 13 ) 
            @assert precision == 25;
        elseif ( 13 <= p && p < 17 ) 
            @assert precision == 25;
        elseif ( 17 <= p  ) 
            @assert precision == 21;
        end
    elseif ( n == 2 && d == 9 ) 
        if ( p < 5) 
            @assert precision == 41;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 33;
        elseif ( 7 <= p && p < 11 ) 
            @assert precision == 32;
        elseif ( 11 <= p && p < 17 ) 
            @assert precision == 31;
        elseif ( 17 <= p  ) 
            @assert precision == 29;
        end
    elseif ( n == 3 && d == 4 ) 
        if ( p < 5) 
            @assert precision == 16;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 9;
        elseif ( 7 <= p && p < 23 ) 
            @assert precision == 6;
        elseif ( 23 <= p && p < 43 ) 
            @assert precision == 5;
        elseif ( 43 <= p  ) 
            @assert precision == 4;
        end
    elseif ( n == 3 && d == 5 ) 
        if ( p < 5) 
            @assert precision == 21;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 17;
        elseif ( 7 <= p && p < 11 ) 
            @assert precision == 14;
        elseif ( 11 <= p && p < 23 ) 
            @assert precision == 12;
        elseif ( 23 <= p && p < 29 ) 
            @assert precision == 11;
        elseif ( 29 <= p  ) 
            @assert precision == 10;
        end
    elseif ( n == 3 && d == 6 ) 
        if ( p < 5) 
            @assert precision == 38;
        elseif ( 5 <= p && p < 7 ) 
            @assert precision == 29;
        elseif ( 7 <= p && p < 11 ) 
            @assert precision == 28;
        elseif ( 11 <= p && p < 17 ) 
            @assert precision == 26;
        elseif ( 17 <= p && p < 23 ) 
            @assert precision == 24;
        elseif ( 23 <= p  ) 
            @assert precision == 22;
        end
    end
    return precision
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

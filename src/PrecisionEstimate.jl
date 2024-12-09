# module PrecisionEstimate

#p = 41
#r = 7 # p-adic precision
#n = 2

#TODO: parts of this file needs to be dual licensed under the GPL because I copied hardcoded values from controlledreduction


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
    impreciselog(p,a)

log with any base, via the log rules,
though I think there should be floating point error.
"""
impreciselog(p,a) = log(a) / log(p)

"""
    ibm_derham_bound(p,m,n)

ibm stands for integral basis multiplier

Calculates the smallest integer t such that
p^t* ω is in the integral latticee of the cohomology,
for all ω in the griffiths-dwork basis of cohomology.
Here, `m` is the pole order, 
and the bound we have depends only on p, m, and n.

This is called ϕ in Costa's thesis, and it is 
called f in Abbott-Kedlaya-Roe.
This method implements theorem 3.4.6 only of 
Abbott-Kedlaya-Roe.
"""
function ibm_derham_bound(p,m,n)

    sum = 0
    for i in 1:n
        sum += floor(Int,impreciselog(p,max(1,m-i)))
    end

    sum
    #floor(Int,(n)*impreciselog(p,m))
end

"""
    ibm_constant(p,m,ibm_bound)

Calcuates the `N` in theorem 3.4.9 of Abbott-Kedlaya-Roe
"""
function ibm_constant(p,m,ibm_bound)
    g(m,i) = valuation(binomial(-m,i),p)

    NN(l) = ibm_bound[(m+l)*p] - l - g(m,l)

    l = 1
    N = NN(1)
    while 0 < NN(l)
        N = max(N,NN(l))
        l += 1
    end
    #verbose && println("In AKR notation, f((m+l)p) - l - g(m,l) < 0 at l = $l")

    N
end

"""
"""
function ibm_crank_bound(p,n,M)
    nChunks = ceil(Int,M/p)
    nBounds = nChunks * p
    bounds = ibm_derham_bound.(p,1:(nBounds*p^2),n)

    for m = 1:nChunks
#        println(bounds)

        # find the best bound for f(p*m)
        c = 1
        while true
            N = ibm_constant(p,m,bounds)
            #println("N: $N")
            if N ≤ n - 1 + bounds[m]
                if n - 1 + bounds[m] ≤ bounds[p*m] 
                    bounds[p*m] = n - 1 + bounds[m]
                end
                #bounds[p*m] = n - 1 + bounds[m]
                #print("reached log bound n - 1 + f(m) = $(bounds[p*m])")
                #println(" at for chunk $m, i.e. f($(p*m))")
                break
            elseif bounds[p*m] ≤ N
                #println("crank stabilized at N=$N for chunk $m, i.e. f($(p*m))")
                if bounds[p*m] < N
                    #println("original bound is better than N!")
                end
                break
            else
                # go another round!
                bounds[p*m] = N
                #println("go another round!")
            end
            #println("finished crank number $c")
        end

        # update the elements before
        for k = (p*m)-p+1:p*m
            if bounds[m*p] < bounds[k]
                bounds[k] = bounds[m*p]
            end
        end

    end

    #all(bounds .<= ibm_derham_bound.(p,1:nBounds*p^2,n)) || error("assertion failure")

    bounds[1:nBounds]
end

#"""
#This algorithm is not proven to be correct,
#please don't use it.
#"""
#function ibm_crank_bound_longwise(p,n,M;MM=nothing)
#    if MM == nothing
#        MM = M * p^2
#    end
#
#    # start with de rham bound
#    bounds = ibm_derham_bound.(p,1:MM,n)
#
#    numCranks = 3
#    
#    for i = 1:numCranks
#
#        m = 1
#        while true
#            N = ibm_constant(p,m,bounds)
#            if N < n - 1 + bounds[m]
#                #bounds[p*m] = n - 1 + bounds[m]
#                for mm = p*m:p:MM
#                    bounds[mm] = n - 1 + bounds[divexact(mm,p)]
#                end
#                println("log bound takes over at m = $(m*p)")
#                break
#            else
#                bounds[p*m] = N
#            end
#            m = m + 1
#        end
#
#        # adjust bounds down
#        for l = p:p:MM
#            for k = l-p+1:l
#                if bounds[p] < bounds[k]
#                    bounds[k] = bounds[p]
#                end
#            end
#        end
#
#    end
#
#    bounds
#end


"""
    calculate_relative_precision(HP, slope, hodge_numbers, weight, p)

Calculates the vector of relative precisions r_m 

INPUTS: 
* "polygon" -- A SlopesPolygon struct describing a polygon whose values describe the
divisibility of the roots. Note: it is traditional to use the hodge polygon here,
but if you know the Newton Polygon (e.g. you have a K3 surface and you already have the
Artin-Mazur height) then the Newton Polygon will do just as well.
* "weight" -- integer, the motivic weight of the cohomology group, equals to the dimension of the hypersurface
* "p" -- integer, prime number 
"""
function calculate_relative_precision(polygon, weight, p) 

#   * "HP" -- list, the i-th item corresponds to the height of above i in the Hodge polygon
#   * "slope" -- list, the i-th item corresponds to the slope of the i-th segment in the Hodge polygon
#   * "hodge_numbers" -- list, the list of Hodge numbers 
    HP = polygon.values
    slope = Int.(polygon.slopesbefore)
    hodge_numbers = polygon.slopelengths

    r_vector = [0 for i in 1:length(hodge_numbers)]
    Pdeg = length(HP) - 1  # degree of the L-polynomial P_n(X,T)
    max_digits = 0
    
    for i in 1:(div(Pdeg, 2) + 1)
        r = ceil(Int, log(2*Pdeg/i)/log(p) + i*weight*0.5) - HP[i+1]

        if r >= max_digits
            max_digits = r 
            for j in 0:slope[i+1]
                r_vector[j+1] = r 
            end 

            for j in slope[i+1]+1:length(hodge_numbers)-1
                r = r-1
                r_vector[j+1] = r
            end 
        end 
    end 
    println(max_digits)

    return reverse(r_vector)
end 

function calculate_series_precision(p,n,r_m)

    #println("Relative precision: $r_m")
    #if !(p < 2*(n) + maximum(r_m))
    #    #println("Unsupported p: $p. Returning nonsense.")
    #    #return zeros(Int,n)
    #    N = [(r_m[m] == 0 ? 0 : n+1 + r_m[m] - (m+1)) for m in 1:n]
    #    return N
    #end

    # TODO: prove that this is correct
    #
    # We may assume in (1.8) of Costa's thesis that the worst case is i=0

    bounds = ibm_crank_bound(p,n,p*n)
    #bounds = ibm_derham_bound.(p,1:p*n,n)

    #TODO: why does ibm_derham_bound sometimes give answers that are better than the crank?
    
    N = [(r_m[m] == 0 ? 0 : r_m[m] - m + 1 + bounds[p*m]) for m in 1:n]

    println("Series precision: $N")
    N
end


#"""
#    relative_precision()
#
#Determines the relative p-adic precision [r_m] for the basis of cohomology 
#"""
#function relative_precision(k, q)
#    if k == 2 && q == 7
#        r_m = [1,2]
#    elseif k == 21 && q == 7
#        r_m = [2,3,3]
#    elseif k == 12 && q == 7
#        r_m = [3,4]
#    end
#    return r_m
#end 

#function series_precision(r_m, p, n)
function series_precision(p, n, d, r_m)
    #if r_m == [2,3] && p == 7 && n == 2
    #    N_m = [2,2]
    #elseif r_m == [2,3,3] && p == 7 && n == 3
    #    N_m = [4,4,3]
    #end 
    #return N_m
    N = [0,0]

    # begin copypasta
    if ( n == 2 && d == 3 ) 
        if ( p < 5) 
            N = [2,3]
        elseif ( 5 <= p && p < 17 ) 
            N = [2,2]
        elseif ( 17 <= p  ) 
            N = [1,1]#[0,1]
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
    # end copypasta

    if N == [0,0]
        N = calculate_series_precision(p,n,r_m)
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
    s_m_valuation = [valuation(ZZ(factorial(big(p*s-1))), ZZ(p)) for s in s_m]

    maximum([r_m[m] + s_m_valuation[m] - m + 1 for m in 1:length(r_m)])
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

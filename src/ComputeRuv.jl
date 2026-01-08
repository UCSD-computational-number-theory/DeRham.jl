
"""
    computeRuvS
Returns a list of n+2 matrices for the Ruv map as in Proposition 1.15
"""
function computeRuvS(V,S,f,pseudoInverseMat,cache,params)
    vars_reversed = params.vars_reversed
    termorder = params.termorder
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    #R = coefficient_ring(parent(f))
    R = base_ring(pseudoInverseMat)

    g_length = binomial(n*d,n*d-n)


    ev1 = (cache[n*d - n])#gen_exp_vec(n+1,n*d-n,termorder)
    ev2 = (cache[n*d-n+d-length(S)])#gen_exp_vec(n+1,n*d-n+d-length(S),termorder)
    ev3 = (cache[n*d-n-length(S)+1])#gen_exp_vec(n+1,n*d-n-length(S)+1,termorder)
    ev4 = (cache[n*d-n-length(S)])#gen_exp_vec(n+1,n*d-n-length(S),termorder)
    explookup = cache[n*d - n,:reverse]
    explookup2 = cache[n*d-n+d-length(S),:reverse]
    
    temp = Vector{Int64}(undef, n+1)
    #MS2 = matrix_space(R, length(ev2),1)
    matrixtype = typeof(pseudoInverseMat) #eltype(valtype(Ruvs))

    result = Vector{matrixtype}(undef, n+2)
    for i in 1:n+2
        result[i] = zero_matrix(R,g_length,g_length)
    end

    gJS = zeros(R,size(pseudoInverseMat,1))

    Stilda = zeros(Int, n+1)
    for i in S
        Stilda[i+1] = 1
    end
    distances = Vector{Int64}(undef, n+1)
    for i in 1:(n+1)
        if Stilda[i] == 1
            distances[i] = length(ev3)
        else
            distances[i] = length(ev4)
        end
    end
    distance = 0
    for i in 1:length(ev1)  # indexing over the columns of the Ruv matrices 
        mon = Vector{Int64}(undef, n+1)
        for m in 1:(n+1)  # forming the polynomial x^v*g / x^S
            mon[m] = ev1[i][m] + V[m] - Stilda[m]
        end

        #gVec = zero_matrix(R,length(ev2),1)#MS2()
        #for j in 1:length(ev2)
        #    if ev2[j] == mon
        #        gVec[j] = one(R)#R(1)
        #    else
        #        gVec[j] = zero(R)#R(0)
        #    end
        #end
        #gJS = pseudoInverseMat*gVec   # writing x^v*g/x^S as \sum_{i\in S} g_i*\partial_i f + \sum_{i\notin S} g_i*x_i*\partial_if
        #println("After LingAlg problem: $gJS")

        # gVec above has a one at coordinate j, so we may do the matmul
        # by taking the column.
        #ind = findfirst(==(mon),ev2)
        ind = get(explookup2,mon,-1)
        if ind == -1
            #gJS = zeros(base_ring(pseudoInverseMat),size(pseudoInverseMat,1))
            zero!(gJS)
        else
            #gJS = pseudoInverseMat[:,ind]
            for h in 1:size(pseudoInverseMat,1)
                gJS[h] = pseudoInverseMat[h,ind]
            end
        end
        #println(gJS)

        distance = 0
        for j in 1:(n+1)
            if Stilda[j] == 1
                for k in 1:length(ev3)
                    for m in 1:(n+1)
                        if m == j
                            temp[m] = ev3[k][m] + Stilda[m] - 1
                        else
                            temp[m] = ev3[k][m] + Stilda[m]
                        end
                    end
                    l = get(explookup,temp,-1)
                    
                    result[j+1][l,i] = result[j+1][l,i] + gJS[distance+k]
                    result[1][l,i] = result[1][l,i] + (ev3[k][j])*gJS[distance+k]
                end
            else
                for k in 1:length(ev4)
                    for m in 1:(n+1)
                        temp[m] = ev4[k][m] + Stilda[m]
                    end
                    l = get(explookup,temp,-1)
                    result[j+1][l,i] = result[j+1][l,i] + gJS[distance+k]
                    result[1][l,i] = result[1][l,i] + (ev4[k][j])*gJS[distance+k] + gJS[distance+k]
                end
            end
            distance = distance + distances[j]
        
        end
    end
    #TODO: there is some sort of race condition on 
    # our dictionary, and putting this print statement here 
    # fixes it. Later, we'll need to fix this for reals
    (2 < params.verbose && 1 < Threads.nthreads()) && println("New key: $V")
    (2 < params.verbose) && begin
        nnzs = count.(!=(0),result)
        total_entries = g_length^2
        percentages = nnzs * 100 ./ total_entries 
        println("Ruv calculated for V = $V with densities")
        for i in 1:length(result)
            println("    Matrix $i: $(nnzs[i]) of $total_entries entreis, $(percentages[i]) %")
        end
    end
    
    return result
end


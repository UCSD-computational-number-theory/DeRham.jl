
#function precomputeRuvs(S,f,pseudoInverseMat,Ruvs,cache,params)
#    d = total_degree(f)
#    n = nvars(parent(f)) - 1
#
#    evs = gen_exp_vec(n+1,d,params.termorder)
#    Threads.@threads for V in evs
#    #for V in evs
#        computeRuvS(V,S,f,pseudoInverseMat,Ruvs,cache,params)
#    end
#
#    
#    #(6 < params.verbose) && begin 
#    #    println("V = [1,2,0] : $(Ruvs[[1,2,0]])")
#    #    println("V = [0,2,1] : $(Ruvs[[0,2,1]])")
#    #end
#end

#function computeRuvOld(V,S,f,pseudoInverseMat,Ruvs,cache,params)
#    vars_reversed = params.vars_reversed
#    termorder = params.termorder
#    n = nvars(parent(f)) - 1
#    d = total_degree(f)
#    R = coefficient_ring(parent(f))
#    MS1 = matrix_space(R, binomial(n*d,n*d-n), binomial(n*d,n*d-n))
#    if S != collect(0:n-1)
#        throw(ArgumentError("computeRuvOld only works for S={0...n-1}"))
#    end
#    if haskey(Ruvs, V)
#        return get(Ruvs, V, 0)
#    else
#        (4 < params.verbose) && println("New key: $V")
#    end
#    ev1 = cache[n*d-n]#gen_exp_vec(n+1,n*d-n,termorder)
#    ev2 = cache[n*d-n+d-length(S)]#gen_exp_vec(n+1,n*d-n+d-length(S),termorder)
#    ev3 = cache[n*d-n-length(S)+1]#gen_exp_vec(n+1,n*d-n-length(S)+1,termorder)
#    explookup = cache[n*d - n,:reverse]
#    temp = Vector{Int64}(undef, n+1)
#    MS2 = matrix_space(R, length(ev2),1)
#    result = Vector{typeof(MS1())}(undef, n+2)
#    for i in 1:n+2
#        result[i] = MS1(0)
#    end
#    Stilda = zeros(Int, length(S))
#    for i in S
#        Stilda[n+1-i] = 1
#    end
#    for i in 1:length(ev1)
#        mon = Vector{Int64}(undef, n+1)
#        for m in 1:(n+1)
#            mon[m] = ev1[i][m] + V[m] - Stilda[m]
#        end
#        gVec = MS2()
#        for j in 1:length(ev2)
#            if ev2[j] == mon
#                gVec[j] = R(1)
#            else
#                gVec[j] = R(0)
#            end
#        end
#        gJS = pseudoInverseMat*gVec
#        #println("After LingAlg problem: $gJS")
#        for j in 1:(n+1)
#            for k in 1:length(ev3)
#                for m in 1:(n+1)
#                    if m == n+1-j+1
#                        temp[m] = ev3[k][m] + Stilda[m] - 1
#                    else
#                        temp[m] = ev3[k][m] + Stilda[m]
#                    end
#                end
#                #print("ev1[l]: $((ev1[l],typeof(ev1[l])));")
#                #print("ev3[k]: $((ev3[k],typeof(ev3[k])));") 
#                #println(" $(ev1[l] == ev3[k])")
#                l = get(explookup,temp,-1)
#
#                if params.vars_reversed == true
#                    result[j+1][l,i] = gJS[Int((j-1)*(length(gJS)/(n+1))+1)+k-1,1]
#                    result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((j-1)*(length(gJS)/(n+1))+1)+k-1,1]
#                else
#                    result[j+1][l,i] = gJS[Int((n+1-j)*(length(gJS)/(n+1))+1)+k-1,1]
#                    result[1][l,i] = result[1][l,i] + (ev3[k][n+1-j+1])*gJS[Int((n+1-j)*(length(gJS)/(n+1))+1)+k-1,1]
#                end
#            end
#        end
#    end
#    get!(Ruvs, V, result)
#    return result
#end

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
    # if the following isn't matched in Ruvs, this method will error.
    R = base_ring(pseudoInverseMat)

    g_length = binomial(n*d,n*d-n)

    #if haskey(Ruvs, V)
    #    return get(Ruvs, V, 0)
    #else
    #    (4 < params.verbose) && println("New key: $V")
    #end

    # instead of doing this, lets reverse and reverse back:
     
    ev1 = (cache[n*d - n])#gen_exp_vec(n+1,n*d-n,termorder)
    ev2 = (cache[n*d-n+d-length(S)])#gen_exp_vec(n+1,n*d-n+d-length(S),termorder)
    ev3 = (cache[n*d-n-length(S)+1])#gen_exp_vec(n+1,n*d-n-length(S)+1,termorder)
    ev4 = (cache[n*d-n-length(S)])#gen_exp_vec(n+1,n*d-n-length(S),termorder)
    explookup = cache[n*d - n,:reverse]
    explookup2 = cache[n*d-n+d-length(S),:reverse]
    # if vars_reversed && !cache.vars_reversed
    #     for vec in ev1
    #         reverse!(vec)
    #     end 

    #     for vec in ev2
    #         reverse!(vec)
    #     end 

    #     for vec in ev3
    #         reverse!(vec)
    #     end 

    #     for vec in ev4
    #         reverse!(vec)
    #     end     
    # end 
    
    #ev1 = cache[n*d - n]#gen_exp_vec(n+1,n*d-n,termorder)
    #ev2 = cache[n*d-n+d-length(S)]#gen_exp_vec(n+1,n*d-n+d-length(S),termorder)
    #ev3 = cache[n*d-n-length(S)+1]#gen_exp_vec(n+1,n*d-n-length(S)+1,termorder)
    #ev4 = cache[n*d-n-length(S)]#gen_exp_vec(n+1,n*d-n-length(S),termorder)
    #explookup = cache[n*d - n,:reverse]

    #if vars_reversed
    #    reverse!.(ev1)
    #    reverse!.(ev2)
    #    reverse!.(ev3)
    #    reverse!.(ev4)
    #end
    
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
        # if vars_reversed
        #     Stilda[n+1-i] = 1
        # else
            Stilda[i+1] = 1
        # end
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
                    #print("ev1[l]: $((ev1[l],typeof(ev1[l])));")
                    #print("ev3[k]: $((ev3[k],typeof(ev3[k])));") 
                    #println(" $(ev1[l] == ev3[k])")
                    # if vars_reversed && !cache.vars_reversed
                    #     reverse!(temp)
                    #     l = get(explookup, temp, -1)
                    #     reverse!(temp)
                    # else 
                        l = get(explookup,temp,-1)
                    # end 
                    #println(j, " ", l, " ", i, " ", k, " ", distance)
                    
                    result[j+1][l,i] = result[j+1][l,i] + gJS[distance+k]
                    #if V == [0,3,0] && j == 3 && l == 9 && i == 5 && k == 9
                    #    error()
                    #end
                    result[1][l,i] = result[1][l,i] + (ev3[k][j])*gJS[distance+k]
                end
            else
                for k in 1:length(ev4)
                    for m in 1:(n+1)
                        temp[m] = ev4[k][m] + Stilda[m]
                    end
                    #print("ev1[l]: $((ev1[l],typeof(ev1[l])));")
                    #print("ev3[k]: $((ev3[k],typeof(ev3[k])));") 
                    #println(" $(ev1[l] == ev3[k])")
                    # if vars_reversed && !cache.vars_reversed
                    #     reverse!(temp)
                    #     l = get(explookup, temp, -1)
                    #     reverse!(temp)
                    # else 
                        l = get(explookup,temp,-1)
                    # end
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
    
    # if vars_reversed && !cache.vars_reversed
    #     for vec in ev1
    #         reverse!(vec)
    #     end 

    #     for vec in ev2
    #         reverse!(vec)
    #     end 

    #     for vec in ev3
    #         reverse!(vec)
    #     end 

    #     for vec in ev4
    #         reverse!(vec)
    #     end     
    # end 

    return result
end


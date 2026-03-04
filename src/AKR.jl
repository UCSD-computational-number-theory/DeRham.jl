
function reducechain_akr(g,m,f,p,picache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))

    while m > n
        gtemp = picache[m]*g
        len = Int(length(gtemp)/(n+1))
        g = derivative(vector_to_polynomial(gtemp[1:len],n,d*(m-1)-n,PR,params.termorder),1)
        for i in 1:n
            gpolytemp = derivative(vector_to_polynomial(gtemp[(i*len+1):((i+1)*len)],n,d*(m-1)-n,PR,params.termorder),i+1)
            g = g + gpolytemp
        end
        g = polynomial_to_vector(g,n+1,params.termorder)
        if size(g)[1] == 0
            g = fill(R(0), binomial(d*(m-1)-1,n))
        end
        m = m - 1
    end


    return g
end

function reducepoly_akr(pol,S,f,p,picache,params)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    PR = parent(f)
    R = coefficient_ring(parent(f))
    result = PR()
    for term in pol
        #println("Reducing terms of order $(term[2])")
        terms = termsoforder(pol,term[2])
        #println(terms)
        for t in terms
            g = polynomial_to_vector(t[1],n+1,:invlex)
            reduced = reducechain_akr(g,t[2],f,p,picache,params)
            reduced_poly = vector_to_polynomial(reduced,n,d*n-n-1,PR,:invlex)
            result += reduced_poly
        end
    end

    return [[result, n]]
end

function reducetransform_akr(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    M = factor(modulus(base_ring(parent(f))))[7]

    result = similar(FT)

    highm = FT[length(FT)][length(FT[length(FT)])][2]
    
    picache = Vector{typeof(pseudoInverseMat)}(undef, highm)
    picache[n] = pseudoInverseMat
    for i in (n+1):highm
        l = d*i - n - 1
        generate_degree_forward(cache,l)
        generate_degree_forward(cache,l-(d-1))
        generate_degree_forward(cache,l-d)
        pi_new = pseudo_inverse_controlled_lifted(f,S,l,M,params,cache)
        MS = matrix_space(base_ring(parent(f)), nrows(pi_new), ncols(pi_new))
        pseudo_inverse_mat = MS()
        for i in 1:nrows(pi_new)
            for j in 1:ncols(pi_new)
                pseudo_inverse_mat[i,j] = ZZ(pi_new[i,j])
            end
        end
        picache[i] = pseudo_inverse_mat
    end

    for i in 1:length(FT) #pol in FT

        pol = FT[i]
        if (0 < params.verbose)
            println("Reducing vector $i")
            @time reduction = reducepoly_akr(pol,S,f,p,picache,params)
        else
            reduction = reducepoly_akr(pol,S,f,p,picache,params)
        end
        result[i] = reduction

    end

    return result
end
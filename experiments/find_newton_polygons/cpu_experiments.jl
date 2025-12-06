using Combinatorics

"""
For this one, the example is fast, probably
less than 10s,
so we use course-grained parallelism
"""
function cpu_example_fast_random(n,d,p,N,df)

    if d < n
        T = collect(0:d-1)
    else
        T = collect(0:n-1)
    end

    l = ReentrantLock()

    Threads.@threads for i = 1:N
        f = DeRham.random_hypersurface(n,d,p)
        np = DeRham.newton_polygon(f,S=T,fastevaluation=true,algorithm=:naive)
        
        
        if np != false # f is smooth
            #np_key = tuple(np.slopes, np.slopelengths)
            
            @lock l begin 
                update_df(df, n, d, p, np, f)
            end
        end
    end

    return df
end

"""
In this experiment, reducing a single vector is fast,
so we use threading so that each example can be finished in a
shorter period of time.
"""
function cpu_vector_fast_random(n,d,p,N,df)

    if d < n
        T = collect(0:d-1)
    else
        T = collect(0:n-1)
    end


    for i = 1:N
        f = DeRham.random_hypersurface(n,d,p)
        np = DeRham.newton_polygon(f,S=T,fastevaluation=true,algorithm=:naive,use_threads=true)
        
        
        if np != false # f is smooth
            #np_key = tuple(np.slopes, np.slopelengths)
            
            update_df(df, n, d, p, np, f)
        end
    end

    return df
end

# MARK - pasted from MMPSingularities.jl

"""
Calculates if the hypersuface defined by the 
polynomial poly is F-split

note that p must be prime for this to have mathematical meaning
"""
function isFSplit(p,poly)
  #maybe TODO: check that p is prime

  !inPowerOfVariableIdeal(p,p,poly^(p-1))

end#function

"""
Returns true if the polynomial f 
is in the "frobenius power" \\frak{m}^[m],
where \\frak{m} is the ideal of variables of the ring.

"""
function inPowerOfVariableIdeal(p,m,f)
  # don't need this because exponent_vectors will have 
  # no elements for the zero polynomial
  f == zero(f) && return true


  for i in 1:length(f)
    ev = exponent_vector(f,i)

    if all(ev .< m)
      #println("Found term not in the Frob power of the maximal ideal: " * string(exponent_vector))
      
      # We not in the power of the maximal ideal, we don't have any
      # powers that are big enough
      return false
    end
  end

  true
end#function

# MARK - end pasted from MMPSingularities.jl

"""
In this experiment, reducing a single vector is fast,
so we use threading so that each example can be finished in a
shorter period of time.
"""
function cpu_vector_fast_random_K3(n,d,p,N,df)

    if d < n
        T = collect(0:d-1)
    else
        T = collect(0:n-1)
    end


    i = 1
    while i ≤ N
        f = DeRham.random_hypersurface(n,d,p)

        if isFSplit(p,f)
            continue
        end

        println("Found non-F-split example!")
        np = DeRham.newton_polygon(f,S=T,fastevaluation=true,algorithm=:naive,use_threads=true)
        
        println("Completed Zeta function.")
        
        if np != false # f is smooth
            #np_key = tuple(np.slopes, np.slopelengths)
            
            update_df(df, n, d, p, np, f)
        end
        
        i = i + 1
    end

    return df
end


function all_monomials(n,d,p)

    exp_vecs = DeRham.gen_exp_vec(n,d)
    tuples = collect(Iterators.product(0:p-1,exp_vecs))

    tuples[:]
end

"""

R - the polynomial ring which is to be the parent
d - degree
exp_vecs - the output of gen_exp_vec(n,d)
context - an MPolyBuildCtx
"""
function random_monomial(R,p,exp_vecs,context)
    i = rand(1:length(exp_vecs))
    F = base_ring(R)
    push_term!(context,F(rand(0:p-1)),exp_vecs[i])
    finish(context)
end

function cpu_example_fast_example(n,d,p,N,f,df; skip_single_monomials=false)

    R = parent(f)
    F = base_ring(R)

    exp_vecs = DeRham.gen_exp_vec(n,d)

    if d < n
        T = collect(0:d-1)
    else
        T = collect(0:n-1)
    end

    l = ReentrantLock()

    nT = Threads.nthreads()
    ctx_channel = Channel{MPolyBuildCtx}(nT)
    for i in 1:nT
        put!(ctx_channel, MPolyBuildCtx(R))
    end

    function run_example(g)
        np = false 
        # try 
            np = DeRham.newton_polygon(g,S=T,fastevaluation=true,algorithm=:naive)
            
        # catch e

        #     println(g)
        #     throw(e)
        # end
        
        if np != false # f is smooth
            #np_key = tuple(np.slopes, np.slopelengths)
            
            @lock l begin 
                update_df(df, n, d, p, np, g)
            end
        end
    end

    nMons = length(exp_vecs)*p

    # if N is big, just do all of the single monomials
    if nMons < N
        if !skip_single_monomials

            println("N=$N exceeds the number of monomials of degree $d in $n variables. Computing f + mon for all such `mon`")
            Threads.@threads for (coef, exp_vec) in all_monomials(n,d,p)
                ctx = take!(ctx_channel)
                push_term!(ctx,F(coef),exp_vec)
                mon = finish(ctx)
                put!(ctx_channel, ctx)
                ff = f + mon

                run_example(ff)
            end
        end

        if nMons^2 < N
            println("N is more than the number of pairs of monomials. This experiment only considers random pairs of monomials, so this test is suboptimal.")
        end

        println("Randomly trying pairs of monomials.")
        Threads.@threads for i = 1:N
            ctx = take!(ctx_channel)
            mon1 = random_monomial(R,p,exp_vecs,ctx)
            mon2 = random_monomial(R,p,exp_vecs,ctx)
            put!(ctx_channel, ctx)

            ff = f + mon1 + mon2

            run_example(ff)
        end
    else
        
        println("Randomly trying monomials.")
        Threads.@threads for i = 1:N
            ctx = take!(ctx_channel)
            mon = random_monomial(R,p,exp_vecs,ctx)
            put!(ctx_channel, ctx)

            ff = f + mon

            run_example(ff)
        end

    end

    return df
end

function cpu_vector_fast_example(n,d,p,N,f,df)
    R = parent(f)
    F = base_ring(R)

    exp_vecs = DeRham.gen_exp_vec(n,d)

    if d < n
        T = collect(0:d-1)
    else
        T = collect(0:n-1)
    end

    ctx = MPolyBuildCtx(R)

    function run_example(g)
        np = false 

        np = DeRham.newton_polygon(g,S=T,fastevaluation=true,algorithm=:naive,use_threads=true)
            
        if np != false # f is smooth
            update_df(df, n, d, p, np, g)
        end
    end

    nMons = length(exp_vecs)*p

    # if N is big, just do all of the single monomials
    if nMons < N
        if !skip_single_monomials

            println("N=$N exceeds the number of monomials of degree $d in $n variables. Computing f + mon for all such `mon`")
            for (coef, exp_vec) in all_monomials(n,d,p)
                push_term!(ctx,F(coef),exp_vec)
                mon = finish(ctx)
                ff = f + mon

                run_example(ff)
            end
        end

        if nMons^2 < N
            println("N is more than the number of pairs of monomials. This experiment only considers random pairs of monomials, so this test is suboptimal.")
        end

        println("Randomly trying pairs of monomials.")
        for i = 1:N
            mon1 = random_monomial(R,p,exp_vecs,ctx)
            mon2 = random_monomial(R,p,exp_vecs,ctx)

            ff = f + mon1 + mon2

            run_example(ff)
        end
    else
        println("Randomly trying monomials.")
        for i = 1:N
            mon = random_monomial(R,p,exp_vecs,ctx)

            ff = f + mon

            run_example(ff)
        end

    end

    return df
end

"""
Takes an example, randomly weights 
"""
function cpu_vector_fast_weighted_example(n,d,p,N,f,df; num_monomials=1,nonzero_weights=false)
    R = parent(f)
    F = base_ring(R)

    exp_vecs = DeRham.gen_exp_vec(n,d)

    if d < n
        T = collect(0:d-1)
    else
        T = collect(0:n-1)
    end

    ctx = MPolyBuildCtx(R)

    function run_example(g)
        np = false 

        np = DeRham.newton_polygon(g,S=T,fastevaluation=true,algorithm=:naive,use_threads=true)
            
        if np != false # f is smooth
            update_df(df, n, d, p, np, g)
        end
    end

    nMons = length(exp_vecs)*p
    nWeights = p^length(terms(f))
    nExamples = nMons*nWeights

    # if N is big, just do all of the single monomials
    if nExamples < N

        println("Warning: N=$N exceeds the number of monomial-weight combinations of degree $d in $n variables. Consider using a different test that enumerates these combinations.")
    end

    println("Randomly scaling monomials of f by weights and adding monomials.")
    for i = 1:N

        ts = terms(f)

        if nonzero_weights
            weights = rand(1:p-1,length(ts))
        else
            weights = rand(0:p-1,length(ts))
        end

        ff = sum(weights .* ts)

        for i in 1:num_monomials
            mon = random_monomial(R,p,exp_vecs,ctx)
            ff = f + mon
        end

        run_example(ff)
    end

    return df
end

"""
Returns the integer partitions of n
with less than or equal to k parts.

The result is returned as a vector of vectors,
with each inner vector having length k.
Each inner vector is padded with zeros so it has length k
"""
function partitions_leq(n,k)

    function padwithzeros!(partition,size) 
        padding_length = size - length(partition)
        if 0 < padding_length
            append!(partition,zeros(eltype(partition), padding_length))
        end
    end

    result = collect.(Oscar.partitions(n,1))# = [[n]]
    padwithzeros!.(result,k)

    for i in 2:k
        new_partitions = collect.(Oscar.partitions(n,i))
        padwithzeros!.(new_partitions,k)
        result = vcat(result,new_partitions)
    end
    
    result
end

"""
TODO: PR this to Oscar? I don't think Oscar has such a method...
"""
function monomial_symmetric_polynomial(R,partition; ctx=nothing)
    d = sum(partition)
    vars = gens(R)
    l = length(partition)
    nDistinct = length(unique(partition))
    n = length(vars)

    n != l && throw("Need $l variables for partition $partition but have $n.")

    if ctx == nothing
        ctx = MPolyBuildCtx(R)
    end

    for p in Combinatorics.multiset_permutations(partition,l)
        push_term!(ctx, one(base_ring(R)), p)
    end

    finish(ctx)
end
        
function all_monomial_symmetric_polynomials(n,d,p)
    R, _ = polynomial_ring(GF(p),n)

    partitions = partitions_leq(d,n)

    ctx = MPolyBuildCtx(R)
    
    monomial_symmetric_polynomial.((R,),partitions,ctx=ctx)
end

function all_tuples_of_length(n,p)
    collect(Iterators.product(fill(0:p-1,n)...))[:]
end

"""
Takes random linear combinations of symmetric polynomials and computes
their newton polygons
"""
function cpu_vector_fast_symmetric(n,d,p,N,df)

    symmpolys = all_monomial_symmetric_polynomials(n,d,p)

    # R = parent(symmpolys[1])

    if d < n
        T = collect(0:d-1)
    else
        T = collect(0:n-1)
    end

    function run_example(g)
        np = false 

        np = DeRham.newton_polygon(g,S=T,fastevaluation=true,algorithm=:naive,use_threads=true)
            
        if np != false # f is smooth
            update_df(df, n, d, p, np, g)
        end
    end



    # formula comes from p^l possible vectors, but scaling
    # by an element in the base gives the same surface
    # So we fix the first coordiate to be either 1 
    # (p^(l-1) possibilities) or 0 (p^(l-1) possibilites)
    l = length(symmpolys)
    nExamples = 2*p^(l - 1)

    # if N is big, just do all of the single monomials
    if nExamples < N

        println("N=$N exceeds the number of symmetric polynomials of degree $d in $n variables. Computing the newton polygon for all such examples")

        for weights in all_tuples_of_length(l-1,p)
            interesting_term = sum(weights .* symmpolys[2:end])

            run_example(interesting_term)
            run_example(symmpolys[1] + ineresting_term)
        end

    else
        println("Randomly trying linear combinations of symmetric polynomials.")
        for i = 1:N
            weights = rand(0:p-1,l)

            ff = sum(weights .* symmpolys)

            run_example(ff)
        end

    end

    return df
end

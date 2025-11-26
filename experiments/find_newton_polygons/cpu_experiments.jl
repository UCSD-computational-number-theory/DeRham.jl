
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

function cpu_vector_fast_random(n,d,p,N)

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

function cpu_vector_fast_fermat(n,d,p,N)

end


function is_Ssmooth(f,S,params)
    p = characteristic(parent(f))
    n = length(gens(parent(f))) - 1
    d = total_degree(f)
    l = d * n - n + d - length(S)

    PR = parent(f)
    R = base_ring(PR)

    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    fLift = liftCoefficients(ZZ,PRZZ,f)

    params = default_params()
    cache = controlled_reduction_cache(n,d,S,params)

    M_ZZ = compute_controlled_matrix(fLift,l,S,ZZ,PRZZ,params,cache)

    #M = zeros(Float64,size(M_ZZ)...)
    #float_entries!(M,M_ZZ)
    #M_gpu = CuModMatrix{Float64}(M,Int(p))
    #flag, B = GPUFiniteFieldMatrices.is_invertible_with_inverse(M_gpu)
    
    M = matrix(R,M_ZZ)
    if params.use_gpu
        flag, B = GPUFiniteFieldMatrices.is_invertible_with_inverse(
            CuModMatrix(M, N=p, rows=rows(M), cols=cols(M))
        )
    else
        flag, B = Nemo.is_invertible_with_inverse(M, side=:right)
    end
    #is_invertible(M; side=:right)
    flag
end

function issmooth_linalg(f)
    n = length(gens(parent(f))) - 1
    is_Ssmooth(f,collect(0:n))
end

function random_change_of_variables(f)

    vars = gens(parent(f))
    n = length(vars)-1 # we don't need it here, but use the standard convention
    p = characteristic(parent(f))

    SLn = special_linear_group(n + 1, GF(p))
    mat = rand(SLn)
    new_vars = matrix(mat) * vars
    f_transformed = evaluate(f, new_vars)

    f_transformed
end

"""
    find_Ssmooth_model(f, M, S_target, params)
Find a change of variable that makes f S_target-smooth 

Inputs: 
* "f" -- Oscar polynomial (should be homogeneous)
* "M" -- integer, algorithm precision 
* "S_target" -- list, target S-smoothness 
"""
function find_Ssmooth_model(f, M, S_target, params, changef, cache)
    p = Int64(characteristic(parent(f)))
    q = p
    PR = parent(f)
    R = coefficient_ring(PR)
    n = nvars(PR) - 1
    SLn = special_linear_group(n + 1, GF(q))
    d = total_degree(f)
    vars = gens(PR)
    l = d * n - n + d - length(S_target)

    #bool = true 
    f_transformed = f
    #println(f_transformed)
    num_iter = min(2*p,100)
    f_changed = false 
    if !changef
        pseudo_inverse_mat_new = pseudo_inverse_controlled_lifted(f_transformed,S_target,l,M,params,cache)
            #bool = false
            
        return f_changed, f_transformed, pseudo_inverse_mat_new
    end

    for i in 1:num_iter 
        try 
            pseudo_inverse_mat_new = pseudo_inverse_controlled_lifted(f_transformed,S_target,l,M,params,cache)
            #bool = false
            
            return f_changed, f_transformed, pseudo_inverse_mat_new
        catch e
            if !changef
                throw(ArgumentError("f is not $S_target smooth"))
                return false, false, false 
            end 
            (0 < params.verbose) && println("This f is not S-smooth, changing to one that is")
            if isa(e, ArgumentError) && e.msg == "f is not smooth"
                #throw(ArgumentError("f is not smooth"))
                #bool = false  
                return false, false, false  
            else
                mat = rand(SLn)
                new_vars = matrix(mat) * vars
                f_transformed = evaluate(f, new_vars)
                f_changed = true 
            end
        end 
    end
    return (false,false,false)
    
end 

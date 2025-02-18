"""
    find_Ssmooth_model(f, M, S_target, params)
Find a change of variable that makes f S_target-smooth 

Inputs: 
* "f" -- Oscar polynomial (should be homogeneous)
* "M" -- integer, algorithm precision 
* "S_target" -- list, target S-smoothness 
"""
function find_Ssmooth_model(f, M, S_target, params)
    p = Int64(characteristic(parent(f)))
    q = p
    PR = parent(f)
    R = coefficient_ring(PR)
    n = nvars(PR) - 1
    SLn = special_linear_group(n + 1, GF(q))
    d = total_degree(f)
    vars = gens(PR)
    l = l = d * n - n + d - length(S_target)
    
    bool = true 
    while bool
        mat = rand(SLn)
        new_vars = mat * vars
        f_transformed = evaluate(f, new_vars)
        try 
            pseudo_inverse_mat_new = pseudo_inverse_controlled_lifted(f_transformed,S_target,l,M,params)
            bool = false
            return f_transformed, pseudo_inverse_mat_new
        catch e
        end 
    end 
    
end 

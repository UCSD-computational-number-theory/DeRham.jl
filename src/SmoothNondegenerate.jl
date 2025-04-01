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
    l = d * n - n + d - length(S_target)
    
    #TODO: fix the anti-pattern below and remove this band-aid
    #if !check_smoothness(f)
    #    error("f is not smooth")
    #end

    #bool = true 
    f_transformed = f
    #while bool
    num_iter = 100
    for i in 1:num_iter 
        try 
            pseudo_inverse_mat_new = pseudo_inverse_controlled_lifted(f_transformed,S_target,l,M,params)
            bool = false
            return f_transformed, pseudo_inverse_mat_new
        catch e
            (0 < verbose) && println("This f is not S smooth, changing to one that is")
            if isa(e, ArgumentError) && e.msg == "f is not smooth"
                throw(ArgumentError("f is not smooth"))
                bool = false  
                return false 
            else
                mat = rand(SLn)
                new_vars = matrix(mat) * vars
                f_transformed = evaluate(f, new_vars)
            end
        end 
    end 
    
end 

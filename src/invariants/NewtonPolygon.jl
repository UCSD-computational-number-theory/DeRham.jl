"""
    newton_polygon(f)

    Computes the newton polygon of f by first computing its zeta function
"""
function newton_polygon(f; S=[-1], verbose=0, changef=true, algorithm=:default, termorder=:invlex, vars_reversed=false, fastevaluation=true, always_use_bigints=false, use_gpu=false)
    PR = parent(f)
    R = coefficient_ring(PR)
    p = Int64(characteristic(PR))
    zf = zeta_coefficients(f; S=S, verbose=verbose, changef=changef, algorithm=algorithm, termorder=termorder, vars_reversed=vars_reversed, fastevaluation=fastevaluation, always_use_bigints=always_use_bigints, use_gpu=use_gpu)
    
    if zf == false # f isn't smooth
        return false
    end 

    return newton_polygon(p, zf)
end 

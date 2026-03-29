"""
    artin_mazur_height
    computes the Artin-Mazur height of a quartic K3 surface given the coefficients of its L-polynomial 
"""
function artin_mazur_height(coeffs, q)
    @assert length(coeffs) == 22
    if (coeffs[1] == q) && (coeffs[22] == q) # if the L-polynomial is in the boat shape format
        for i in 2:11
            if mod(coeffs[i], q) != 0
                return i-1
            end 
        end 
        return Inf
    else
        boatshape_coeffs = boat_shape_Lpoly(coeffs, 21, q)
        for i in 2:11
            if mod(boatshape_coeffs[i], q) != 0
                return i-1
            end 
        end 
        return Inf
    end 
end
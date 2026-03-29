"""
    LPolynomial

    converts the output of zeta_coefficients to an Oscar polynomial
"""
function LPolynomial(zeta_coeffs)
    P, T = polynomial_ring(ZZ, "T")
    deg = length(zeta_coeffs) - 1
    return sum(T^(deg+1-i) * ZZ(zeta_coeffs[i]) for i in 1:deg+1)
end 

"""
    boat_shape_Lpoly(coeffs, deg, q)
Converts the coefficients of the L-polynomial into the boat-shaped L-polynomial, via the change of variable q*P(T/q)
"""
function boat_shape_Lpoly(coeffs, deg, q)
    n = deg + 1
    coeffs_new = [ZZ(0) for i in 1:n]
    for i in 1:n
        coeffs_new[i] = div(q*ZZ(coeffs[i]), ZZ(q)^(n-i))
    end 
    
    return coeffs_new
end 
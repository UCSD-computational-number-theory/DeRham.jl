"""
    prec_vec(polygon, relative_precision)

returns the vector by adding the height of the Hodge polygon and appropriate relative precision for each position

INPUTS: 
* "polygon" -- A SlopesPolygon struct describing a polygon whose values describe the
divisibility of the roots.
* "relative_precision" -- list, list of relative precisions, can be computed by "calculate_relative_precision"
"""
function prec_vec(polygon, relative_precision)
    vals = Int.(polygon.values)
    hodge_numbers = polygon.slopelengths
    prec_vec = reverse(vals)
    i = 1
    num_term = 1 
    for j in 1:length(hodge_numbers)
        num_term = num_term + hodge_numbers[j]
        while i < num_term
            prec_vec[i] = prec_vec[i] + relative_precision[j]
            i = i + 1
        end 
    end 
    prec_vec[length(prec_vec)] = prec_vec[length(prec_vec)] + relative_precision[length(hodge_numbers)]

    return prec_vec 
end 

"""
    sign_fe(n, d, prec_vec, cp_coeffs)

returns the sign of the functional equation. The sign is always 1 if the dimension of the hypersurface is odd. 

INPUTS: 
* "n" -- integer, dimension of the ambient projective space 
* "d" -- integer, degree of characteristic polygon = dimension of cohomology group 
* "prec_vec" -- list, prec_vec[i] = maximum p-adic valuation of the i-th coefficient of the Frobenius char poly? 
* "cp_coeffs" -- list, cp_coeffs[i] = i-th coefficient of the Frobenius char poly, which a priori may be incorrect 
"""
function sign_fe(n, d, prec_vec, cp_coeffs)
    sign = 1
    dimension = n-1 # dimension of hypersurface = motivic weight = n-1 
    if dimension % 2 == 0 
        for i in 1:Int(floor(d/2))
            # NOTE-TO-SELF: check if the symmetric index is d-i or d+1-i
            p_power = p^(min(prec_vec[i], prec_vec[d-i]+((d-2*i)*dimension)//2))
            if (mod(cp[i], p_power) != 0) && (mod(cp[d-i], p_pwer) != 0)
                if 0 == mod(cp[i] + cp[d-i] * p^(((d-2*i)*dimension)//2), p_power)
                    sign = -1
                else
                    sign = 1 
                    @assert 0 == mod(cp[i] - cp[d-i] * p^(((d-2*i)*dimension)//2), p_power)
                end 
            end 
        end 
    end 

    return sign 
end 

"""
    apply_symmetry(n, d, p, prec_vec, cp_coeffs, mod, sign)

Update the coefficients in the Frobenius characteristic polynomial using symmetry in the coefficients coming from 
the Weil conjecture/Poincare duality 

INPUTS: 
* "n" -- integer, dimension of the ambient projective space 
* "d" -- integer, degree of characteristic polygon = dimension of cohomology group 
* "p" -- integer, prime number 
* "prec_vec" -- list, prec_vec[i] = maximum p-adic valuation of the i-th coefficient of the Frobenius char poly? 
* "cp_coeffs" -- list, cp_coeffs[i] = i-th coefficient of the Frobenius char poly, which a priori may be incorrect 
* "mod" -- list, mod[i] = p^(prec_vec[i])
* "sign" -- integer, sign of the functional equation
"""
function apply_symmetry(n, d, p, prec_vec, cp_coeffs, mod, sign)
    dimension = n-1
    for i in 1:(div(Pdeg, 2) + 1)
        k = ((d-2*i) * dimension)//2 
        if prec_vec[i] >= prec_vec[d-i] + k 
            prec_vec[d-i] = prec_vec[i] - k 
            mod[d-i] = p^(prec_vec[d-i])
            cp_coeffs[d-i] = (sign * cp_coeffs[i])//p^k  # ?
            cp_coeffs[d-i] = mod(cp_coeffs[d-i], mod[d-i])
        else
            prec_vec[i] = prec[d-i] + k 
            mod[i] = p^(prec_vec[i])
            cp_coeffs[i] = sign * cp_coeffs[d-i] * p^k
            cp_coeffs[i] = mod(cp_coeffs[i], mod[i]) 
        end 
    end 

    return cp_coeffs 
end 
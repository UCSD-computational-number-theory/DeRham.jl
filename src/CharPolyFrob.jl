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
    cp = cp_coeffs #TODO: why the naming discrepancy?
    if dimension % 2 == 0 
        for i in 2:Int(floor(d/2))
            # NOTE-TO-SELF: check if the symmetric index is d-i or d+1-i
            k = div((d-2*(i-1))*dimension,2)
            p_power = p^(min(prec_vec[i], prec_vec[d+1-i]+k))
            if (mod(cp[i], p_power) != 0) && (mod(cp[d+1-i], p_power) != 0)
                if 0 == mod(cp[i] + cp[d+1-i] * p^k, p_power)
                    sign = -1
                else
                    sign = 1 
                    @assert 0 == mod(cp[i] - cp[d+1-i] * p^k, p_power)
                end 
            end 
        end 
    end 

    return sign 
end 

"""
    apply_symmetry(n, d, p, prec_vec, cp_coeffs, modulus, sign)

Update the coefficients in the Frobenius characteristic polynomial using symmetry in the coefficients coming from 
the Weil conjecture/Poincare duality 

INPUTS: 
* "n" -- integer, dimension of the ambient projective space 
* "d" -- integer, degree of characteristic polygon = dimension of cohomology group 
* "p" -- integer, prime number 
* "prec_vec" -- list, prec_vec[i] = maximum p-adic valuation of the i-th coefficient of the Frobenius char poly? 
* "cp_coeffs" -- list, cp_coeffs[i] = i-th coefficient of the Frobenius char poly, which a priori may be incorrect 
* "modulus" -- list, modulus[i] = p^(prec_vec[i])
* "sign" -- integer, sign of the functional equation
"""
function apply_symmetry(n, d, p, prec_vec, cp_coeffs, modulus, sign)
    dimension = n-1
    for i in 1:div(d, 2)+1
        k = div((d-2*(i-1)) * dimension, 2) 
        if prec_vec[i] >= (prec_vec[d+2-i] + k) 
            prec_vec[d+2-i] = prec_vec[i] - k 
            modulus[d+2-i] = p^(prec_vec[d+2-i])
            cp_coeffs[d+2-i] = div(sign * cp_coeffs[i], p^k)  
            cp_coeffs[d+2-i] = mod(cp_coeffs[d+2-i], modulus[d+2-i])
        else
            prec_vec[i] = prec_vec[d+2-i] + k 
            modulus[i] = p^(prec_vec[i])
            cp_coeffs[i] = sign * cp_coeffs[d+2-i] * p^k
            cp_coeffs[i] = mod(cp_coeffs[i], modulus[i]) 
        end 
    end 
    #println("modulus=$modulus")
    return cp_coeffs 
end 

"""
    apply_newton_identity(n, d, p, prec_vec, cp_coeffs, modulus)
calculate the i-th power sum of the roots and correct the coefficients along the way 

INPUTS: 
* "n" -- integer, dimension of the ambient projective space 
* "d" -- integer, degree of characteristic polygon = dimension of cohomology group 
* "p" -- integer, prime number 
* "prec_vec" -- list, prec_vec[i] = maximum p-adic valuation of the i-th coefficient of the Frobenius char poly? 
* "cp_coeffs" -- list, cp_coeffs[i] = i-th coefficient of the Frobenius char poly, which a priori may be incorrect 
* "modulus" -- list, modulus[i] = p^(prec_vec[i])
"""
function apply_newton_identity(n, d, p, prec_vec, cp_coeffs, modulus)
    e = [0 for i in 1:d+1]  # e[k] = (k-1)-th elementary symmetric polynomial evaluated at the roots
    s = [0 for i in 1:d+1]  # s[k] = k-th power sum of the roots 
    for k in 1:d+1 
        if (k-1)%2 == 0
            e[k] = cp_coeffs[d+2-k]  # FIXME: check the index 
        else
            e[k] = -cp_coeffs[d+2-k]
        end
    end 

    for k in 1:d-1
        sum = 0 
        for i in 1:k-1
            if i%2 == 1
                sum = sum + e[k+1-i] * s[i+1]
            else
                sum = sum - e[k+1-i] * s[i+1]
            end  
        end 
        if k%2 == 0
            s[k+1] = sum - k * e[k+1]
        else
            s[k+1] = -1 * (sum - k * e[k+1])
        end 
        if k == 12
            #println("k=12, sum=$sum")
        end 

        pN = k * modulus[d+1-k]
        s[k+1] = mod(s[k+1], pN)
        #println(s[k+1])

        if ZZ(s[k+1]^2) > d*d*ZZ(p^((n-1)*k))
            s[k+1] = -mod(-s[k+1], pN)
        end 
        #println(s[k+1])

        if k%2 == 0
            e[k+1] = div(sum-s[k+1], k)
            cp_coeffs[d+1-k] = e[k+1]
        else
            e[k+1] = div(sum+s[k+1], k)
            cp_coeffs[d+1-k] = -e[k+1]
        end
        #println("ek=$e")
        #println("sk=$s")
    end 
    return cp_coeffs
end 

"""
    compute_Lpolynomial(n, d, p, polygon, relative_precision, cp_coeffs)
INPUTS: 
* "n" -- integer, dimension of the ambient projective space 
* "p" -- integer, prime number 
* "polygon" -- A SlopesPolygon struct describing a polygon whose values describe the
divisibility of the roots.
* "relative_precision" -- list, list of relative precisions, can be computed by "calculate_relative_precision"
* "cp_coeffs" -- list, list of integers representing the coefficient of the characteristic polynomial of the 
                        Frobenius matrix computed over the integers 

"""
function compute_Lpolynomial(n, p, polygon, relative_precision, cp_coeffs)
    dimension = n - 1  # dimension of the projective space  
    Lpoly_coeffs = cp_coeffs
    #println("initial coefficients = $Lpoly_coeffs")
    prec_vec = DeRham.prec_vec(polygon, relative_precision)
    #println("prec_vec = $prec_vec")
    d = length(cp_coeffs) - 1  # degree of characteristic polynomial 
    modulus = [0 for i in 1:d+1]
    for i in 1:d+1
        modulus[i] = ZZ(p^prec_vec[i])
        Lpoly_coeffs[i] = mod(Lpoly_coeffs[i], modulus[i])
    end 
    #println("modulus=$modulus")
    #println("coefficients after moding by prec = $Lpoly_coeffs")

    Lpoly_coeffs[d+1] = 1  # set leading coefficient to 1 
    sign = sign_fe(n, d, prec_vec, Lpoly_coeffs)  # determine the sign of the functional equation 
    #println("sign of functional equation = $sign")
    Lpoly_coeffs[1] = sign * ZZ(p^((d*dimension)/2))
    Lpoly_coeffs = DeRham.apply_symmetry(n,d,p,prec_vec,Lpoly_coeffs,modulus,sign)
    #println("coefficients after apply symmetry = $Lpoly_coeffs")
    #println("modulus=$modulus")
    Lpoly_coeffs = DeRham.apply_newton_identity(n,d,p,prec_vec,Lpoly_coeffs,modulus)
    #println("coeffficients after applying Newton's identity = $Lpoly_coeffs")

    return Lpoly_coeffs 
end


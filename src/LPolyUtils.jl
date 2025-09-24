"""
    vp(a, p)
    Returns the p-adic valuation of a 
    
    INPUTS: 
    * "a" -- integer
    * "p" -- integer, a prime number 
"""
function vp(a, p)
    @assert is_prime(p)
    if a == 0 
        return Inf 
    end 

    e = 0 
    while mod(a, p) == 0
        a = div(a, p)
        e = e + 1
    end

    return e
end 

"""
    newton_polygon(f, p)
    Returns the Newton polygon of f 

    INPUTS: 
    * "f" -- list, coefficients of a polynomial, constant term first
    * "p" -- integer, a prime number
"""
function newton_polygon(f, p)
    @assert is_prime(p)

    pts = []
    for (i, a) in enumerate(f)
        if a != 0
            push!(pts, (i-1, vp(a, p)))
        end 
    end 

    # compute the lower convex hull of points
    hull = Tuple{Int64, Int64}[]
    for pt in pts
        while length(hull) >= 2
            (x1, y1) = hull[end-1]
            (x2, y2) = hull[end]
            (x3, y3) = pt
            # Compute cross product to decide whether pt is "below" the last segment.
            cp = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
            if cp <= 0
                pop!(hull)
            else
                break
            end
        end
        push!(hull, pt)
    end
    return SlopesPolygon(hull)
   #return hull
end 

"""
zeta_to_pointcount(Z,n)

    Returns the first n point counts given the zeta function Z
"""
function zeta_to_pointcount(Z, n)
    zeta_log = log(Z)
    return [coeff(zeta_log, i) * i for i in 1:n]
end 

"""
k3boat_shape_Lpoly_to_pointcount(coeffs, q, n)
    Returns the first n point counts (up to X(F_q^n)) given the boat-shaped L-polynomial of a K3 surface over F_q 
    Z(X,T) = 1/(1-T)(1-qT)(1-q^2T)q^(-1)L(qT) 

    INPUTS: 
    * "coeffs" -- list, the coefficients of the boat-shaped L-polynomial 
    * "q" -- integer, the cardinality of the base field 
    * "n" -- integer
"""
function k3boat_shape_Lpoly_to_pointcount(coeffs, q, n)
    @assert length(coeffs) == 22
    @assert (coeffs[1] == q) and (coeffs[22] == q) 

    R,T = power_series_ring(QQ, max(n+10, 25), :T)
    L = sum([coeffs[i] * (q*T)^(i-1) for i in 1:22])
    zeta = q/((1-T) * (1-q*T) * (1-q^2*T) * L)
    return zeta_to_pointcount(zeta, n)
end

"""
k3_Lpoly_to_pointcount(coeffs, q, n)
    Returns the first n point counts (up to X(F_q^n)) given the L-polynomial of a K3 surface over F_q 
    Z(X,T) = 1/(1-T)(1-qT)(1-q^2T)L(T) 

    INPUTS: 
    * "coeffs" -- list, the coefficients of the boat-shaped L-polynomial 
    * "q" -- integer, the cardinality of the base field 
    * "n" -- integer
"""
function k3_Lpoly_to_pointcount(coeffs, q, n)
    @assert length(coeffs) == 22
    @assert coeffs[22] == 1

    R,T = power_series_ring(QQ, max(n+10, 25), :T)
    L = sum([coeffs[i] * T^(22-i) for i in 1:22])
    zeta = 1/((1-T) * (1-q*T) * (1-q^2*T) * L)
    return zeta_to_pointcount(zeta, n)
end


"""
curve_Lpoly_to_pointcount(coeffs, q, n)
    Returns the first n point counts (up to X(F_q^n)) given the L-polynomial of a curve X 
    Z(X,T) = L/(1-T)(1-qT)

    INPUTS: 
    * "coeffs" -- list, the coefficients of the boat-shaped L-polynomial 
    * "q" -- integer, the cardinality of the base field 
    * "n" -- integer
"""
function curve_Lpoly_to_pointcount(coeffs, q, n)
    R,T = power_series_ring(QQ, max(n+10, length(coeffs)+1), :T)
    d = length(coeffs)
    L = sum([coeffs[i] * T^(d-i) for i in 1:d])
    zeta = L/((1-T) * (1-q*T))
    return zeta_to_pointcount(zeta, n)
end 


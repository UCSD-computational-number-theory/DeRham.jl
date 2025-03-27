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




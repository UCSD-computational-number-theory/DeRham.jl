"""
    hodgepolygon(basis::Array)

Calculates the hodge polygon of the cohomology module with 
griffiths-dwork basis basis

basis -- an array of "polynomials with pole" as descirbed in PolynomialWithPole.jl
"""
function hodgepolygon(basis::Vector, n)
    #WRONG: n = highestpoleorder(basis)
    hodgenumbers = zeros(Int,n)
    for i in 0:n-1
        h = length(termsoforder(basis,n-i))
        hodgenumbers[i+1] = h
    end

    SlopesPolygon(hodgenumbers)
end

"""
    hodgepolygon(f; basis=nothing, params=default_params())

Calculates the hodge polygon of f

f - the polynomial to get the hodge polygon of
"""
function hodgepolygon(f::RingElem; basis=nothing, params=default_params())
    n = nvars(parent(f)) - 1
    PR = parent(f)
    R = coefficient_ring(parent(f))

    if basis == nothing
        d = total_degree(f)
        S = collect(0:n) #collect(0:d-1)
        cache = controlled_reduction_cache(n,d,S,params)
        basis = compute_monomial_bases(f, params, cache) # basis of cohomology 
    end
    
    Basis = []
    for i in 1:n
        for j in basis[i]
            push!(Basis,[j,i])
        end
    end

    hodgepolygon(Basis,n)
end
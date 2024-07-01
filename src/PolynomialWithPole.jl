module PolynomialWithPole

using Oscar

"""
A polynomial with (a) pole is mathematically
an element of k[x_1, ..., x_n]_f, the localization
of the polynomial right at an element f. 
For us, f is predetermined.

Such an element is determined by the order of the pole
and the numerator, i.e. a polynomial.

Right now, this is modeled in our code by 
arrays of length 2. 
The first term is an oscar polynomial and 
the second is an integer, the pole order.

Another common type is an array of polynomials
with pole, ordered by pole order. 

I suspect in the future we may change this to
an array of tuples or a custom struct. 
If/when that happens, hopefully all of the necessary 
changes will be in this file.

NOTE: these typealiases aren't used anywhere yet, 
it mostly exists so that this docstring can exist
right now.
"""
const PolyWithPole = Array{Any}
const PolysWithPole = Array{Array{Any}}

"""
gets all of the terms of polpol, 
including their pole information

polpol - a PolynomialWithPole
"""
function terms(polpol)
    t = Oscar.terms(polpol[1])
    result = []
    for i in t
        push!(result,[i,polpol[2]])
    end
    return result
end

"""
gets all of the terms of polpol,
forgetting their pole information
"""
numeratorterms(polpol) = Oscar.terms(polpol[1])


"""
Get the PolynomialWithPole which has the highest
pole order in the list

polpols - array of polynomials with pole
"""
function highestpoleorder(polpols)
    function cmp_poleorder(p1,p2)
        if p1[2] ≤ p2[2]
            p2
        else
            p1
        end
    end 
    
    reduce(cmp_poleorder,polpols)[2]
end

"""
Gets all terms of a certain pole order in an
array of polynomials with pole

polpols - array of polynomials with pole
"""
function termsoforder(polpols,order)
    result = []
    for polpol in polpols
        if polpol[2] == order
            append!(result,terms(polpol)) 
        end
    end

    result
end

end

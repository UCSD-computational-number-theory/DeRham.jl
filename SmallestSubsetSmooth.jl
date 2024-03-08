module SmallestSubsetSmooth

using Oscar
using Combinatorics
using Printf

#=
n = 2
d = 3
p = 7
precision = 30
# k = GF(p)

R, vars = polynomial_ring(PadicField(p, precision), ["x$i" for i in 0:n])
(x,y,z) = vars
polynomial = z * y^2 - x^3 - x^2 * z
=#

function smallest_subset_s_smooth(polynomial,n)
    R = parent(polynomial)
    vars = gens(R)
    for k in 0:n
        for combo in combinations(0:n, k)
            combo_c = setdiff(collect(0:n), combo)
            I = ideal(R, [derivative(polynomial, i + 1) for i in combo])
            J = ideal(R, [vars[i + 1] * derivative(polynomial, i + 1) for i in combo_c])
            jacobian = I + J
    
            entire = ideal(R, [vars[i + 1] for i in 0:n])
            if radical(jacobian) == entire
                return combo
            end
        end
    end
    return nothing
end

#=
res = smallest_subset_s_smooth()
if isnothing(res)
    @printf "Nothing satisfies!"
else
    @printf "Found answer of length %i:\n" length(res)
    return res
end
=#

end
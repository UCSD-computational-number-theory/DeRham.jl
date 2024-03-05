module Utils

using Oscar

macro assert(ex)
    return :( $ex ? nothing : throw(AssertionError($(string(ex)))) )
end


# Convert polynomial to vector form with specified order. Default is lexicographic.
function polynomial_to_vector(f, n, R, PR; order=:lex)
    vars = gens(PR)

    if order == :lex
        d = total_degree(f)

        mon = compute_monomials(n, d, PR,vars)
        res = fill(R(0), length(mon))
        for i in eachindex(mon)
            res[i] = coeff(f, mon[i])
        end
        return res
    else
        throw(ArgumentError("Invalid option '$order'"))
    end
end

# Convert vector to polynomial with specified order. Default is lexicographic.
function vector_to_polynomial(vect, n, d, PR, order=:lex)
    if order == :lex
        res = PR()
        mon = compute_monomials(n + 1, d, PR)
        for i in eachindex(vect)
            res += PR(vect[i]) * mon[i]
        end
    else
        throw(ArgumentError("Invalid option '$order'"))
    end
    return res
end

# Given an Oscar Matrix M, return a the indices of those columns whcih are pivots after being particular
# into reduced row echelon form.
function pivot_columns(M)
    rank, M = rref(M)
    res = fill(0, rank)
    ncols = size(M, 2)
    j = 1
    for i in 1:rank
        while j <= ncols && M[i, j] == 0
            j += 1
        end
        res[i] = j
        j+=1
    end
    return res
end

# Computes all monomials of degree `d` in `n` variables in the polynomial ring `PR`.
function compute_monomials(n, d, PR, vars, order=:lex)
    if d <= 0
        return []
    end

    if order == :lex
        result = []
        
        function backtrack(start, current)
            if length(current) == d
                push!(result, prod(PR(var) for var in current))
                return
            end

            for i in start:n
                backtrack(i, [current..., vars[i]])
            end
        end

        backtrack(1, [])

        return result
    else
        throw(ArgumentError("Invalid option '$order'"))
    end
end

end
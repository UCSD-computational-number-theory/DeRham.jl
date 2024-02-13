module PrecisionEstimate

p = 7
r = 12 # precision
n = 4

function compute_N(p, r, m, n)
    num = n * log(n + r)
    denom = (n + r) * log(p) - n
    return ceil(Int, -m + (n + r) * (1 + num / denom)) # Convert the result to an integer
end

function compute_precisions_each(p, r, n)
    for i = 1:n
        println(compute_N(p, r, i, n)) # Change Printf to println
    end
end

compute_precisions_each(p, r, n)
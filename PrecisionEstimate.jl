module PrecisionEstimate

p = 41
r = 7 # precision
n = 2

function compute_N(p, r, m, n)
    num = n * log(n + r)
    denom = (n + r) * log(p) - n
    return ceil(Int, -m + (n + r) * (1 + num / denom))
end

function compute_precisions_each(p, r, n)
    for m = 1:n
        println(compute_N(p, r, m, n))
    end
end

compute_precisions_each(p, r, n)
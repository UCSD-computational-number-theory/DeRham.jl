
"""
    controlled_reduction_cache

Give the minimal PolyExpCache needed for controlled reduction
as in Prop 1.15 of Costa's thesis.

n in this function is the number of variables minus one,
i.e. the weight of the affine Monsky-Washnitzer homology module.

the d we need are

FORWARD:
h*d - n - 1 for h in 1:n (if such values are positive)
d*n - n 
d*n - n - 1
d*n - n - length(S),
d*n - n - length(S) + 1,
d*n - n + d - length(S) + 1, #TODO: is this +1 wrong??
d

REVERSE:
d
n*d - n

fields
------
n - int
d - int
S - vector of ints
params - the ControlledReductionParamaters
"""
function controlled_reduction_cache(n,d,S,params)
    degsforward = [d*n - n,
                   d*n - n - 1,
                   d*n - n - length(S),
                   d*n - n - length(S) + 1,
                   d*n - n + d - length(S),
                   d]
    for h in 1:n
        if 0 ≤ h*d - n - 1
            append!(degsforward,h*d - n - 1)
        end
        if 0 ≤ h*d - n - d
            append!(degsforward,h*d - n - d)
        end
    end

    degsreverse = [d,
                   n*d - n,
                   d*n - n + d - length(S)]

    PolyExpCache(n+1,
                 params.termorder,
                 degsforward,
                 degsreverse,
                 vars_reversed=false)

end
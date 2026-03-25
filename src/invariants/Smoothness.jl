"""
    check_smoothness(f)

Using Oscar, check whether f defines a smooth hypersurface.

f is assumed to be homogeneous already. Otherwise Oscar
will throw an error.

note: the name issmooth / is_smooth is already taken by oscar,
and really we're just wrapping that method.
"""
function check_smoothness(f)
    p = characteristic(parent(f))
    nVars = length(gens(parent(f)))

    graded, _ = grade(parent(f))

    R, _ = quo(graded, ideal(graded, [graded(f)]))

    V = proj(R)

    is_smooth(V)
end
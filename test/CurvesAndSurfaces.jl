



"""
Works only for p = 7, 11, 13
"""
function test_ellipticcurve_1(p)
    n = 2
    d = 3

    F = GF(p)
    R, (x,y,z) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = y^2*z - x^3 - x*z^2 - z^3

    zeta = ZetaFunction.compute_all(f,precision,false,true)[2]

    t = gen(parent(zeta))

    correctzeta = if p == 7
        7t^2 - 3t + 1
    elseif p == 11
        11t^2 + 2t + 1
    elseif p == 13
        13t^2 + 4t + 1
    else
        @assert false "unsupported prime for testing"
    end

    @test zeta == correctzeta
end

function test_ellipticcurve_2(p)
    n = 2
    d = 3

    F = GF(p)
    R, (x,y,z) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = y^2*z - x^3 - x*z^2

    zeta = ZetaFunction.compute_all(f,precision,false,true)[2]

    t = gen(parent(zeta))

    correctzeta = if p == 7
        7t^2 + 1
    elseif p == 11
        11t^2 + 1
    elseif p == 13
        13t^2 - 6t + 1
    else
        @assert false "unsupported prime for testing"
    end

    @test zeta == correctzeta
end

function test_fermat_k3(p)

end

function test_fermatdeform_k3(p)

end

function test_highergenus_1(p)

end

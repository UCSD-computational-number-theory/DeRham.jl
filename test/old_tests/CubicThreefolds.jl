
function test_fermat_cubic(p)

    if p != 13
        error("testing not implemented for this prime")
    end

    R, (x,y,z,v,w) = polynomial_ring(GF(p),5)
    f = x^3 + y^3 + z^3 + v^3 + w^3

    zeta = DeRham.zeta_function(f)

    correctzeta = if p == 13
        [ZZ(51185893014090757),-7571877664806325,564530524121655,-27041473401150,913932150105,-22639590675,415990965,-5602350,53235,-325,1]
    else
        error("unreachable: are you trying to add support for a new prime?")
    end

    @test zeta == correctzeta
end


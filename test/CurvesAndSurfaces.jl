



"""
Works only for p = 7, 11, 13
"""
function test_ellipticcurve_1(p)
    n = 2
    d = 3

    F = GF(p)
    R, (x,y,z) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = y^2*z - x^3 - x*z^2 - z^3

    zeta = DeRham.compute_all(f,precision,false,true)[2]

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

    zeta = DeRham.compute_all(f,precision,false,true)[2]

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

"""
Data computed using Edgar Costa's `controlledreduction` library.
"""
function test_fermat_k3(p)
    n = 3
    d = 4

    abs_precision = 6
    series_precision = (4,4,3)

    F = GF(p)
    R, (x,y,z,w) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = x^4 + y^4 + z^4 + w^4 

    zeta = DeRham.compute_all(f,false,true)[2]

    t = gen(parent(zeta))
    ts = (t .^ (21:-1:0))

    convert_to_poly(row_of_coefs) = sum(ts .* transpose(row_of_coefs))
    
    correctzeta = if p == 7
        
        # controlledreduction takes about 30 seconds, mostly spent on  controlled reduction
        convert_to_poly([-558545864083284007 79792266297612001 113988951853731430 -16284135979104490 -10468373129424315 1495481875632045 569707381193160 -81386768741880 -20346692185470 2906670312210 498286339236 -71183762748 -8474257470 1210608210 98825160 -14117880 -756315 108045 3430 -490 -7 1])
    elseif p == 11

        #controlledreduction takes about 55 seconds
        convert_to_poly([-7400249944258160101211 672749994932560009201 611590904484145462910 -55599173134922314810 -22745116282468219695 2067737843860747245 501269780329878120 -45569980029988920 -7249769550225510 659069959111410 71898540993972 -6536230999452 -495169015110 45015365010 2338460520 -212587320 -7247295 658845 13310 -1210 -11 1])
      
    elseif p == 13

        #controlledreduction takes about 1 mintue
        convert_to_poly([-247064529073450392704413 -80405615970649536087235 -224910813903914786258 2508620616620588000570 188312900398839895003 -29841375627214911331 -3496390230500968632 172033060544399704 27724721295752390 -827978102045094 -111941095381388 8610853490876 376867593102 -74670735230 -2741627512 329708184 16651063 -621751 -49010 26 55 1])

    else
        @assert false "unsupported prime for testing"
    end

    @test zeta == correctzeta

end

"""
Data computed using Edgar Costa's `controlledreduction` library.
"""
function test_fermatdeform_k3(p)
    n = 3
    d = 4

    abs_precision = 6
    series_precision = (4,4,3)

    F = GF(p)
    R, (x,y,z,w) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = x^4 + y^4 + z^4 + w^4 + 2x*y*z*w

    zeta = DeRham.compute_all(f,false,true)[2]

    t = gen(parent(zeta))
    ts = (t .^ (21:-1:0))

    convert_to_poly(row_of_coefs) = sum(ts .* transpose(row_of_coefs))
    
    correctzeta = if p == 7
        
        convert_to_poly([-558545864083284007 -11398895185373143 100961643070447838 465261027974414 -8075602128413043 128184160768461 374379136212648 -16277353748376 -11045347186398 913524955266 213551288244 -30507326892 -2663338062 657187314 19765032 -9277464 -64827 83349 -98 -434 1 1])
    elseif p == 11

        convert_to_poly([-7400249944258160101211 -1284340899416705472111 433673550452394055518 89969571072874291238 -9511594081759437327 -2744452047306082707 63797972041984488 47227070212897608 1186325926400538 -491306696792142 -32681154997260 2971014090660 369125993082 -7366150638 -2423495448 -27056568 9619137 275517 -21538 -858 21 1])
      
    elseif p == 13

        convert_to_poly([-247064529073450392704413 247064529073450392704413 -108182101487783012190098 26314565226758029992186 -3538020111026967214597 147057070629482744861 34806407249581714760 -6961281449916342952 483994420334420294 10297753624136602 -4752809364986124 365600720383548 -4687188722866 -1303537692158 110939378056 -3282230120 -82055753 11681449 -514098 12506 -169 1])

    else
        @assert false "unsupported prime for testing"
    end

    @test zeta == correctzeta

end

"""
Data computed using Edgar Costa's `controlledreduction` library.
"""
function test_highergenus_1(p)
    n = 2
    d = 3
    series_precision = (4,4)
    absolute_precision = 7

    F = GF(p)
    R, (x,y,z) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = x^5 + y^5 + z^5 +x*z*y^3

    zeta = DeRham.compute_all(f,4,absolute_precision,false,true)[2]

    t = gen(parent(zeta))
    ts = (t .^ (21:-1:0))

    convert_to_poly(row_of_coefs) = sum(ts .* transpose(row_of_coefs))

    correctzeta = if p == 7
        # takes about 0.6 sec
        convert_to_poly([117649 -16807 -4802 -2401 2303 14 4 2 47 -7 -2 -1 1])
    elseif p == 11
        # takes about 1.1 sec
        convert_to_poly([1771561 2898918 2283996 1051490 290400 43230 5210 3930 2400 790 156 18 1])
    elseif p == 13
        # takes about 1.2 sec
        convert_to_poly([4826809 742586 628342 57122 68783 6188 5236 476 407 26 22 2 1])
    else
        @assert false "unsupported prime for testing"
    end

    @test zeta == correctzeta

end

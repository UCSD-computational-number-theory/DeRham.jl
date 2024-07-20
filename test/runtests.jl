
using Test
using Oscar

# source files, in the future replace this with `using DeRham`
include("../src/DeRham.jl")

include("FirstEllipticCurveExample.jl")
include("CurvesAndSurfaces.jl")
include("Orderings.jl")

@testset "The curve y^2 - x^3 - x - 1 = 0, reproducing Costa's results" begin
    testEllCurve1_7()
    testMonomialBasis()
    testLinAlgProb()
    testFrobTrans()
    testRedOfTerms()
    testT()
    #testFrobMat()
end


#@testset "Curves and Surfaces" begin # the "whiteboard picture"
#    test_ellipticcurve_1(7)
##    test_ellipticcurve_1(11)
##    test_ellipticcurve_1(13)
##
##    test_ellipticcurve_2(7)
##    test_ellipticcurve_2(11)
##    test_ellipticcurve_2(13)
##
##    test_fermat_k3(7)
##    test_fermat_k3(11)
##    test_fermat_k3(13)
##
##    test_fermatdeform_k3(7)
##    test_fermatdeform_k3(11)
##    test_fermatdeform_k3(13)
##    
##    test_highergenus_1(7)
##    test_highergenus_1(11)
##    test_highergenus_1(13)
#end

#@testset "Bigger primes" begin
#
#    test_ellipticcurve_1(next_prime(20))
#    test_ellipticcurve_1(next_prime(50))
#    test_ellipticcurve_1(next_prime(100))
#    test_ellipticcurve_1(next_prime(1000))
#    test_ellipticcurve_1(next_prime(10000))
#    test_ellipticcurve_1(next_prime(100000))
#
#    # maybe do it with another variety?
#
#end
#

@testset "Monomial orderings" begin
    test_supported_monomial_orderings()

    test_reversing_variables()

end

#@testset "Elliptic curves" begin
#
#    test_lmfdb_elliptic_curves()
#
#end



using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using Primes
#inclute("../src/DeRham.jl")

include("FirstEllipticCurveExample.jl")
#include("CurvesAndSurfaces.jl")
include("Orderings.jl")
include("Precision.jl")
include("NaivePointCounts.jl")
include("ManageCSVTests.jl")

@testset "The curve y^2 - x^3 - x - 1 = 0, reproducing Costa's results" begin
    #testEllCurve1_7() TODO: why is this failing?
    testMonomialBasis()
    testLinAlgProb()
    testFrobTrans()
    testRedOfTerms()
    testT()
    #testFrobMat()
end


@testset "Curves and Surfaces" begin # the "whiteboard picture"

    runcsvtest("ellipticcurves.csv")
    runcsvtest("highergenus.csv")
    
    zf = f -> DeRham.zeta_function(f,algorithm=:naive,fastevaluation=true)
    runcsvtest("k3surfaces.csv",zeta_function=zf)
    #runcsvtest("othersurfaces.csv")
end

#@testset "Threefolds" begin
#    runcsvtest("threefolds.csv")
#end

#TODO: we don't have this data rn
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
    #test_naive_algorithm()
    #TODO we need to re-enable this when fastevaluation gets fixed for costachunks
    #test_fastevaluation()

    #test_reversing_variables()

end

@testset "Precision" begin
    test_hodge_polygon_values()
    test_hodge_polygon_examples()
    test_algorithm_precision()
    #test_series_precision()
end

#TODO: we need to re-enable this
#@testset "Naive Point Counts" begin
#    test_fermat_cubic_naive()
#end

#@testset "Elliptic curves" begin
#
#    test_lmfdb_elliptic_curves()
#
#end


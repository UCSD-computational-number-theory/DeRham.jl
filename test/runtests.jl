
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

function larger_tests()
    zf = f -> DeRham.zeta_function(f,algorithm=:naive,fastevaluation=true)

    # Keep track of which files to run with (dim,degree) tuples
    fermatfiles = [(1,3), (1,4), (1,5), (1,6), (1,7), (1,8),
                   (2,4), (2,5)]
    randomfiles = [(1,3), (1,4), (1,5), (1,6), (2,4)]
    manyfiles = [(1,3),(1,4),(1,5)]
    
    for (dim,deg) in fermatfiles
        filename = "dim_$(dim)_deg_$(deg)_fermat.csv"
        runcsvtest(filename,zeta_function=zf)
    end

    for (dim,deg) in randomfiles
        filename = "dim_$(dim)_deg_$(deg)_random.csv"
        runcsvtest(filename,zeta_function=zf)
    end

    for (dim,deg) in manyfiles
        filename = "dim_$(dim)_deg_$(deg)_many.csv"
        runcsvtest(filename,zeta_function=zf)
    end
end

@testset "CPU Fast Evaluation + Naive Strategy" begin 

    zf = f -> DeRham.zeta_function(f,algorithm=:naive,fastevaluation=true)

    # Keep track of which files to run with (dim,degree) tuples
    fermatfiles = [(1,3), (1,4), (1,5), (1,6), (1,7), (1,8),
                   (2,4)] # (2,5)
    randomfiles = [(1,3), (1,4), (1,5), (1,6), (2,4)]
    manyfiles = [(1,3),(1,4),(1,5)]
    
    for (dim,deg) in fermatfiles
        filename = "dim_$(dim)_deg_$(deg)_fermat.csv"
        runcsvtest(filename,zeta_function=zf)
    end

    for (dim,deg) in randomfiles
        filename = "dim_$(dim)_deg_$(deg)_random.csv"
        runcsvtest(filename,zeta_function=zf)
    end

    for (dim,deg) in manyfiles
        filename = "dim_$(dim)_deg_$(deg)_many.csv"
        runcsvtest(filename,zeta_function=zf)
    end

    #runcsvtest("ellipticcurves.csv")
    #runcsvtest("highergenus.csv")
    
    #runcsvtest("k3surfaces.csv",zeta_function=zf)
    #runcsvtest("othersurfaces.csv")
end

@testset "CPU S=[0,1,2], Fast Evaluation + Naive Strategy" begin

    zf = f -> DeRham.zeta_function(f,S=[0,1,2],algorithm=:naive,fastevaluation=true)

    fermatfiles = [(2,3),(2,4),(3,3)]
    randomfiles = [(2,3),(2,4)]

    for (dim,deg) in fermatfiles
        filename = "dim_$(dim)_deg_$(deg)_fermat.csv"
        runcsvtest(filename,zeta_function=zf)
    end

    for (dim,deg) in randomfiles
        filename = "dim_$(dim)_deg_$(deg)_random.csv"
        runcsvtest(filename,zeta_function=zf)
    end
end

@testset "CPU S=[0,1], Fast Evaluation + Naive Strategy" begin

    # do dimension 1, degree 3 and 4
    zf = f -> DeRham.zeta_function(f,S=[0,1],algorithm=:naive,fastevaluation=true)

    for i in 3:4
        runcsvtest("dim_1_deg_$(i)_fermat.csv",zeta_function=zf)
        runcsvtest("dim_1_deg_$(i)_random.csv",zeta_function=zf)
        runcsvtest("dim_1_deg_$(i)_many.csv",zeta_function=zf)
    end
end

@testset "CPU Bigints + Fast Evaluatoin + Naive Strategy" begin
    # do dimension 1, degree 3 and 4
    zf = f -> DeRham.zeta_function(f,algorithm=:naive,fastevaluation=true,always_use_bigints=true)

    for i in 3:4
        runcsvtest("dim_1_deg_$(i)_fermat.csv",zeta_function=zf)
        runcsvtest("dim_1_deg_$(i)_random.csv",zeta_function=zf)
        runcsvtest("dim_1_deg_$(i)_many.csv",zeta_function=zf)
    end

    #TODO: insert specific examples that need the use of bigints here
end

@testset "CPU Costachunks Strategy" begin

    # do dimension 1, degree 3 and 4, and K3 surfaces
    zf = f -> DeRham.zeta_function(f,algorithm=:costachunks)

    for i in 3:4
        runcsvtest("dim_1_deg_$(i)_fermat.csv",zeta_function=zf)
        runcsvtest("dim_1_deg_$(i)_random.csv",zeta_function=zf)
        runcsvtest("dim_1_deg_$(i)_many.csv",zeta_function=zf)
    end

    runcsvtest("dim_1_deg_4_fermat.csv",zeta_function=zf)
    runcsvtest("dim_1_deg_4_random.csv",zeta_function=zf)
    runcsvtest("dim_1_deg_4_many.csv",zeta_function=zf)
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


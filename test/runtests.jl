
using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using Primes
using CUDA
#inclute("../src/DeRham.jl")

include("FirstEllipticCurveExample.jl")
#include("CurvesAndSurfaces.jl")
include("Orderings.jl")
include("Precision.jl")
include("NaivePointCounts.jl")
include("ManageCSVTests.jl")

@testset "The curve y^2 - x^3 - x - 1 = 0, reproducing Costa's results" begin
    testEllCurve1_7() #TODO: why is this failing?
    testMonomialBasis()
    testLinAlgProb()
    testFrobTrans()
    #testRedOfTerms()
    testT()
    #testFrobMat()
end

# currently, this runs on all the examples that can be done with a full S
function larger_tests(zf)

    # Keep track of which files to run with (dim,degree) tuples
    fermatfiles = [(1,3), (1,4), (1,5), (1,6), (1,7), (1,8),
                   (2,4), (2,5)]
    randomfiles = [(1,3), (1,4), (1,5), (1,6), (1,7), (2,4)]
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

    # No BigInts necessary
    teststring("3;47;x1^6 + x2^6 + x3^6;ZZRingElem[52599132235830049,0,11191304731027670,0,1071507899779245,0,60794774455560,0,2263635219090,0,57794941764,0,1024733010,0,12458760,0,99405,0,470,0,1]",zf)

    # Needs UInt64, not Int64
    teststring("3;53;x1^6 + x2^6 + x3^6;ZZRingElem[174887470365513049,0,32997635918021330,0,2801686068511245,0,140965336780440,0,4654515837090,0,105385264236,0,1657001010,0,17865240,0,126405,0,530,0,1
]",zf)

    # Needs BigInts
    teststring("3;59;x1^6 + x2^6 + x3^6;ZZRingElem[511116753300641401,0,86629958186549390,0,6607369692194445,0,298638178178280,0,8857912064610,0,180160923348,0,2544645810,0,24645480,0,156645,0,590,0,1]",zf)
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

@testset "CPU Cubic Surfaces" begin
    # do dimension 2, degree 3 and 
    zf = f -> DeRham.zeta_function(f,S=[1,2],algorithm=:naive,fastevaluation=true)

    runcsvtest("dim_2_deg_3_many.csv", zeta_function=zf)
end 

if CUDA.functional()
    @testset "GPU (CUDA) S=[0,1,2], Fast Evaluation + Naive Strategy" begin
    
        zf = f -> DeRham.zeta_function(f,S=[0,1,2],algorithm=:naive,fastevaluation=true,use_gpu=true)
    
        fermatfiles = [(2,3)]#,(2,4),(3,3)]
        randomfiles = [(2,3)]#,(2,4)]
    
        for (dim,deg) in fermatfiles
            filename = "dim_$(dim)_deg_$(deg)_fermat.csv"
            runcsvtest(filename,zeta_function=zf)
        end
    
        #for (dim,deg) in randomfiles
        #    filename = "dim_$(dim)_deg_$(deg)_random.csv"
        #    runcsvtest(filename,zeta_function=zf)
        #end
    end
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


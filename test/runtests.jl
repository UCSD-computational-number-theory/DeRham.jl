
using Test
using Oscar

# source files, in the future replace this with `using DeRham`
using DeRham
using CUDA

include("FirstEllipticCurveExample.jl")
include("Precision.jl")

@testset "DeRham.jl" begin

    @testset "The curve y^2 - x^3 - x - 1 = 0, reproducing Costa's results" begin
        testEllCurve1_7() #TODO: why is this failing?
        testMonomialBasis()
        testLinAlgProb()
        testFrobTrans()
        testT()
    end

    @testset "Hodge polygon" begin
        test_hodge_polygon_values()
        test_hodge_polygon_examples()
    end

end

# Declarative E2E catalogue (test/e2e/) — see test/e2e/README.md. Bare
# Pkg.test() runs conventional tests above plus the ci E2E workflow;
# test_args (--workflow=/--name=/--tag=) select only within the E2E suite.
include("e2e/E2E.jl")
using .E2E

include("e2e/test/runtests.jl")

@testset "E2E" begin
    E2E.run!(ARGS)
end

include("quality_gates.jl")

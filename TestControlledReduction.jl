module TestControlledReduction

using Tests
using Oscar

include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("AutomatedScript.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")

function runTests()
    testEllCurve1_7()
end

function testEllCurve1_7()
    n = 2
    d = 3
    p = 7
    Fp = GF(p,1)

    R = Fp
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    @test ControlledReduction.computeFrobeniusMatrix(n,d,f,precision,7,R,PR,vars) == [84 48; 294 17]
end
end
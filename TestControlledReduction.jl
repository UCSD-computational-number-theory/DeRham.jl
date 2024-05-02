module TestControlledReduction

using Test
using Oscar

include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("FindMonomialBasis.jl")
include("AutomatedScript.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")

function runTests()
    @testset "All tests" begin
        testEllCurve1_7()
        testMonomialBasis()
        testLinAlgProb()
        testFrobTrans()
        testRedOfTerms()
        testT()
        testFrobMat()
    end
end

function testEllCurve1_7()
    n = 2
    d = 3
    p = 7
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    @test ControlledReduction.computeAll(n,d,f,precision,7,R,PR,vars) == [84 48; 294 17]
end

function testMonomialBasis()
    n = 2
    d = 3
    p = 7
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    @test CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR) == 1
end

function testLinAlgProb()
    n = 2
    d = 3
    p = 7
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    S = [0,1,2]
    @test CopiedFindMonomialBasis.psuedo_inverse_controlled(f,S,R,PR) == 1
end

function testFrobTrans()
    n = 2
    d = 3
    p = 7
    N = 6
    M = 15
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    PrecisionRing, = residue_ring(ZZ,p^M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR)
    fLift = ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    M = 15
    @test ControlledReduction.applyFrobeniusToBasis(Basis,n,d,fLift,N,p,PrecisionRing,PrecisionRingPoly) == 1
end

function testRedOfTerms()
    n = 2
    d = 3
    p = 7
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    N = 6
    M = 15
    PrecisionRing, = residue_ring(ZZ,p^M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR)
    fLift = ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    S = [0,1,2]
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = ControlledReduction.applyFrobeniusToBasis(Basis,n,d,fLift,N,p,PrecisionRing,PrecisionRingPoly)
    psuedoInverseMatTemp = CopiedFindMonomialBasis.psuedo_inverse_controlled(f,S,R,PR)
    psuedoInverseMat = zeros(PrecisionRing,nrows(psuedoInverseMatTemp),ncols(psuedoInverseMatTemp))
    for i in 1:nrows(psuedoInverseMat)
        for j in 1:ncols(psuedoInverseMat)
            psuedoInverseMat[i,j] = PrecisionRing(lift(ZZ,psuedoInverseMatTemp[i,j]))
        end
    end
    @test ControlledReduction.computeReductionOfTransformLA(FBasis,n,d,p,N,S,fLift,psuedoInverseMat,PrecisionRing,PrecisionRingPoly) == 1
end

function testT()
    n = 2
    d = 3
    p = 7
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    M = 15
    PrecisionRing, = residue_ring(ZZ,p^M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR)
    fLift = ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    @test ControlledReduction.computeT(BasisTLift,fLift,n,d,PrecisionRing,PrecisionRingPoly) == 1
end

function testFrobMat()
    n = 2
    d = 3
    p = 7
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    N = 6
    M = 15
    PrecisionRing, = residue_ring(ZZ,p^M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR)
    fLift = ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,ControlledReduction.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    T = ControlledReduction.computeT(BasisTLift,fLift,n,d,PrecisionRing,PrecisionRingPoly)
    S = [0,1,2]
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = ControlledReduction.applyFrobeniusToBasis(Basis,n,d,fLift,N,p,PrecisionRing,PrecisionRingPoly)
    psuedoInverseMatTemp = CopiedFindMonomialBasis.psuedo_inverse_controlled(f,S,R,PR)
    psuedoInverseMat = zeros(PrecisionRing,nrows(psuedoInverseMatTemp),ncols(psuedoInverseMatTemp))
    for i in 1:nrows(psuedoInverseMat)
        for j in 1:ncols(psuedoInverseMat)
            psuedoInverseMat[i,j] = PrecisionRing(lift(ZZ,psuedoInverseMatTemp[i,j]))
        end
    end
    Reductions = ControlledReduction.computeReductionOfTransformLA(FBasis,n,d,p,N,S,fLift,psuedoInverseMat,PrecisionRing,PrecisionRingPoly)
    @test ControlledReduction.computeFrobeniusMatrix(n,d,Reductions,T) == 1
end

end
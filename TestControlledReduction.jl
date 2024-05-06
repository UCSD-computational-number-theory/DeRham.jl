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
include("ZetaFunction.jl")

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
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    @test ZetaFunction.computeAll(n,d,f,precision,7,R,PR,vars) == [84 48; 294 17]
end

function testMonomialBasis()
    n = 2
    d = 3
    p = 7
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    @test CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR) == [[1],[z^3]]
end

function testLinAlgProb()
    n = 2
    d = 3
    p = 7
    R = GF(p)
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
    N = 2 # the series precision
    M = 3 # tge absolute precision
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    PrecisionRing = residue_ring(ZZ,p^M)
    println(eltype(PrecisionRing))
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

    frobterms = ControlledReduction.applyFrobeniusToBasis(Basis,n,d,fLift,N,p,PrecisionRing,PrecisionRingPoly)

    x0,x1,x2 = PVars

    #TODO:test failing
    @test frobterms[1][1] == [133*x0^6*x1^6*x2^6, 7]

    #TODO: test failing
    @test frobterms[1][2] == [1*x0^27*x1^6*x2^6 + 
                              1*x0^13*x1^6*x2^20 + 
                              342*x0^6*x1^20*x2^13 + 
                              1*x0^6*x1^6*x2^27, 14]
    
    #TODO: test failing
    @test frobterms[2][1] == [56*x0^6*x1^6*x2^27, 14]
    
    @test frobterms[2][2][2] == 21
    bigpolyterms = terms(frobterms[2][2][1])
   
    coefficients = leading_coefficient.(bigpolyterms)
    exp_vecs = leading_exponent_vector.(bigpolyterms)

    # Costa's code shows:
    #
    # [4 1 4] --> 2
    # [4 3 2] --> 341
    # [5 1 3] --> 2
    # [7 1 1] --> 2
    #
    # To get the monomial from 
    # 
    # key --> value
    # 
    # I think you need to do
    #
    # prod([x,y,z] .^ (p .* key)) * value
    
    #TODO:test failing
    @test coefficients == [2,2,341,2]
    #TODO:test failing
    @test exp_vecs == [[27 6 27],
                       [27 20 13],
                       [35 6 20],
                       [48 6 6]]

end

function testRedOfTerms()
    n = 2
    d = 3
    p = 7
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    N = 6
    M = 15
    PrecisionRing = residue_ring(ZZ,p^M)
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
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    M = 15
    PrecisionRing = residue_ring(ZZ,p^M)
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
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    N = 6
    M = 15
    PrecisionRing = residue_ring(ZZ,p^M)
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
    @test ZetaFunction.computeFrobeniusMatrix(n,d,Reductions,T) == 1
end

end

module TestControlledReduction

using Test
using Oscar

include("ControlledReduction.jl")
include("PrecisionEstimate.jl")
include("CopiedFindMonomialBasis.jl")
include("Utils.jl")
include("SmallestSubsetSmooth.jl")
include("ZetaFunction.jl")
include("Frobenius.jl")
include("FinalReduction.jl")

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
    l = d * n - n + d - length(S)
    M = 3
    @test Array(CopiedFindMonomialBasis.pseudo_inverse_controlled_lifted(f,S,l,R,PR,M)) == 
      [155 0 0 0 0 11 0 0 0 221 0 0 310 0 22; 
       0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 
       0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 
       225 0 0 0 0 155 0 0 0 11 0 0 221 0 310; 
       0 0 0 0 0 0 0 0 0 0 0 0 0 114 0; 
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 114; 
       0 0 0 1 0 0 172 0 0 0 0 0 0 0 0; 
       59 0 0 0 1 94 0 172 0 166 0 0 61 0 188; 
       0 0 0 0 0 0 0 0 172 0 0 0 0 286 0; 
       0 57 0 173 0 0 0 0 0 0 172 0 0 114 0; 
       83 0 57 0 173 59 0 0 0 94 0 172 166 0 61; 
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
       114 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
       0 114 0 0 0 0 0 0 0 0 0 0 0 0 0; 
       166 0 114 0 0 118 0 0 0 188 0 0 332 0 236; 
       11 0 0 0 0 221 0 0 0 310 0 0 22 0 214; 
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] 
end

function testFrobTrans()
    n = 2
    d = 3
    p = 7
    N = 2 # the series precision
    M = 3 # the absolute precision
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    #x,y,z = Vars
    #f = y^2*z - x^3 - x*z^2 - z^3
    x0,x1,x2 = Vars
    f = x1^2*x2 - x0^3 - x0*x2^2 - x2^3
    PrecisionRing, = residue_ring(ZZ, p^M)
    println(eltype(PrecisionRing))
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f, R, PR)
    fLift = Utils.liftCoefficients(PrecisionRing, PrecisionRingPoly, f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,Utils.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    #M = 15

    frobterms = Frobenius.applyFrobeniusToBasis(Basis,n,d,fLift,N,p,PrecisionRing,PrecisionRingPoly)

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
    @test coefficients == [2,341,2,2]
    #TODO:test failing
    @test exp_vecs == [[27, 6, 27],
                       [6, 20, 34],
                       [13, 6, 41],
                       [6, 6, 48]]

end

function testRedOfTerms()
    n = 2
    d = 3
    p = 7
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    N = 2
    M = 3
    
    # TODO: change other instances of pseudo_inverse_controlled in this file and ZetaFunction.jl to this method
    S = [0,1,2]
    pseudoInverseMat = Array(CopiedFindMonomialBasis.pseudo_inverse_controlled_lifted(f,S,R,PR,M))

    PrecisionRing = parent(pseudoInverseMat[1])
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR)
    fLift = Utils.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,Utils.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = Frobenius.applyFrobeniusToBasis(Basis,n,d,fLift,N,p,PrecisionRing,PrecisionRingPoly)
    #pseudoInverseMat = zeros(PrecisionRing,nrows(pseudoInverseMatTemp),ncols(pseudoInverseMatTemp))
    #for i in 1:nrows(pseudoInverseMat)
    #    for j in 1:ncols(pseudoInverseMat)
    #        pseudoInverseMat[i,j] = PrecisionRing(lift(ZZ,pseudoInverseMatTemp[i,j]))
    #    end
    #end
    @test ControlledReduction.reducetransform_LA_descending(FBasis,n,d,p,N,S,fLift,pseudoInverseMat,PrecisionRing,PrecisionRingPoly) == 1
end

function testT()
    n = 2
    d = 3
    p = 7
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    M = 3
    PrecisionRing, = residue_ring(ZZ,p^M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = CopiedFindMonomialBasis.compute_monomial_bases(f,R,PR)
    fLift = Utils.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,Utils.liftCoefficients(ZZ,PRZZ,j,false))
        end
        push!(BasisTLift,temp)
    end
    @test FinalReduction.computeT(BasisTLift,fLift,n,d,PrecisionRing,PrecisionRingPoly) == 1
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
    fLift = Utils.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp, Utils.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    T = FinalReduction.computeT(BasisTLift,fLift,n,d,PrecisionRing,PrecisionRingPoly)
    S = [0,1,2]
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = Frobenius.applyFrobeniusToBasis(Basis,n,d,fLift,N,p,PrecisionRing,PrecisionRingPoly)
    pseudoInverseMatTemp = CopiedFindMonomialBasis.pseudo_inverse_controlled(f,S,R,PR)
    pseudoInverseMat = zeros(PrecisionRing,nrows(pseudoInverseMatTemp),ncols(pseudoInverseMatTemp))
    for i in 1:nrows(pseudoInverseMat)
        for j in 1:ncols(pseudoInverseMat)
            pseudoInverseMat[i,j] = PrecisionRing(lift(ZZ,pseudoInverseMatTemp[i,j]))
        end
    end
    Reductions = ControlledReduction.reducetransform_LA(FBasis,n,d,p,N,S,fLift,pseudoInverseMat,PrecisionRing,PrecisionRingPoly)
    @test ZetaFunction.computeFrobeniusMatrix(n,d,Reductions,T) == 1
end

end

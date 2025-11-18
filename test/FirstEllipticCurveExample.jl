#module TestControlledReduction


# vars_reversed =true, do not use fast evaluation, bigints, or the gpu
first_ellcurve_params() = DeRham.ZetaFunctionParams(0,true,:costachunks,:invlex,false,false,false,false)
function first_ellcurve_cache()
    n=2
    d = 3 #total_degree(f)
    S = collect(0:n)
    DeRham.controlled_reduction_cache(n,d,S,first_ellcurve_params())
end

function testEllCurve1_7()
    n = 2
    #d = 3
    p = 7
    #series_precision = 2
    #absolute_precision = 3
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    frobmat = DeRham.zeta_function(f,givefrobmat=true)[1]
    R = parent(frobmat[1,1])
    @test frobmat == R[231 11; 294 17]
end

function testMonomialBasis()
    n = 2
    #d = 3
    p = 7
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    params = first_ellcurve_params()
    cache = first_ellcurve_cache()
    @test DeRham.compute_monomial_bases(f,params,cache) == [[1],[z^3]]
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
    params = first_ellcurve_params()
    cache = first_ellcurve_cache()
    psinv = DeRham.pseudo_inverse_controlled_lifted(f,S,l,M,params,cache)
    println(psinv)
    @test Array(psinv) == 
        [114 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
         0 114 0 0 0 0 0 0 0 0 0 0 0 0 0; 
         166 0 114 0 0 118 0 0 0 188 0 0 332 0 122; 
         11 0 0 0 0 221 0 0 0 310 0 0 22 0 99; 
         0 0 0 0 0 0 0 0 0 0 0 0 0 342 0; 
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 342; 
         0 0 0 1 0 0 172 0 0 0 0 0 0 170 0; 
         59 0 0 0 1 94 0 172 0 166 0 0 61 0 188; 
         0 0 0 0 0 0 0 0 172 0 0 0 0 0 0; 
         0 57 0 173 0 0 0 0 0 0 172 0 0 0 0; 
         83 0 57 0 173 59 0 0 0 94 0 172 166 0 61; 
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
         155 0 0 0 0 11 0 0 0 221 0 0 310 0 23; 
         0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 
         0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 
         225 0 0 0 0 155 0 0 0 11 0 0 221 0 310; 
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]

      # This is the matrix when the variables are reversed:
      #
      # [155 0 0 0 0 11 0 0 0 221 0 0 310 0 22; 
      #  0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 
      #  0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 
      #  225 0 0 0 0 155 0 0 0 11 0 0 221 0 310; 
      #  0 0 0 0 0 0 0 0 0 0 0 0 0 114 0; 
      #  0 0 0 0 0 0 0 0 0 0 0 0 0 0 114; 
      #  0 0 0 1 0 0 172 0 0 0 0 0 0 0 0; 
      #  59 0 0 0 1 94 0 172 0 166 0 0 61 0 188; 
      #  0 0 0 0 0 0 0 0 172 0 0 0 0 286 0; 
      #  0 57 0 173 0 0 0 0 0 0 172 0 0 114 0; 
      #  83 0 57 0 173 59 0 0 0 94 0 172 166 0 61; 
      #  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      #  114 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      #  0 114 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      #  166 0 114 0 0 118 0 0 0 188 0 0 332 0 236; 
      #  11 0 0 0 0 221 0 0 0 310 0 0 22 0 214; 
      #  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      #  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] 
end

function testFrobTrans()
    n = 2
    d = 3
    p = 7
    N = [2,2] # the series precision
    M = 3 # the absolute precision
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    #x,y,z = Vars
    #f = y^2*z - x^3 - x*z^2 - z^3
    x0,x1,x2 = Vars
    f = x1^2*x2 - x0^3 - x0*x2^2 - x2^3
    PrecisionRing, = residue_ring(ZZ, p^M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])

    params = first_ellcurve_params()
    cache = first_ellcurve_cache()

    BasisT = DeRham.compute_monomial_bases(f, params, cache)
    fLift = DeRham.liftCoefficients(PrecisionRing, PrecisionRingPoly, f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,DeRham.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
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

    frobterms = DeRham.applyFrobeniusToBasis(Basis,fLift,N,p,params)

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

#=
function testRedOfTerms()
    n = 2
    d = 3
    p = 7
    R = GF(p,1)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x0,x1,x2 = Vars
    f = x1^2*x2 - x0^3 - x0*x2^2 - x2^3
    S = [0,1,2]
    N = [2,2]
    M = 3
    params = first_ellcurve_params()
    
    # TODO: change other instances of pseudo_inverse_controlled in this file and ZetaFunction.jl to this method
    S = [0,1,2]
    #pseudoInverseMat = Array(CopiedFindMonomialBasis.pseudo_inverse_controlled_lifted(f,S,R,PR,M))
    pseudo_inverse_mat = 
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
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;]

    PrecisionRing, = residue_ring(ZZ, p^M)

    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = DeRham.compute_monomial_bases(f,R,PR,:invlex)
    fLift = DeRham.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp,DeRham.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = DeRham.applyFrobeniusToBasis(Basis, fLift, N, p, params)
    #pseudoInverseMat = zeros(PrecisionRing,nrows(pseudoInverseMatTemp),ncols(pseudoInverseMatTemp))
    #for i in 1:nrows(pseudoInverseMat)
    #    for j in 1:ncols(pseudoInverseMat) 
    #        pseudoInverseMat[i,j] = PrecisionRing(lift(ZZ,pseudoInverseMatTemp[i,j]))
    #    end
    #end
    cache = DeRham.controlled_reduction_cache(n,d,[0,1,2],params.termorder)

    Reductions = DeRham.reducetransform(FBasis, N, S, fLift, matrix(PrecisionRing,pseudo_inverse_mat), p, params,cache)
    ev = DeRham.gen_exp_vec(n+1,n*d-n-1,:invlex)
    reductions_as_rows = DeRham.convert_p_to_m([Reductions[1][1][1],Reductions[2][1][1]],ev)
    RR = parent(reductions_as_rows[1,1])
    @test reductions_as_rows == RR[86 0 98 0 226 0 329 236 0 272; 
                                   133 0 224 0 203 0 238 91 0 322]
end
=#

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
    params = first_ellcurve_params()
    cache = first_ellcurve_cache()
    precisionring, pi = residue_ring(ZZ,p^M)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i in 0:n])
    basis = DeRham.compute_monomial_bases(f,params,cache)
    @test Array(DeRham.computeT(f,basis,M,params,cache)) == 
    [257 0 85 0 0 0 172 257 0 0;
     172 0 52 0 114 0 0 170 0 1]
end

#=
function testFrobMat()
    n = 2
    d = 3
    p = 7
    R = GF(p)
    PR, Vars = polynomial_ring(R, ["x$i" for i in 0:n])
    x,y,z = Vars
    f = y^2*z - x^3 - x*z^2 - z^3
    N = [2,2]
    M = 3
    PrecisionRing, = residue_ring(ZZ,p^M)
    PrecisionRingPoly, PVars = polynomial_ring(PrecisionRing, ["x$i" for i in 0:n])
    BasisT = DeRham.compute_monomial_bases(f,R,PR,:invlex)
    fLift = DeRham.liftCoefficients(PrecisionRing,PrecisionRingPoly,f)
    BasisTLift = []
    for i in BasisT
        temp = []
        for j in i
            push!(temp, DeRham.liftCoefficients(PrecisionRing,PrecisionRingPoly,j))
        end
        push!(BasisTLift,temp)
    end
    T = DeRham.computeT(f,BasisT,M)
    Basis = []
    for i in 1:n
        for j in BasisTLift[i]
            push!(Basis,[j,i])
        end
    end
    FBasis = DeRham.applyFrobeniusToBasis(Basis,fLift,N,p)

    l = d * n - n + d - length(S)
    pseudo_inverse_mat_new = DeRham.pseudo_inverse_controlled_lifted(f,S,l,M)
    MS = matrix_space(ZZ, nrows(pseudo_inverse_mat_new), ncols(pseudo_inverse_mat_new))
    pseudo_inverse_mat = MS()
    for i in 1:nrows(pseudo_inverse_mat_new)
        for j in 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i,j] = ZZ(pseudo_inverse_mat_new[i,j])
        end
    end 

    Reductions = DeRham.reducetransform_LA_descending(FBasis,N,S,fLift,pseudo_inverse_mat,p)
    frobMat = DeRham.compute_frobenius_matrix(n,p,d,N,Reductions,T,Basis)
    RR = parent(frobMat[1,1])
    @test frobMat == RR[231 11; 294 17]
end
=#
#end

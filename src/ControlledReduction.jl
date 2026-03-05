
"""
    ControlledReductionContext{MatrixType,VectorType}

A controlled reduction context is used for in-place evaluation
chunks of controlled reduction.
It contains all of the major allocations necessary, so you can
create one context at the beginning (or one per thread)
and then not have to worry about allocating anything else.

MatrixType - the type of A and B, for example zzModMatrix
VectorType - the type of g, for exampel Vector{UInt64}

fields
------
Ruvs - an AbstractPEP which containes the R_uv necessary for the reduction.
  This is intended to be shared across all contexts (i.e. all threads)
  for a single example.
A - the matrix which will store A in Ay + B
B - the matrix which will store B in Ay + B
temp - temporary matrix used for intermediate calculations
g - the vector g which contains the vector which is being reduced
g_temp - temporary vector used for intermediate calculations

Current examples that are in use:

ControlledReductionContext{zzModMatrix,Vector{UInt64}}
ControlledReductionContext{ZZModMatrix,Vector{ZZRingElem}}
ControlledReductionContext{CuModMatrix,CuModVector}

"""
struct ControlledReductionContext{MatrixType,VectorType}
    Ruvs::AbstractPEP{MatrixType}
    A::MatrixType
    B::MatrixType
    temp::MatrixType
    g::VectorType
    g_temp::VectorType
end

"""
    chooseV(I, d, S)
Choose the direction of reduction V following Edgar's method

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f
* "S" -- list 
"""
function chooseV(I, d, S)
    v = zeros(Int,length(I))
    n = length(I)

    sum = 0
    
    for i in S 
        if I[n-i] > 0
            v[n-i] = v[n-i] + 1
            sum = sum + 1
            if sum >= d 
                break
            end
        end 
    end 
    
    while sum < d
        for i in 0:(n-1) 
            if I[n-i] > v[n-i]
                v[n-i] = v[n-i] + 1
                sum = sum + 1
                break
            end 
        end 
    end 
    return v
end 

function varbyvar_chooseV(I, d)
    v = zeros(Int,length(I))
    n = length(I)
    J = copy(I)
    i = 0
    while i < d
        for j in eachindex(J)
            if J[n-j+1] > 0
                J[n-j+1] = J[n-j+1] - 1
                v[n-j+1] = v[n-j+1] + 1
                i = i + 1
                break
            elseif j == n
                throw(error)
            end
        end
    end
    return v
end

"""
    skim_chooseV(I, d)

choose direction of reduction V by scanning from right to left
Note: this differs from Edgar's method, see chooseV

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f 
"""
function skim_chooseV(I, d)
    V = zeros(Int,length(I))
    i = 0
    s = length(I)
    foundNonZero = false
    while i < d
        #FIXME reversed to match Costa's
        if s == 0 && foundNonZero == false
            return V
        elseif s == 0
            s = length(I)
            foundNonZero = false
        end
        if (I - V)[s] > 0
            V[s] = V[s] + 1
            i = i + 1
            foundNonZero = true
        end
        s = s-1
    end
    return V
end

"""
    rev_chooseV(I, d)

choose direction of reduction in the same way as Costa's code

INPUTS: 
* "I" -- list/tuple, exponents of monomials
* "d" -- integer, degree of f 
"""
function rev_chooseV(I, d, S)
    reverse!(I)

    V = chooseV(I,d,S)

    reverse!(I)

    reverse!(V)

    return V
end

"""
    tweak(I,m)

If for a vectors of ints, I, we let |I| = sum(I[i]). This function returns I after removing an integer vector, J, with |J| = m
from it

Removes from the "front" of the vector"

INPUTS
* "I" -- vector of nonnegative integers 
* "m" -- nonnegative integer
"""
function tweak(J,m)
    count = 0
    o = m
    I = copy(J)
    while m > 0
        for i in axes(I,1)
            if (I[i] > 0)
                I[i] = I[i] - 1
                m = m - 1
                break
            end
        end
        if count > length(I)*o
            return I
        end
        count = count+1
    end
    return I
end


"""
    rev_tweak(J,m)

If for a vectors of ints, I, we let |I| = sum(I[i]). This function returns I after removing an integer vector, J, with |J| = m
from it

unlike tweak, this, removes from the "back" of the vector instead of the front.

INPUTS
* "I" -- vector of nonnegative integers 
* "m" -- nonnegative integer
"""
function rev_tweak(J,m)
    count = 0
    o = m
    I = copy(J)
    while m > 0
        for i in reverse(axes(I,1))
            if (I[i] > 0)
                I[i] = I[i] - 1
                m = m - 1
                break
            end
        end
        if count > length(I)*o
            return I
        end
        count = count+1
    end
    return I
end

"""
    undo_rev_tweak(I,p)

gives the unique vector J containing only
multiples of p such that tweak(J, |J-I|) = I

We're assuming that |J-I| < p here.

"""
function undo_rev_tweak(I,p)

    # is this correct for large d and n?
    J = copy(I)

    k = 1
    #accumulator = 0
    while k ≤ length(J)
        while (J[k] % p) != 0
            J[k] = J[k] + 1
            #accumulator += 1
        end
        k = k + 1
    end

    # @assert accumulator == m

    J
end

"""
    costadata_of_intial_term
Returns the data used by costa's code given a polynomial term.
More specifically, given a polynomial h, it writes h in the form 
h = x^{rev_tweak(U)}g

Note: this only works with a single term, so it should
only be used at the beginning of reduction
"""
function costadata_of_initial_term!(term,g,n,d,p,S,cache,params)
    termorder = params.termorder
    vars_reversed = params.vars_reversed

    R = base_ring(parent(term[1])) 
    i = term

    Stilda = zeros(Int, n+1)
    for j in S
        Stilda[j+1] = 1
    end

    U = exponent_vector(i[1],1) + Stilda
    
    gCoeff = coeff(i[1],1)

    g_exps = U - tweak(U, n*d-n)  # TODO: make our code work with tweak 
    
    ev = cache[n*d-n]
    
    my_zero!(g)

    for j in axes(g,1)
        if g_exps == ev[j]
            CUDA.@allowscalar g[j] = lift(gCoeff)
            break
        end 
    end

    return (U,g)
end

"""
Takes an array of Costa's data tuples, and a new term
to be added. If the new term already exists, it is 
added, if not then it concatenates.

Note: perhaps this would be better served by a custom struct?

Note: right now, we don't have a method that just
takes an array of Costa's data tuples and consolodates it. 
Perhaps this also would be well suited by a custom struct?

^^ such concerns have speed implications, so better to wait
until we're really trying to optimize this.

"""
function incorporate_initial_term!(costadata_arr,costadata)
    ind_already = false

    # if the u is already there, just add the vectors
    for i in eachindex(costadata_arr)
        (u,g) = costadata_arr[i]
        if all(u .== costadata[1])
            costadata_arr[i] = (u, g .+ costadata[2])
            ind_already = true
        end
    end

    # otherwise, push on a new term
    if !ind_already
        #(9 < verbose) && println("incorporation: u has type $(typeof(costadata[1]))")
        push!(costadata_arr,costadata)
    end
end

function remove_duplicates_gpu!(costadata_arr)
    i = 1
    while i <= (length(costadata_arr)-1)
        j = i+1
        while j <= length(costadata_arr)
            if all(costadata_arr[i][1][1] .== costadata_arr[j][1][1])
                my_add!(costadata_arr[i][1][2],costadata_arr[i][1][2],costadata_arr[j][1][2])
                deleteat!(costadata_arr,j)
            else
                j = j + 1
            end
        end
        i = i + 1
    end
end

function remove_duplicates!(costadata_arr)
    i = 1
    while i <= (length(costadata_arr)-1)
        j = i+1
        while j <= length(costadata_arr)
            if all(costadata_arr[i][1][1] .== costadata_arr[j][1][1])
                costadata_arr[i] = ((costadata_arr[i][1][1], costadata_arr[i][1][2] + costadata_arr[j][1][2]),costadata_arr[i][2][1])
                deleteat!(costadata_arr,j)
            else
                j = j + 1
            end
        end
        i = i + 1
    end
end


"""
Converts a Costa's data tuple to a polynomial with pole
after reduction.

Thus, the pole order is n.
"""
function poly_of_end_costadata(costadata,PR,p,d,n,params)
    (u,g_vec) = costadata
    vars = gens(PR)
    cpu_g = Array(g_vec)
    CUDA.@allowscalar g = vector_to_polynomial(cpu_g,n,d*n-n,PR,params.termorder)

    return [prod(vars .^ u) * g, n]
end


"""
Converts an array of Costa's data tuples to an array of polynomials with pole
after reduction.

Thus, the pole order is always n.

"""
function poly_of_end_costadatas(costadatas,PR,p,d,n,S,params)
    res = PR(0)
    vars = gens(PR)
    for costadata in costadatas
        res += poly_of_end_costadata(costadata,PR,p,d,n,params)[1]
    end
    XS =  prod(PR(vars[i+1]) for i in S; init = PR(1))
    [[div(res,XS), n]]
end

function cuMod(A::zzModMatrix)
    m = modulus(base_ring(parent(A)))
    # if this conversion is illegal, 
    #then we're not allowed to make a CuModMatrix
    M = Int(m) 

    float_A = float_entries(A)
    CuModMatrix(float_A,M,elem_type=Float64)
end

function KaratsubaMat(A::zzModMatrix,A_temp=nothing)
    m = modulus(base_ring(parent(A)))
    temp = collect(factor(m))[1]
    d = temp[2]
    p = temp[1]

    N1 = Int(p)^Int(round(d/2))
    N2 = Int(p)^Int(d - round(d/2))

    GPUFiniteFieldMatrices.KaratsubaMatrix(Float64,cuMod(A),N1,N2,N1*N2)
end

function Karatsuba_copyto!(A_gpu::KaratsubaMatrix,A::Matrix{Float64},A_temp=nothing)
    
    temp = CuModMatrix(A,A_gpu.N1*A_gpu.N2,elem_type=Float64)

    GPUFiniteFieldMatrices.divide_elements!(A_gpu.data2,temp,A_gpu.N1)
    GPUFiniteFieldMatrices.mod_elements!(A_gpu.data2,A_gpu.N2)

    LinearAlgebra.mul!(A_gpu.data1,A_gpu.data2,Int(A_gpu.N1))
    GPUFiniteFieldMatrices.sub!(A_gpu.data1,temp,A_gpu.data1; mod_N=A_gpu.N1)
    GPUFiniteFieldMatrices.mod_elements!(A_gpu.data1,A_gpu.N1)

    A_gpu
end

#TODO: incorporate this idiomatically using adapt.jl
function adapt_to_gpu(Ruvs)
    gpu_Ruvs = Dict{Vector{Int64},Vector{CuModMatrix{Float64}}}()

    for (v,linpoly) in Ruvs

        gpu_Ruvs[v] = Vector(undef,length(linpoly)) 
        for i in 1:length(linpoly)
            gpu_Ruvs[v][i] = cuMod(linpoly[i])
        end
    end

    gpu_Ruvs
end

function default_context(matspace,Ruv,params)
    g_length = number_of_rows(matspace)
    B = base_ring(matspace)
    m = modulus(B)
    if params.use_gpu == true && (ZZ(2)^25 < m < ZZ(2)^106)
        println("using karatsuba")
        temp = collect(factor(m))[1]
        d = temp[2]
        p = temp[1]

        N1 = Int(p)^Int(round(d/2))
        N2 = Int(p)^Int((d - round(d/2)))

        A = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1*N2,true)
        B = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1*N2,true)
        temp = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1*N2,true)
        GPUFiniteFieldMatrices.initialize_plan!(A)
        GPUFiniteFieldMatrices.initialize_plan!(B)
        GPUFiniteFieldMatrices.initialize_plan!(temp)

        g = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,N1,N2,N1*N2,true)
        g_temp = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,N1,N2,N1*N2,true)
        GPUFiniteFieldMatrices.initialize_plan!(g)
        GPUFiniteFieldMatrices.initialize_plan!(g_temp)

        return ControlledReductionContext{KaratsubaMatrix{Float64},KaratsubaVector{Float64}}(Ruv,A,B,temp,g,g_temp)
    elseif params.use_gpu == true

        M = Int(m)
        A = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)
        B = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)
        temp = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)

        g = GPUFiniteFieldMatrices.zeros(Float64,g_length,M)
        g_temp = GPUFiniteFieldMatrices.zeros(Float64,g_length,M)

        return ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}(Ruv,A,B,temp,g,g_temp)
    else # use the CPU
        A = matspace()
        B = matspace()
        temp = matspace()

        if params.always_use_bigints || ZZ(2)^64 < ZZ(m)
            # Big modulus
            # NOTE: can't use `similar` for a pointer type
            g = [ZZ(0) for i in 1:g_length]
            g_temp = [ZZ(0) for i in 1:g_length]
        else
            # Small modulus
            g = zeros(UInt,g_length) 
            g_temp = similar(g)
        end

        return ControlledReductionContext(Ruv,A,B,temp,g,g_temp)
    end
end

"""
pregen_default_context
Note that n is (number of variables - 1) in this case
"""
function pregen_default_context(n,d,p,S;verbose=0, givefrobmat=false, algorithm=:naive, termorder=:invlex, vars_reversed=false, fastevaluation=false, always_use_bigints=false, use_gpu=false, use_threads=false,lazy=false)
    params = ZetaFunctionParams(verbose,givefrobmat,algorithm,termorder,vars_reversed,fastevaluation,always_use_bigints,use_gpu,use_threads)

    m = pregen_precision_info(n,d,p)
    M = Integer(p^m)
    g_length = binomial(d*n,d*n-n)
    matspace = matrix_space(residue_ring(ZZ,M)[1], g_length, g_length)

    Ruv = pregen_select_Ruv_PEP(n,d,S,params,matspace)

    if params.use_gpu == true && (ZZ(2)^25 < M < ZZ(2)^106)
        println("using karatsuba")

        N1 = Int(p)^Int(round(m/2))
        N2 = Int(p)^Int((m - round(m/2)))

        A = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1*N2,true)
        B = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1*N2,true)
        temp = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,g_length,N1,N2,N1*N2,true)
        GPUFiniteFieldMatrices.initialize_plan!(A)
        GPUFiniteFieldMatrices.initialize_plan!(B)
        GPUFiniteFieldMatrices.initialize_plan!(temp)

        g = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,N1,N2,N1*N2,true)
        g_temp = GPUFiniteFieldMatrices.KaratsubaZeros(Float64,g_length,N1,N2,N1*N2,true)
        GPUFiniteFieldMatrices.initialize_plan!(g)
        GPUFiniteFieldMatrices.initialize_plan!(g_temp)

        return ControlledReductionContext{KaratsubaMatrix{Float64},KaratsubaVector{Float64}}(Ruv,A,B,temp,g,g_temp)
    elseif params.use_gpu == true

        A = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)
        B = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)
        temp = GPUFiniteFieldMatrices.zeros(Float64,g_length,g_length,M)

        g = GPUFiniteFieldMatrices.zeros(Float64,g_length,M)
        g_temp = GPUFiniteFieldMatrices.zeros(Float64,g_length,M)

        return ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}(Ruv,A,B,temp,g,g_temp)
    else # use the CPU
        A = matspace()
        B = matspace()
        temp = matspace()

        if params.always_use_bigints || ZZ(2)^64 < ZZ(M)
            # Big modulus
            # NOTE: can't use `similar` for a pointer type
            g = [ZZ(0) for i in 1:g_length]
            g_temp = [ZZ(0) for i in 1:g_length]
        else
            # Small modulus
            g = zeros(UInt,g_length) 
            g_temp = similar(g)
        end

        return ControlledReductionContext(Ruv,A,B,temp,g,g_temp)
    end
end

function default_context_type(matspace,params)
    B = base_ring(matspace)
    m = modulus(B)

    if params.use_gpu == true
        ControlledReductionContext{CuModMatrix{Float64},CuModVector{Float64}}
    elseif params.always_use_bigints || ZZ(2)^64 < ZZ(m)
        ControlledReductionContext{ZZModMatrix,Vector{ZZRingElem}}
    else
        ControlledReductionContext{zzModMatrix,Vector{UInt}}
    end
end

"""
Gives a special PEP for when S = [0] and d=3

This should be all the Ruv ever needed in this situation.
"""
function cubic_S_zero_Vs(n)
    Vs = Vector{Vector{Int}}()

    for i in 1:n+1
        V = zeros(Int,n+1)
        V[i] = 3
        push!(Vs,V)
    end

    for i in 1:n
        V1 = zeros(Int,n+1)
        V1[i] = 1
        V1[i+1] = 2
        push!(Vs,V1)

        V2 = zeros(Int,n+1)
        V2[i] = 2
        V2[i+1] = 1
        push!(Vs,V2)
    end

    Vs
end

"""
This function exists so that we can have the option to time it
in verbose mode.
"""
function select_Ruv_PEP(n,d,S,params,compute,lazy,oscar_matspace,cache)

    m = Integer(modulus(base_ring(oscar_matspace)))
    
    compute_float = V -> float_entries.(compute(V))

    function compute_gpu(V; copyto=nothing) 

        if (3 < params.verbose) 
            println("Moving Ruv for $V to the gpu")
        end

        if (ZZ(2)^25 < m < ZZ(2)^106)
            if copyto != nothing
                Karatsuba_copyto!.(copyto, float_entries.(compute(V)))
            else
                KaratsubaMat.(compute(V))
            end
        else
            if copyto != nothing
                copyto!.(copyto, float_entries.(compute(V)))
                copyto
            else
                cuMod.(compute(V))
            end
        end 
    end


    EXTRA_MEMORY = false 

    if 5 == n && d == 3 && S == [n] && params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106) && !EXTRA_MEMORY
        println("hey howdy hey")
        
        Ruv = LRULazyPEP{KaratsubaMatrix{Float64}}(compute_gpu,d,usethreads=false)
    elseif 5 == n && d == 3 && S == [n] && params.use_gpu && !EXTRA_MEMORY
        println("howdy")

        Ruv = LRULazyPEP{CuModMatrix{Float64}}(compute_gpu,d,usethreads=false)

    elseif 3 < n && d == 3 && S == [n] && params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106) #4 < n
        println("hey")

        eager_Vs = cubic_S_zero_Vs(n)

        cpu_Ruv = LazyPEP{Matrix{Float64}}(compute_float,eagerVs=eager_Vs,usethreads=false)
        m = M = Int(modulus(base_ring(oscar_matspace)))
        temp = collect(factor(m))[1]
        d = temp[2]
        p = temp[1]

        N1 = Int(p)^Int(round(d/2))
        N2 = Int(p)^Int(d - round(d/2))

        if (3 < params.verbose)
            create_gpu = matrices -> begin
                println("Moving an Ruv to the GPU (first time creation!)")
                CUDA.@time KaratsubaMatrix.(Float64,CuModMatrix.(matrices,N1*N2,elem_type=Float64),N1,N2,N1*N2)
            end

            convert_gpu = (dest,matrices) -> begin
                println("Moving an Ruv to the GPU (already allocated)")

                CUDA.@time begin
                    for i in 1:length(dest)
                        Karatsuba_copyto!(dest[i],matrices[i])
                    end
                end
            end
        else
            # should it really be N1??? or actually N1 * N2 TODO
            create_gpu = matrices -> KaratsubaMatrix.(Float64,CuModMatrix.(matrices,N1*N2,elem_type=Float64),N1,N2,N1*N2)
            convert_gpu = (dest,matrices) -> begin
                for i in 1:length(dest)
                    Karatsuba_copyto!(dest[i],matrices[i])
                end
            end
        end
        s = size(oscar_matspace(),1)

        # memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
        # # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
        # memory = s^2 * 8 * (n+2)

        # maxsize = div(memory_cap,memory)

        #println(maxsize)
        #testing
        if n == 5 # (cubic) fourfold
            maxsize = 3 
        elseif n == 4 #(cubic) threefold
            maxsize = 3#20 
        else
            memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
            # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
            memory = s^2 * 8 * (n+2)

            maxsize = div(memory_cap,memory)
        end

        Ruv = CachePEP{Matrix{Float64},KaratsubaMatrix{Float64}}(cpu_Ruv,create_gpu,convert_gpu,maxsize)

    elseif 3 < n && d == 3 && S == [n] && params.use_gpu #4 < n
        println("hello")

        eager_Vs = cubic_S_zero_Vs(n)

        cpu_Ruv = LazyPEP{Matrix{Float64}}(compute_float,eagerVs=eager_Vs,usethreads=false)
        m = M = Int(modulus(base_ring(oscar_matspace)))

        if (3 < params.verbose)
            create_gpu = matrices -> begin
                println("Moving an Ruv to the GPU (first time creation!)")
                CUDA.@time CuModMatrix.(matrices,M,elem_type=Float64)
            end

            convert_gpu = (dest,matrices) -> begin
                println("Moving an Ruv to the GPU (already allocated)")

                CUDA.@time begin
                    for i in 1:length(dest)
                        copyto!(dest[i],matrices[i])
                    end
                end
            end
        else
            create_gpu = matrices -> CuModMatrix.(matrices,M,elem_type=Float64)
            convert_gpu = (dest,matrices) -> begin
                for i in 1:length(dest)
                    copyto!(dest[i],matrices[i])
                end
            end
        end
        s = size(oscar_matspace(),1)

        if n == 5 # (cubic) fourfold
            maxsize = 3 # really should probably make this an LRU cache for varbyvar
        elseif n == 4 #(cubic) threefold
            maxsize = 20 
        else
            memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
            # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
            memory = s^2 * 8 * (n+2)

            maxsize = div(memory_cap,memory)
        end
        #maxsize = 13 

        Ruv = CachePEP{Matrix{Float64},CuModMatrix{Float64}}(cpu_Ruv,create_gpu,convert_gpu,maxsize)

        #(1 < params.verbose) && println("Initial cache info: $(cache_info(Ruv.Ucomponent))")
    elseif params.use_gpu && lazy && (ZZ(2)^25 < m < ZZ(2)^106)
        compute_gpu_karatsuba = V -> @. KaratsubaMat(float_entries(compute(V)))
        Ruv = LazyPEP{KaratsubaMatrix{Float64}}(compute_gpu_karatsuba)

    elseif params.use_gpu && lazy
        Ruv = LazyPEP{CuModMatrix{Float64}}(compute_gpu)

    elseif params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106)
    compute_gpu_karatsuba = V -> @. KaratsubaMat(float_entries(compute(V)))
        Ruv = EagerPEP{KaratsubaMatrix{Float64}}(cache[d],compute_gpu_karatsuba,usethreads=false)

    elseif params.use_gpu
        Ruv = EagerPEP{CuModMatrix{Float64}}(cache[d],compute_gpu,usethreads=params.use_threads)

    elseif lazy
        Ruv = LazyPEP{typeof(oscar_matspace())}(compute)

    else
        Ruv = EagerPEP{typeof(oscar_matspace())}(cache[d],compute,usethreads=params.use_threads)
    end

    Ruv
end

function pregen_select_Ruv_PEP(n,d,S,params,oscar_matspace)

    m = Integer(modulus(base_ring(oscar_matspace)))

    if 3 < n && d == 3 && S == [n] && params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106) #4 < n

        cpu_Ruv = PregenLazyPEP{Matrix{Float64}}(nothing)
        m = M = Int(modulus(base_ring(oscar_matspace)))
        temp = collect(factor(m))[1]
        d = temp[2]
        p = temp[1]

        N1 = Int(p)^Int(round(d/2))
        N2 = Int(p)^Int(d - round(d/2))

        create_gpu = matrices -> KaratsubaMatrix.(Float64,CuModMatrix.(matrices,N1,elem_type=Float64),N1,N2,N1*N2)
        convert_gpu = (dest,matrices) -> begin
            for i in 1:length(dest)
                copyto!(dest[i],matrices[i])
            end
        end
        s = size(oscar_matspace(),1)

        memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
        # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
        memory = s^2 * 8 * (n+2)

        maxsize = div(memory_cap,memory)

        #println(maxsize)
        #testing
        #maxsize = 13 

        Ruv = CachePEP{Matrix{Float64},KaratsubaMatrix{Float64}}(cpu_Ruv,create_gpu,convert_gpu,maxsize)

    elseif 3 < n && d == 3 && S == [n] && params.use_gpu #4 < n

        cpu_Ruv = PregenLazyPEP{Matrix{Float64}}(nothing)
        m = M = Int(modulus(base_ring(oscar_matspace)))

        create_gpu = matrices -> CuModMatrix.(matrices,M,elem_type=Float64)
        convert_gpu = (dest,matrices) -> begin
            for i in 1:length(dest)
                copyto!(dest[i],matrices[i])
            end
        end
        s = size(oscar_matspace(),1)

        memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
        # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
        memory = s^2 * 8 * (n+2)

        maxsize = div(memory_cap,memory)

        #println(maxsize)
        #testing
        #maxsize = 13 

        Ruv = CachePEP{Matrix{Float64},CuModMatrix{Float64}}(cpu_Ruv,create_gpu,convert_gpu,maxsize)

        #(1 < params.verbose) && println("Initial cache info: $(cache_info(Ruv.Ucomponent))")
    
    elseif params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106)
        Ruv = PregenLazyPEP{KaratsubaMatrix{Float64}}(nothing)
    elseif params.use_gpu 
        Ruv = PregenLazyPEP{CuModMatrix{Float64}}(nothing)
    else
        Ruv = PregenLazyPEP{typeof(oscar_matspace())}(nothing)
    end

    Ruv
end

function reducetransform(FT,N_m,S,f,pseudoInverseMat,p,params,cache,context)
    if params.algorithm == :pchunk
        reducetransform_pchunk(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    elseif params.algorithm == :depthfirst
        reducetransform_depthfirst(FT,N_m,S,f,pseudoInverseMat,p,cache,params,context)
    elseif params.algorithm == :varbyvar
        reducetransform_varbyvar(FT,N_m,S,f,pseudoInverseMat,p,cache,params,context)
    elseif params.algorithm == :akr
        reducetransform_akr(FT,N_m,S,f,pseudoInverseMat,p,cache,params)
    else
        println(params.algorithm)
        throw(ArgumentError("Unsupported Algorithm: "))
    end
end

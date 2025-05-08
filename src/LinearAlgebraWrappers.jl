
"""
KeyType - the type of u and v, for example Vector{Int64}
MatrixType - the type of A and B, for example zzModMatrix
VectorType - the type of g, for exampel Vector{UInt64}


#const OscarSmallReductionContext = ControlledReductionContext{Vector{Int64},
#                                                         zzModMatrix,
#                                                         Vector{UInt64}}
#const OscarBigReductionContext = ControlledReductionContext{Vector{Int64},
#                                                         ZZModMatrix,
#                                                         Vector{ZZRingElem}}

In the future: make a GPUReductionContext

"""
struct ControlledReductionContext{KeyType,MatrixType,VectorType}
    Ruvs::Dict{KeyType, Vector{MatrixType}} #TODO: this type might be too rigid
    A::MatrixType
    B::MatrixType
    temp::MatrixType
    g::VectorType
    g_temp::VectorType
end

### The following is derived from Nemo.jl, and needs to be licensed with the GPL 

#TODO: this isn't used right now, but it might be faster to use it
function my_addmul!(A::zzModMatrix, B::zzModMatrix, C::UInt, D::zzModMatrix)
  @ccall Oscar.Nemo.libflint.nmod_mat_scalar_addmul_ui(A::Ref{zzModMatrix}, 
                                                       B::Ref{zzModMatrix}, 
                                                       D::Ref{zzModMatrix}, 
                                                       C::UInt)::Nothing
  return A
end

function my_mul!(A::zzModMatrix, B::zzModMatrix, c)
    n = characteristic(base_ring(parent(A)))
    ui(i) = 0 ≤ i ? UInt(i % n) : UInt((i % n) + n)
    c = ui(c)
    @ccall Oscar.Nemo.libflint.nmod_mat_scalar_mul(A::Ref{zzModMatrix},
                                                   B::Ref{zzModMatrix},
                                                   c::UInt)::Nothing
    return A
end

function my_mul!(A::ZZModMatrix, B::ZZModMatrix, c::Int)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_scalar_mul_si(A::Ref{ZZModMatrix}, 
                                               B::Ref{ZZModMatrix}, 
                                               c::Int, 
                                               base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return A 
end

function my_mul!(A::ZZModMatrix, B::ZZModMatrix, c::ZZRingElem)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_scalar_mul_fmpz(A::Ref{ZZModMatrix}, 
                                                 B::Ref{ZZModMatrix},
                                                 C::Ref{ZZRingElem},
                                                 base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return A
end

function my_matvecmul!(z::Vector{UInt},A::zzModMatrix,b::Vector{UInt})
    @ccall Oscar.Nemo.libflint.nmod_mat_mul_nmod_vec(z::Ptr{UInt},
                                                     A::Ref{zzModMatrix},
                                                     b::Ptr{UInt},
                                                     length(b)::Int)::Nothing
    return z
end

function my_matvecmul!(z::Vector{ZZRingElem},A::ZZModMatrix,b::Vector{ZZRingElem})
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_mul_fmpz_vec_ptr(z::Ptr{Ref{ZZRingElem}}, 
                                                  A::Ref{ZZModMatrix}, 
                                                  b::Ptr{Ref{ZZRingElem}}, 
                                                  length(b)::Int, 
                                                  base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return z
end

#function alt_matvecmul!(z::Vector{ZZRingElem},A::ZZModMatrix,b::Vector{ZZRingElem})
#    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_mul_fmpz_vec(z::Ref{ZZRingElem},
#                                                         A::Ref{ZZModMatrix},
#                                                         b::Ref{ZZRingElem},
#                                                         length(b)::Int,
#                                                         base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
#end

### END stuff derived from Nemo.jl

# (scalar) mul

my_mul!(A::CuModMatrix,B::CuModMatrix,c::Number) = GPUFiniteFieldMatrices.mul!(A,B,c) 

# matvecmul

my_matvecmul!(z::CuModVector,A::CuModMatrix,b::CuModVector) = GPUFiniteFieldMatrices.mul!(z,A,b)

# copy

function my_copy!(a,b)
    copy!(a,b)
end

function my_copy!(a::Vector{ZZRingElem},b::Vector{ZZRingElem})
    for i in 1:length(a)
        Oscar.Nemo.set!(a[i],b[i])
    end
end

# zero

function my_zero!(a)
    fill!(a,zero(eltype(a)))
    # apparently, the following is not implemented in Julia, it's an Oscar method
    #zero!(a)
end

my_zero!(a::zzModMatrix) = zero!(a)
my_zero!(a::ZZModMatrix) = zero!(a)
my_zero!(a::CuModArray) = GPUFiniteFieldMatrices.zero!(a)

function my_zero!(a::Vector{ZZRingElem})
    # this is here because of BigInts/ZZRingElem behaving weird in Julia

    for i in 1:length(a)
        Oscar.Nemo.set!(a[i],0)
    end
end

# sub

function my_sub!(A,B,C)
    Oscar.Nemo.sub!(A,B,C)
end

my_sub!(A::CuModMatrix,B::CuModMatrix,C::CuModMatrix) = GPUFiniteFieldMatrices.sub!(A,B,C)

function my_sub!(A::ZZModMatrix,B::ZZModMatrix,C::ZZModMatrix)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_sub(A::Ref{ZZModMatrix},
                                                B::Ref{ZZModMatrix},
                                                C::Ref{ZZModMatrix},
                                                base_ring(B).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct})::Nothing
    return A
end

# add

my_add!(A,B,C) = add!(A,B,C)
my_add!(A::CuModMatrix,B::CuModMatrix,C::CuModMatrix) = GPUFiniteFieldMatrices.add!(A,B,C)

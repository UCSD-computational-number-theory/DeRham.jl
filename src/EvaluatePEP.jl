
"""
    eval_to_linear!(A,B,temp,matrices,U,V)

Takes a list of n+2 matrices and ouputs a list two matrices [A,B] corresponding to R_{(x0,...,xn)+yv, v} = Ay + B

Takes a polynomial F in length(U) variables with matrix coefficients,
and evaluates it at U + y*V, getting a linear polynomial with 
matrix coefficients A*y + B.

INPUTS: 
* "A" -- output matrix A to be modified in place
* "B" -- output matrix B to be modified in place
* "matrices" -- list of size length(u) + 1, the coefficients of F.
    The constant term is matrices[1] and the term corresponding to 
    u[i] and v[i] is matrices[i+1].
* "U" -- vector, U =(x0, ..., xn)
* "V" -- vector, V
* "R" -- ring, base ring of f
* temp -- temporary matrix to be modified in place
"""
function eval_to_linear!(A,B,temp,matrices,U,V)
    my_zero!(B)
    my_add!(B,B,matrices[1])

    my_zero!(A)

    for k in 2:(length(matrices))

        #B = B + matrices[k] * U[k-1]
        #A = A + matrices[k] * V[k-1]

        my_mul!(temp,matrices[k],U[k-1])
        my_add!(B,B,temp)
        my_mul!(temp,matrices[k],V[k-1])
        my_add!(A,A,temp)
    end 

    return (A, B)
end

function eval_to_linear!(A,B,temp,pep::AbstractPEP,U,V)
    matrices = pep[V]
    eval_to_linear!(A,B,temp,matrices,U,V)
end


# use linear indexing since this is component-wise
function eval_to_linear_kernel!(A,B,matrices,U,V,M)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    bb = matrices[1][i]
    aa = 0

    # this assumes we can do n operations without overflow.
    for k in 2:length(matrices)

        m = matrices[k][i]

        bb += m * U[k-1]
        aa += m * V[k-1]
    end

    A[i] = mod(aa, M)
    B[i] = mod(bb, M)

    nothing
end

function eval_to_linear_gpu!(A,B,temp,matrices,U,V)

    tw = GPUFiniteFieldMatrices.TILE_WIDTH

    # threads = (tw,tw)
    threads = tw

    # totalsize = size(A.data)

    blocks = length(A.data) ÷ threads
    # blocks = (totalsize[1] ÷ threads[1], totalsize[2] ÷ threads[2])

    mat_tuple = tuple(map(A -> A.data, matrices)...)

    @cuda threads=threads blocks=blocks eval_to_linear_kernel!(A.data,
                                               B.data,
                                               mat_tuple,
                                               tuple(U...),
                                               tuple(V...),
                                              A.N)

end

function eval_to_linear_karatsuba_kernel!(Adata1,Adata2,Bdata1,Bdata2,matrices1,matrices2,U,V,N1,N2)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    bb1 = matrices1[1][i]
    bb2 = matrices2[1][i]

    aa1 = 0
    aa2 = 0

    for k in 2:length(matrices1)

        m1 = matrices1[k][i]
        m2 = matrices2[k][i]

        temp1, temp2 = GPUFiniteFieldMatrices.karatsuba_scalar_multiply_helper(m1,m2,U[k-1],N1,N2)
        bb1, bb2 = GPUFiniteFieldMatrices.karatsuba_add_helper(bb1,bb2,temp1,temp2,N1,N2)

        temp1, temp2 = GPUFiniteFieldMatrices.karatsuba_scalar_multiply_helper(m1,m2,V[k-1],N1,N2)
        aa1, aa2 = GPUFiniteFieldMatrices.karatsuba_add_helper(aa1,aa2,temp1,temp2,N1,N2)
    end

    Bdata1[i] = mod(bb1, N1)
    Bdata2[i] = mod(bb2, N2)

    Adata1[i] = mod(aa1, N1)
    Adata2[i] = mod(aa2, N2)

    nothing
end

function eval_to_linear_gpu_karatsuba!(A,B,temp,matrices,U,V)

    tw = GPUFiniteFieldMatrices.TILE_WIDTH

    # threads = (tw,tw)
    threads = tw

    # totalsize = size(A.data)

    blocks = length(A.data1.data) ÷ threads
    # blocks = (totalsize[1] ÷ threads[1], totalsize[2] ÷ threads[2])

    mat_tuple_1 = tuple(map(A -> A.data1.data, matrices)...)
    mat_tuple_2 = tuple(map(A -> A.data2.data, matrices)...)

    @cuda threads=threads blocks=blocks eval_to_linear_karatsuba_kernel!(A.data1.data,
                                                                         A.data2.data,
                                                                         B.data1.data,
                                                                         B.data2.data,
                                                                         mat_tuple_1,
                                                                         mat_tuple_2,
                                                                         tuple(U...),
                                                                         tuple(V...),
                                                                         A.N1,
                                                                         A.N2)

end

"""
    finitediff_prodeval_linear!(a,b,start,stop,g,temp)

Generic function that evaluates 
the function F(x) = a*x + b at 
start, ..., stop and 
computes
F(start)*...*F(stop)*g

Since it's a generic function, it'll work if
a and b and g are numbers, but the intended 
use is for a and b to be matrices and g
to be a vector.

My hope is that since this is generic, we'll
be able to replace arrays with CuArrays for large examples
and it'll still work.

a - linear term coefficient
b - constant term coefficient
start - lowest value to evaluate at (should be an integer)
stop - highest value to evaluate at (should be an integer)
g - the value into which the answer is accumulated.
temp - a temporary pre-allocated matrix that can be used in intermediate computations
ui - a function that converts the output of Julia's mod into an unsigned integer 
    (i.e. if it's negative, add the characteristic)

"""
function finitediff_prodeval_linear!(a,b,start,stop,g,temp,g_temp,ui=nothing)

  my_mul!(temp,a,stop)
  my_add!(temp,temp,b)

  if start == stop
    my_matvecmul!(g_temp,temp,g)
    my_copy!(g,g_temp)

    return g
  end
  my_matvecmul!(g_temp,temp,g)
  my_copy!(g,g_temp)


  for k = stop-1:-1:start
    # right now, Fk = F(k+1)
    
    # Fk = Fk - a
    my_sub!(temp,temp,a)
    my_matvecmul!(g_temp,temp,g)
    my_copy!(g,g_temp)
  end
  return g
end

# function finitediff_prodeval_linear!(a::KaratsubaArray,b::KaratsubaArray,start,stop,g::KaratsubaArray,temp::KaratsubaArray,g_temp::KaratsubaArray,ui=nothing)

#   my_mul!(temp,a,stop)
#   my_add!(temp,temp,b)

#   if start == stop
#     my_matvecmul!(g_temp,temp,g)
#     my_copy!(g,g_temp)

#     return g
#   end
#   my_matvecmul!(g_temp,temp,g)
#   my_copy!(g,g_temp)

#   my_copy!(b,temp)
#   i = 1
#   for k = stop-1:-1:start
#     # right now, Fk = F(k+1)
    
#     # Fk = Fk - a
#     if i == 1
#       my_sub!(temp,b,a)
#       my_matvecmul!(g_temp,temp,g)
#     else
#       my_sub!(b,temp,a)
#       my_matvecmul!(g_temp,b,g)
#     end
#     my_copy!(g,g_temp)
#     i = (i + 1)%2
#   end
#   return g
# end

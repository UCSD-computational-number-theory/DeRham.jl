
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

function finitediff_prodeval_linear!(a::KaratsubaArray,b::KaratsubaArray,start,stop,g::KaratsubaArray,temp::KaratsubaArray,g_temp::KaratsubaArray,ui=nothing)

  my_mul!(temp,a,stop)
  my_add!(temp,temp,b)

  if start == stop
    my_matvecmul!(g_temp,temp,g)
    my_copy!(g,g_temp)

    return g
  end
  my_matvecmul!(g_temp,temp,g)
  my_copy!(g,g_temp)

  my_copy!(b,temp)
  i = 1
  for k = stop-1:-1:start
    # right now, Fk = F(k+1)
    
    # Fk = Fk - a
    if i == 1
      my_sub!(temp,b,a)
      my_matvecmul!(g_temp,temp,g)
    else
      my_sub!(b,temp,a)
      my_matvecmul!(g_temp,b,g)
    end
    my_copy!(g,g_temp)
    i = (i + 1)%2
  end
  return g
end
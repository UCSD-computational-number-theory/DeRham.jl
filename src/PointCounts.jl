
"""
Extracts the point counts from the zeta function
"""
function pointcounts(n,d,zeta,p,N)
  if n != 2 || d != 3
    error("not implement for anything but elliptic curves")
  end

  P, t = power_series_ring(ZZ,N,:t)

  #TODO: 
  #
  # * figure out the formula for the extra terms in the 
  #   zeta function
  # * evaluate log(zeta) here
  # * extract the coefficients
  # * multiply by n
  # * return the point counts
  

end

"""
Counts points in the most naive way possible

Returns a tuple containing
a list of the points of the affine
variety k[x_1,...,x_n]/f, 
and the number of projective points
"""
function naivelypointcount(f,q)
  if !is_prime(q)
    error("naive point counting not implemented for non-prime fields yet")
  end

  rationalpoints = []

  n = length(gens(parent(f)))
  iter = Iterators.product(fill(0:q-1,n)...)
  for xs in iter
    val = evaluate(f,collect(xs))
    if val == 0
      push!(rationalpoints,xs)
    end
  end

  n_aff_points = length(rationalpoints)
  (rationalpoints,(n_aff_points - 1)/(q-1))
end

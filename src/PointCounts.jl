
"""
Extracts the point counts from the zeta function
"""
function pointcount(n,d,zeta,p,q)
  if n != 2 || d != 3 || q != p
    error("not implement for anything but elliptic curves")
  end


  # right now this if will always hit but in the future
  # we will implement this for more general things
  if n == 2 && d == 3 && q == p
    t = gens(parent(zeta))[1]
    second_coef = coeff(zeta,1)
    # zeta[2] = -a_p
    # a_p =  p + 1 - #E(F_p)
    return second_coef + p + 1
  end

  # MARK - general case
  #P, t = power_series_ring(ZZ,N,:t)

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

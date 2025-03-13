
function has_good_reduction_fermat(p)
  R, (x,y,z) = polynomial_ring(GF(p),3)
  f = x^3 + y^3 + z^3
  DeRham.check_smoothness(f)
end

using Primes

function get_fermat_bad_primes()
  ps = nextprimes(1,20)
  for p in ps
    if !has_good_reduction_fermat(p)
      println(p)
    end
  end
end

function test_fermat_cubic_naive()
  ps = [3
  5
  7
 11
 13
 17
 19
 23
 29
 31
 37
 41
 43
 47
 53
 59
 61
 67
 71]
  for p in ps
    R, (x,y,z) = polynomial_ring(GF(p),3)
    f = x^3 + y^3 + z^3

    naive = DeRham.naivelypointcount(f,p)


    n = 2
    d = 3
    zeta = DeRham.zeta_function(f)
    println(zeta)
    pc = DeRham.pointcount(n,d,zeta,p,p)
    #pc2 = DeRham.pointcount(n,d,zeta[2],p,p)
    #println()
    #println("p = $p")
    #print("Naive: $(naive[2]), ")
    #println("Controlled: $(pc), $(zeta)")
    ##println("Other controlled: $(pc2), $(zeta[2])")
    #println()

    @test naive[2] == pc

  end
end

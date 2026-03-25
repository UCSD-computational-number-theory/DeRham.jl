"""
    calculate_precision_information(f,basis)

returns all the precision information
related to f

basis -- the basis of cohomology, in the format of
         polynomials with pole
"""
function calculate_precision_information(f, basis, verbose=0)
    p = Int64(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    hodge_polygon = hodgepolygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths

    k = sum(hodge_numbers)
    (0 < verbose) && println("There are $k basis total elements in H^$n to be reduced")

    r_m = calculate_relative_precision(hodge_polygon, n - 1, p)
    N_m = series_precision(p, n, d, r_m)
    N_m[N_m .== 0] .= 1
    M = algorithm_precision(p, n, d, r_m, N_m)

    (hodge_polygon, r_m, N_m, M)
end

"""
    calculate_precision_information(f)

returns the precision information for a polynomial f
"""
function calculate_precision_information(f)
    Basis = get_basis_of_cohomology(f)

    calculate_precision_information(f, Basis)
end

function precision_information_auto(f; basis=nothing, verbose=0)
    if basis === nothing
        return calculate_precision_information(f)
    end
    return calculate_precision_information(f, basis, verbose)
end

function precision_information_auto(f, field; basis=nothing, verbose=0)
    if characteristic(field) != characteristic(parent(f))
        error("Field characteristic $(characteristic(field)) does not match polynomial base field characteristic $(characteristic(parent(f))).")
    end
    return precision_information_auto(f; basis=basis, verbose=verbose)
end

function zeta_function_auto(f; basis=nothing, verbose=0, kwargs...)
    precision_info = precision_information_auto(f; basis=basis, verbose=verbose)
    return zeta_function(f; precision_info=precision_info, verbose=verbose, kwargs...)
end

function precision_information_user_r(f, r; basis=nothing, verbose=0)
    n = nvars(parent(f)) - 1
    r_m = fill(Int(r), n)
    return precision_information_user_r_m(f, r_m; basis=basis, verbose=verbose)
end

function precision_information_user_r_m(f, r_m; basis=nothing, verbose=0)
    r_m_int = [Int(x) for x in r_m]

    p = Int64(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    if basis === nothing
        Basis = get_basis_of_cohomology(f)
    else
        Basis = basis
    end

    hodge_polygon = hodgepolygon(Basis, n)
    hodge_numbers = hodge_polygon.slopelengths
    if length(r_m_int) != length(hodge_numbers)
        error("r_m must have length $(length(hodge_numbers)) for this basis.")
    end

    N_m = series_precision(p, n, d, r_m_int)
    N_m[N_m .== 0] .= 1
    M = algorithm_precision(p, n, d, r_m_int, N_m)

    return (hodge_polygon, r_m_int, N_m, M)
end

function precision_information_max_auto_r(f, r; basis=nothing, verbose=0)
    n = nvars(parent(f)) - 1
    r_m = fill(Int(r), n)
    return precision_information_max_auto_r_m(f, r_m; basis=basis, verbose=verbose)
end

function precision_information_max_auto_r_m(f, r_m; basis=nothing, verbose=0)
    (hodge_polygon, r_m_auto, _, _) = precision_information_auto(f; basis=basis, verbose=verbose)
    (_, r_m_user, _, _) = precision_information_user_r_m(f, r_m; basis=basis, verbose=verbose)

    p = Int64(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    r_m_max = max.(r_m_auto, r_m_user)
    N_m = series_precision(p, n, d, r_m_max)
    N_m[N_m .== 0] .= 1
    M = algorithm_precision(p, n, d, r_m_max, N_m)

    return (hodge_polygon, r_m_max, N_m, M)
end

"""
   precision_information(f,basis)

returns all the precision information
related to f

basis -- the basis of cohomology, in the format of 
         polynomials with pole
"""
function precision_information(f,basis,verbose=0)
    p = Int64(characteristic(parent(f)))
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    hodge_polygon = hodgepolygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths

    k = sum(hodge_numbers) # dimension of H^n
    (0 < verbose) && println("There are $k basis total elements in H^$n to be reduced")

    r_m = calculate_relative_precision(hodge_polygon, n-1, p) 
    N_m = series_precision(p,n,d,r_m)
    N_m[N_m .== 0] .= 1 #TODO: understand the meaning of zero precision better
    M = algorithm_precision(p,n,d,r_m,N_m)

    (hodge_polygon,r_m,N_m,M)
end

"""
    precision_information(f)

returns the precision information for a polynomial f

"""
function precision_information(f)
    n = nvars(parent(f)) - 1
    PR = parent(f)
    R = coefficient_ring(parent(f))

    params = default_params()

    d = total_degree(f)
    S = collect(0:n) #collect(0:d-1)
    cache = controlled_reduction_cache(n,d,S,params)
    basis = compute_monomial_bases(f, params, cache) # basis of cohomology 
    Basis = []
    for i in 1:n
        for j in basis[i]
            push!(Basis,[j,i])
        end
    end

    precision_information(f,Basis)
end

"""
    print_precision_info(n,d,p)

Calcuates how much precision would be needed to compute
the zeta functino of a polynomial with
n+1 variables, of degree d, at the prime p
"""
function print_precision_info(n,d,p)
    R, vars = polynomial_ring(GF(p),n+1)

    # any hypersurface will do, all we care
    # about is the hodge polygon
    fermat_hypersurface = sum(vars .^ d)

    (hp, rel, ser, alg) = calculate_precision_information(fermat_hypersurface)

    println("Hodge Numbers: $(hp.slopelengths)")
    println("Relative precision: $rel")
    println("Series precision: $ser")
    println("Algorithm precision: $alg")
end 

function pregen_precision_info(n,d,p)
    R, vars = polynomial_ring(GF(p),n+1)

    # any hypersurface will do, all we care
    # about is the hodge polygon
    fermat_hypersurface = sum(vars .^ d)

    (hp, rel, ser, alg) = calculate_precision_information(fermat_hypersurface)

    return alg
end 

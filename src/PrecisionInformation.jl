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

function zeta_function_with_precision(f, r::Integer; basis=nothing, verbose=0, kwargs...)
    r_m = fill(Int(r), nvars(parent(f)) - 1)
    return zeta_function_with_precision(f, r_m; basis=basis, verbose=verbose, kwargs...)
end

function zeta_function_with_precision(f, r::AbstractVector; basis=nothing, verbose=0, kwargs...)
    precision_info = precision_information_max_auto_r_m(f, r; basis=basis, verbose=verbose)
    return zeta_function(f; precision_info=precision_info, verbose=verbose, kwargs...)
end

function frobenius_matrix_with_precision(f, r::Integer; basis=nothing, verbose=0, kwargs...)
    r_m = fill(Int(r), nvars(parent(f)) - 1)
    return frobenius_matrix_with_precision(f, r_m; basis=basis, verbose=verbose, kwargs...)
end

function frobenius_matrix_with_precision(f, r::AbstractVector; basis=nothing, verbose=0, kwargs...)
    precision_info = precision_information_user_r_m(f, r; basis=basis, verbose=verbose)
    result = zeta_function(f; precision_info=precision_info, givefrobmat=true, verbose=verbose, kwargs...)
    if result == false
        return false
    end
    return result[1]
end

function hasse_witt_matrix(f; basis=nothing, verbose=0, kwargs...)

    n = nvars(parent(f)) - 1
    d = total_degree(f)
    p = characteristic(parent(f))

    if d < n
        throw(ArgumentError("Hasse-Witt matrices are no supported for varieties of negative Kodaira dimensin (i.e. degree < number of variables"))
    end

    r_m = fill(0,n)
    r_m[1] = 1

    basis = get_basis_of_cohomology(f)

    hodge_polygon = hodgepolygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths
    fh = hodge_numbers[1]

    F = frobenius_matrix_with_precision(f, r_m, basis=basis, verbose=verbose, kwargs...)

    highprec = F[end-fh+1:end,end-fh+1:end]

    K = GF(p)

    map_entries(K, lift.(highprec))
end

function a_number(f; basis=nothing, verbose=0, kwargs...)
    n = nvars(parent(f)) - 1

    if n != 2
        throw(ArgumentError("a-numbers are only supported for curves so far"))
    end

    if 0 < verbose
        println("Calculating Hasse-Witt matrix...")
    end
    HW = hasse_witt_matrix(f; basis=nothing, verbose=0, kwargs...)

    basis = get_basis_of_cohomology(f)

    hodge_polygon = hodgepolygon(basis, n)
    hodge_numbers = hodge_polygon.slopelengths
    g = hodge_numbers[1]

    g - rank(HW)
end

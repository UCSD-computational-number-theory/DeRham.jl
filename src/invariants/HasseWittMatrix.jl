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
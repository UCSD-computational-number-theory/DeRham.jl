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
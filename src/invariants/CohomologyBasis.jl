function get_basis_of_cohomology(f)
    n = nvars(parent(f)) - 1
    d = total_degree(f)

    params = default_params()
    S = collect(0:n)
    cache = controlled_reduction_cache(n, d, S, params)

    get_basis_of_cohomology(f,S,params,cache)
end
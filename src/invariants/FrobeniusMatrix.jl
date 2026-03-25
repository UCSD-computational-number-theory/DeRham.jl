function frobenius_matrix_with_precision(f, r::Integer; basis=nothing, verbose=0, kwargs...)
    r_m = fill(Int(r), nvars(parent(f)) - 1)
    return frobenius_matrix_with_precision(f, r_m; basis=basis, verbose=verbose, kwargs...)
end

function frobenius_matrix_with_precision(f, r::AbstractVector; basis=nothing, verbose=0, kwargs...)
    precision_info = precision_information_user_r_m(f, r; basis=basis, verbose=verbose)
    result = zeta_coefficients(f; precision_info=precision_info, givefrobmat=true, verbose=verbose, kwargs...)
    if result == false
        return false
    end
    return result[1]
end

module Combined
using Oscar
include("ControlledReduction.jl")



function frobenius_matrix_classical(f, n, d, precision, p)
    # Step 1: find monomial basis
    PR = parent(f)
    R = coefficient_ring(PR) # function name?
    
    # Step 2: apply Frobenius

    # Step 3: reduction

    # Step 4: process the answer
end

function frobenius_matrix_controlled(f, n, d, precision, p)

end



end
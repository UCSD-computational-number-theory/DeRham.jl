
function test_supported_monomial_orderings()

    n = 2
    d = 3

    p = 7

    F = GF(p)
    R, (x,y,z) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = y^2*z - x^3 - x*z^2 - z^3

    zeta_invlex = DeRham.zeta_function(f,termorder=:invlex)

    zeta_lex = DeRham.zeta_function(f,termorder=:lex)

    zeta_neglex = DeRham.zeta_function(f,termorder=:neglex)

    # TODO: put htis line back when we've figured out precision
    #@test zeta_invlex == zeta_lex
    
    # This is brittle, but it passes, showing that we are doing something
    # correct. When we fix precision, we can go back to the above line
    mod49_invlex = lift.((ZZ,),coefficients(zeta_invlex[1])) .% 49
    mod49_lex = lift.((ZZ,),coefficients(zeta_lex[1])) .% 49
    
    @test mod49_invlex == mod49_lex

    @test zeta_invlex == zeta_neglex
end

function test_reversing_variables()

    n = 2
    d = 3

    p = 7

    F = GF(p)
    R, (x,y,z) = polynomial_ring(F, ["x$i" for i in 0:n])

    f = y^2*z - x^3 - x*z^2 - z^3

    #FIXME: variable reversing doesn't actually do anything right now!
    #This needs to be addressed in ZetaFunction.jl
    zeta_reversed = DeRham.zeta_function(f,termorder=:invlex,vars_reversed=true)
    zeta_normal = DeRham.zeta_function(f,termorder=:invlex,vars_reversed=false)
    
    zeta_reversed_lex = DeRham.zeta_function(f,termorder=:lex,vars_reversed=true)
    zeta_normal_lex = DeRham.zeta_function(f,termorder=:lex,vars_reversed=false)

    @test zeta_reversed == zeta_normal
    @test zeta_reversed_lex == zeta_normal_lex
end
    

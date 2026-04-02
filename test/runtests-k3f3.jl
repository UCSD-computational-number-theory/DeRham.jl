using Test
using Oscar
using DeRham

include("ManageCSVTests.jl")

@testset "K3 surfaces over F_3, S=[0,1,2]" begin

    @testset "CPU depthfirst, fastevaluation" begin
        zf = f -> DeRham.zeta_coefficients(f, S=[0,1,2], algorithm=:depthfirst, fastevaluation=true, changef=false)
        runcsvtest("dim2_deg_4_k3f3.csv", zeta_coefficients=zf)
    end

end

# Quality gates wired into `Pkg.test` — the real, asserting gate.
# Aqua slice only; JET / ExplicitImports / DispatchDoctor / AllocCheck are
# tracked separately and land here later.
using Test
using Aqua

@testset "quality gates — DeRham" begin
    M = DeRham

    @testset "Aqua — quality assurance" begin
        Aqua.test_all(M)
    end
end

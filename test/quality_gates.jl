# Quality gates wired into `Pkg.test` — the real, asserting gate.
# Aqua + ExplicitImports slices; JET / DispatchDoctor / AllocCheck are tracked
# separately and land here later.
using Test
using Aqua
using ExplicitImports

@testset "quality gates — DeRham" begin
    M = DeRham

    @testset "Aqua — quality assurance" begin
        Aqua.test_all(M)
    end

    @testset "ExplicitImports" begin
        # No umbrella check() exists — call the seven granular checks. Each
        # returns `nothing` on success and throws an informative error on
        # failure. None deleted, none passed `false`.
        @test ExplicitImports.check_no_implicit_imports(M) === nothing
        @test ExplicitImports.check_no_stale_explicit_imports(M) === nothing
        @test ExplicitImports.check_all_explicit_imports_via_owners(M) === nothing
        @test ExplicitImports.check_no_self_qualified_accesses(M) === nothing
        @test ExplicitImports.check_all_explicit_imports_are_public(M) === nothing

        # @allowscalar (owned by GPUArraysCore, reached via CUDA) and libflint
        # (owned by FLINT_jll, reached via Nemo) are not direct deps of
        # DeRham; adding a `_jll` or GPUArraysCore dependency just to satisfy
        # this lint would be worse than the lint.
        @test ExplicitImports.check_all_qualified_accesses_via_owners(
            M;
            ignore = (Symbol("@allowscalar"), :libflint),
        ) === nothing

        # Scoped, temporary allowlist for non-public qualified accesses. The
        # GPUFiniteFieldMatrices names below need to be made public or exported
        # upstream; drop them from here once that lands.
        # This ExplicitImports version's `ignore` for this check only accepts
        # bare symbols or submodules (no `symbol => module` pairs), so a name
        # here is ignored regardless of which module it's accessed from.
        @test ExplicitImports.check_all_qualified_accesses_are_public(
            M;
            ignore = (
                # GPUFiniteFieldMatrices co-developed sibling internals.
                :CuModcopy,
                :KMatMul!,
                :KaratsubaZeros,
                :Karatsubacopy,
                :TILE_WIDTH,
                :add!,
                :divide_elements!,
                :initialize_plan!,
                :karatsuba_add_helper,
                :karatsuba_scalar_multiply_helper,
                :scalar_multiply!,
                :sub!,
                :zero!,
                :zeros,
                # Deep Nemo/FLINT internals reached through Oscar.
                :fmpz_mod_ctx_struct,
                Symbol("@allowscalar"),
                :libflint,
                # Base.RefValue is widely used despite not being declared public.
                :RefValue,
                # OhMyThreads re-exports TaskLocalValue from its TaskLocalValues
                # dependency without declaring it public itself.
                :TaskLocalValue,
                # CUDA.@time is a debug/profiling call already present in
                # src/; removing it would be a logic change out of scope here.
                Symbol("@time"),
            ),
        ) === nothing
    end
end

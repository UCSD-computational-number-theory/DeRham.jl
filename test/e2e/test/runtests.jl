# test/e2e/test/runtests.jl — ordinary automated tests for the E2E harness
# itself (selection, capability handling, and validation logic). These run
# as part of the conventional test suite, independent of --workflow=/
# --name=/--tag= selection of the E2E catalogue's own cases, and do not
# execute any zeta_coefficients computation.

using Test

@testset "E2E harness" begin
    @testset "selection" begin
        @test E2E.parse_selection(String[]) == E2E.Selection(:workflow, :ci)
        @test E2E.parse_selection(["--workflow=presubmit"]) ==
              E2E.Selection(:workflow, :presubmit)
        @test E2E.parse_selection(["--workflow=full"]) == E2E.Selection(:workflow, :full)
        @test E2E.parse_selection(["--name=foo"]) == E2E.Selection(:name, "foo")
        @test E2E.parse_selection(["--tag=gpu"]) == E2E.Selection(:tag, "gpu")

        @test_throws ErrorException E2E.parse_selection(["--workflow=nonsense"])
        @test_throws ErrorException E2E.parse_selection(["--bogus=1"])
        @test_throws ErrorException E2E.parse_selection(["--name=a", "--tag=b"])
        @test_throws ErrorException E2E.parse_selection(["--workflow=ci", "--name=a"])

        cfg = E2E.CatalogueConfig(
            5.0,
            90.0,
            Set(["gpu", "multithreading", "uint", "bigint", "karatsuba"]),
        )
        mk(name, runtime; tags = String[]) = E2E.ExecutionRecord(
            name,
            "zeta_coefficients",
            tags,
            runtime,
            "ex",
            nothing,
            Dict{String,Any}(),
            "synthetic",
        )
        recs = [
            mk("fast", 1.0),
            mk("mid", 10.0),
            mk("slow", 200.0),
            mk("gpuslow", 200.0; tags = ["gpu"]),
        ]

        presubmit = E2E.select_records(recs, cfg, E2E.Selection(:workflow, :presubmit))
        ci = E2E.select_records(recs, cfg, E2E.Selection(:workflow, :ci))
        full = E2E.select_records(recs, cfg, E2E.Selection(:workflow, :full))
        presubmit_names = Set(r.name for r in presubmit.records)
        ci_names = Set(r.name for r in ci.records)
        full_names = Set(r.name for r in full.records)
        @test presubmit_names ⊆ ci_names
        @test ci_names ⊆ full_names
        @test presubmit_names == Set(["fast"])
        @test ci_names == Set(["fast", "mid"])
        @test full_names == Set(["fast", "mid", "slow", "gpuslow"])
        @test presubmit.broad && ci.broad && full.broad

        name_sel = E2E.select_records(recs, cfg, E2E.Selection(:name, "slow"))
        @test length(name_sel.records) == 1 && name_sel.records[1].name == "slow"
        @test !name_sel.broad
        @test_throws ErrorException E2E.select_records(
            recs,
            cfg,
            E2E.Selection(:name, "nope"),
        )

        tag_sel = E2E.select_records(recs, cfg, E2E.Selection(:tag, "gpu"))
        @test Set(r.name for r in tag_sel.records) == Set(["gpuslow"])
        @test tag_sel.broad
        @test_throws ErrorException E2E.select_records(
            recs,
            cfg,
            E2E.Selection(:tag, "not-a-tag"),
        )
        @test_throws ErrorException E2E.select_records(
            recs,
            cfg,
            E2E.Selection(:tag, "multithreading"),
        )
    end

    @testset "zeta_coefficients executor validation" begin
        exec = E2E.ZetaCoefficientsExecutor

        good_variety(kind = "plane_curve", vars = ["x1", "x2", "x3"]) = Dict{String,Any}(
            "coeff_domain" => Dict("kind" => "integer"),
            "dim" => 1,
            "non_middle_factors" => Dict(
                "kind" => "projective_lefschetz",
                "middle_factor_content" => "full",
            ),
            "model" => Dict(
                "kind" => kind,
                "vars" => vars,
                "monomials" => [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                "coeffs" => [1, 1, 1],
            ),
        )
        example(variety, results) =
            E2E.ExampleCase("ex1", "synthetic.json", variety, results, "")
        record(; p = nothing, params = Dict{String,Any}()) = E2E.ExecutionRecord(
            "rec1",
            "zeta_coefficients",
            String[],
            0.1,
            "ex1",
            p,
            params,
            "synthetic.toml",
        )

        # correct ascending orientation round-trips to DeRham's descending
        # convention.
        ex = example(
            good_variety(),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
        )
        @test exec.expected_descending(ex, record()) == BigInt[7, 1, 1]

        # reversed (descending) coeffs_asc is caught as the specific
        # orientation mistake, not a generic error.
        bad = example(
            good_variety(),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [7, 1, 1]))],
        )
        @test_throws r"orientation mistake" exec.expected_descending(bad, record())

        # unsupported model kind / coeff domain / non_middle_factors convention.
        unsupported_kind = example(
            good_variety("hyperelliptic"),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
        )
        @test_throws ErrorException exec.expected_descending(unsupported_kind, record())

        v = good_variety()
        v["coeff_domain"] = Dict("kind" => "number_field")
        unsupported_domain =
            example(v, [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))])
        @test_throws ErrorException exec.expected_descending(unsupported_domain, record())

        v2 = good_variety()
        v2["non_middle_factors"] = Dict("kind" => "explicit")
        unsupported_nmf = example(v2, [Dict("p" => 7, "l_factors" => Dict())])
        @test_throws ErrorException exec.expected_descending(unsupported_nmf, record())

        # multi-result example requires an explicit p on the record.
        multi = example(
            good_variety(),
            [
                Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7])),
                Dict("p" => 11, "Lpoly" => Dict("coeffs_asc" => [1, 0, 11])),
            ],
        )
        @test_throws ErrorException exec.expected_descending(multi, record())
        @test exec.expected_descending(multi, record(; p = 11)) == BigInt[11, 0, 1]

        # parameter validation.
        @test isnothing(
            exec.validate_parameters(
                "rec",
                Dict("algorithm" => "depthfirst", "S" => [0, 1]),
            ),
        )
        @test_throws ErrorException exec.validate_parameters(
            "rec",
            Dict("algorithm" => "not-an-algorithm"),
        )
        @test_throws ErrorException exec.validate_parameters(
            "rec",
            Dict("not_a_real_param" => 1),
        )

        # capability requirement derivation.
        @test :gpu in exec.required_capabilities(record(; params = Dict("use_gpu" => true)))
        @test isempty(exec.required_capabilities(record()))
        gpu_tagged = E2E.ExecutionRecord(
            "rec2",
            "zeta_coefficients",
            ["gpu"],
            0.1,
            "ex1",
            nothing,
            Dict{String,Any}(),
            "synthetic.toml",
        )
        @test :gpu in exec.required_capabilities(gpu_tagged)
    end

    @testset "catalogue validation aggregates errors" begin
        cfg = E2E.CatalogueConfig(
            5.0,
            90.0,
            Set(["gpu", "multithreading", "uint", "bigint", "karatsuba"]),
        )
        variety = Dict{String,Any}(
            "coeff_domain" => Dict("kind" => "integer"),
            "dim" => 1,
            "non_middle_factors" => Dict(
                "kind" => "projective_lefschetz",
                "middle_factor_content" => "full",
            ),
            "model" => Dict(
                "kind" => "plane_curve",
                "vars" => ["x1", "x2", "x3"],
                "monomials" => [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                "coeffs" => [1, 1, 1],
            ),
        )
        good_example = E2E.ExampleCase(
            "ok",
            "s.json",
            variety,
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
            "",
        )
        catalogue = E2E.ExampleCatalogue(Dict("ok" => good_example))

        records = [
            E2E.ExecutionRecord(
                "r_ok",
                "zeta_coefficients",
                String[],
                0.1,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_missing_ref",
                "zeta_coefficients",
                String[],
                0.1,
                "does_not_exist",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_bad_tag",
                "zeta_coefficients",
                ["not-a-tag"],
                0.1,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_bad_type",
                "not-a-type",
                String[],
                0.1,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_bad_runtime",
                "zeta_coefficients",
                String[],
                -1.0,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
        ]

        err = nothing
        try
            E2E.validate_catalogue(catalogue, records, cfg)
            error("expected CatalogueValidationError")
        catch e
            err = e
        end
        @test err isa E2E.CatalogueValidationError
        @test length(err.errors) == 4
        joined = join(err.errors, "\n")
        @test occursin("r_missing_ref", joined)
        @test occursin("r_bad_tag", joined)
        @test occursin("r_bad_type", joined)
        @test occursin("r_bad_runtime", joined)

        expected = E2E.validate_catalogue(catalogue, [records[1]], cfg)
        @test expected["r_ok"] == BigInt[7, 1, 1]
    end
end

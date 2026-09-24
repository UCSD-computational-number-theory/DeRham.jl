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

    @testset "loader — example documents" begin
        good_variety = Dict{String,Any}(
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
        good_case(id) = Dict{String,Any}(
            "id" => id,
            "variety" => good_variety,
            "results" => [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
            "notes" => "",
        )

        mktempdir() do dir
            bad = joinpath(dir, "bad.json")
            write(bad, "{ not valid json")
            @test_throws r"failed to parse as JSON" E2E._read_example_document(bad)

            notobj = joinpath(dir, "notobj.json")
            write(notobj, "[1, 2, 3]")
            @test_throws r"top level must be a JSON object" E2E._read_example_document(
                notobj,
            )

            wrongver = joinpath(dir, "wrongver.json")
            write(wrongver, """{"schema_version": "2", "cases": []}""")
            @test_throws r"schema_version must be" E2E._read_example_document(wrongver)

            nocases = joinpath(dir, "nocases.json")
            write(nocases, """{"schema_version": "3"}""")
            @test_throws r"missing required key 'cases'" E2E._read_example_document(nocases)

            casesnotarr = joinpath(dir, "casesnotarr.json")
            write(casesnotarr, """{"schema_version": "3", "cases": {}}""")
            @test_throws r"'cases' must be an array" E2E._read_example_document(casesnotarr)

            good = joinpath(dir, "good.json")
            write(good, """{"schema_version": "3", "cases": [{"id": "x"}]}""")
            cases = E2E._read_example_document(good)
            @test length(cases) == 1
            @test cases[1]["id"] == "x"
        end

        @test E2E._parse_example_case(good_case("ex1"), "p.json", 1).id == "ex1"

        missing_key = Dict{String,Any}(
            "id" => "ex1",
            "variety" => good_variety,
            "results" => [Dict("p" => 7)],
        )
        @test_throws r"missing required key 'notes'" E2E._parse_example_case(
            missing_key,
            "p.json",
            1,
        )

        empty_id = merge(good_case("ex1"), Dict{String,Any}("id" => ""))
        @test_throws r"empty or non-string id" E2E._parse_example_case(
            empty_id,
            "p.json",
            1,
        )

        bad_variety_type = merge(good_case("ex1"), Dict{String,Any}("variety" => "nope"))
        @test_throws r"'variety' must be an object" E2E._parse_example_case(
            bad_variety_type,
            "p.json",
            1,
        )

        missing_variety_key = merge(
            good_case("ex1"),
            Dict{String,Any}("variety" => Dict{String,Any}("coeff_domain" => Dict())),
        )
        @test_throws r"variety is missing required key" E2E._parse_example_case(
            missing_variety_key,
            "p.json",
            1,
        )

        empty_results = merge(good_case("ex1"), Dict{String,Any}("results" => []))
        @test_throws r"'results' must be a nonempty array" E2E._parse_example_case(
            empty_results,
            "p.json",
            1,
        )

        non_array_results =
            merge(good_case("ex1"), Dict{String,Any}("results" => "notarray"))
        @test_throws r"'results' must be a nonempty array" E2E._parse_example_case(
            non_array_results,
            "p.json",
            1,
        )

        case_json(id) = """
        {"id": "$(id)", "variety": {"coeff_domain": {"kind": "integer"}, "dim": 1, "non_middle_factors": {"kind": "projective_lefschetz", "middle_factor_content": "full"}, "model": {"kind": "plane_curve", "vars": ["x1", "x2", "x3"], "monomials": [[3, 0, 0], [0, 3, 0], [0, 0, 3]], "coeffs": [1, 1, 1]}}, "results": [{"p": 7, "Lpoly": {"coeffs_asc": [1, 1, 7]}}], "notes": ""}
        """
        document_json(cases_json) = """{"schema_version": "3", "cases": [$(cases_json)]}"""

        mktempdir() do dir
            dup_within = joinpath(dir, "dup_within.json")
            write(dup_within, document_json(case_json("dup") * ", " * case_json("dup")))
            into = Dict{String,E2E.ExampleCase}()
            @test_throws r"duplicate example id" E2E._load_examples_into!(into, dup_within)

            f1 = joinpath(dir, "a.json")
            write(f1, document_json(case_json("shared")))
            f2 = joinpath(dir, "b.json")
            write(f2, document_json(case_json("shared")))
            into2 = Dict{String,E2E.ExampleCase}()
            E2E._load_examples_into!(into2, f1)
            @test_throws r"duplicate example id" E2E._load_examples_into!(into2, f2)
        end
    end

    @testset "loader — execution records" begin
        mktempdir() do dir
            records_dir = joinpath(dir, "records")
            mkpath(records_dir)

            bad = joinpath(records_dir, "bad.toml")
            write(bad, "name = \"unterminated")
            @test_throws r"failed to parse as TOML" E2E.load_execution_records(dir)
            rm(bad)

            missing_key = joinpath(records_dir, "missing_key.toml")
            write(missing_key, "[[tests]]\nname = \"x\"\n")
            @test_throws r"missing required key 'test'" E2E.load_execution_records(dir)
            rm(missing_key)

            not_array = joinpath(records_dir, "not_array.toml")
            write(not_array, "test = \"nope\"\n")
            @test_throws r"'test' must be an array of tables" E2E.load_execution_records(
                dir,
            )
            rm(not_array)

            missing_field = joinpath(records_dir, "missing_field.toml")
            write(missing_field, "[[test]]\nname = \"r1\"\ntype = \"zeta_coefficients\"\n")
            @test_throws r"missing required key 'example'" E2E.load_execution_records(dir)
            rm(missing_field)

            empty_name = joinpath(records_dir, "empty_name.toml")
            write(
                empty_name,
                "[[test]]\nname = \"\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\n",
            )
            @test_throws r"empty or non-string name" E2E.load_execution_records(dir)
            rm(empty_name)

            r1 = joinpath(records_dir, "r1.toml")
            write(
                r1,
                "[[test]]\nname = \"dup\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\n",
            )
            r2 = joinpath(records_dir, "r2.toml")
            write(
                r2,
                "[[test]]\nname = \"dup\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\n",
            )
            @test_throws r"duplicate execution record name" E2E.load_execution_records(dir)
            rm(r1)
            rm(r2)

            ok = joinpath(records_dir, "ok.toml")
            write(
                ok,
                "[[test]]\nname = \"ok1\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\nruntime = 1.0\n",
            )
            recs = E2E.load_execution_records(dir)
            @test length(recs) == 1 && recs[1].name == "ok1"
        end
    end
end

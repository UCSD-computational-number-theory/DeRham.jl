# Run.jl — executes the selected records serially (no concurrency), against
# the already-validated expected values, and reports aggregate estimated
# runtime plus each case's expected and actual runtime. Runtime metadata is
# workflow guidance only: no timing sweep, drift warning, or performance
# assertion is derived from it.

struct CaseOutcome
    name::String
    passed::Bool
    skipped::Bool
    skip_reason::Union{Nothing,String}
    expected_runtime::Float64
    actual_runtime::Float64
end

function run_selected(
    catalogue::ExampleCatalogue,
    selected::NamedTuple,
    expected::Dict{String,Any},
)
    records = selected.records
    broad = selected.broad

    estimated_total = isempty(records) ? 0.0 : sum(r.runtime for r in records)
    println(
        "[e2e] selected $(length(records)) case(s), estimated total runtime $(round(estimated_total, digits = 3))s",
    )

    outcomes = CaseOutcome[]
    for record in records
        executor = EXECUTOR_REGISTRY[record.type]
        caps = executor.required_capabilities(record)
        unavailable = [c for c in caps if !executor.capability_available(c)]

        if !isempty(unavailable)
            reason = "requires unavailable capabilities: $(sort(collect(unavailable)))"
            if broad
                println("[e2e] SKIP $(record.name): $(reason)")
                push!(
                    outcomes,
                    CaseOutcome(record.name, false, true, reason, record.runtime, 0.0),
                )
                continue
            else
                error("execution record $(repr(record.name)) $(reason)")
            end
        end

        example = catalogue.examples[record.example]
        F = executor.build_input(example, record)
        actual_runtime = @elapsed raw = executor.execute(F, record.parameters)
        actual = executor.normalize(raw)
        passed = executor.compare(expected[record.name], actual)
        status = passed ? "PASS" : "FAIL"
        println(
            "[e2e] $(status) $(record.name)  expected_runtime=$(record.runtime)s actual_runtime=$(round(actual_runtime, digits = 3))s",
        )
        @test passed
        push!(
            outcomes,
            CaseOutcome(
                record.name,
                passed,
                false,
                nothing,
                record.runtime,
                actual_runtime,
            ),
        )
    end
    return outcomes
end

function run_e2e(e2e_root::AbstractString, args::AbstractVector{<:AbstractString})
    config = load_config(joinpath(e2e_root, "config.toml"))
    catalogue = load_example_catalogue(e2e_root)
    records = load_execution_records(e2e_root)
    expected = validate_catalogue(catalogue, records, config)
    selection = parse_selection(args)
    selected = select_records(records, config, selection)
    return run_selected(catalogue, selected, expected)
end

# Validate.jl — validates the entire combined example + execution-record
# catalogue before any selected case runs. Cheap: no Oscar objects are built
# and no zeta computation happens here (see executors/ZetaCoefficients.jl).
# Errors are aggregated across the whole catalogue and reported together,
# each tagged with its source file and record/example id, rather than
# failing on the first problem found.

struct CatalogueValidationError <: Exception
    errors::Vector{String}
end

function Base.showerror(io::IO, e::CatalogueValidationError)
    println(io, "E2E catalogue validation failed with $(length(e.errors)) error(s):")
    for msg in e.errors
        println(io, "  - $(msg)")
    end
end

# Validates `catalogue` + `records` against `config` and the executor
# registry. Returns a Dict{String,Any} mapping each execution record's name
# to its resolved expected value (computed once here, reused at run time),
# or throws a CatalogueValidationError listing every problem found.
function validate_catalogue(
    catalogue::ExampleCatalogue,
    records::Vector{ExecutionRecord},
    config::CatalogueConfig,
)
    errors = String[]
    expected = Dict{String,Any}()

    for record in records
        loc = "$(record.source_path) (record $(repr(record.name)))"

        if !isfinite(record.runtime) || record.runtime < 0
            push!(errors, "$(loc): runtime must be a finite, nonnegative number, got $(record.runtime)")
        end

        for tag in record.tags
            tag in config.allowed_tags ||
                push!(errors, "$(loc): unknown tag $(repr(tag)); allowed tags are $(sort(collect(config.allowed_tags)))")
        end

        executor = get(EXECUTOR_REGISTRY, record.type, nothing)
        if executor === nothing
            push!(errors, "$(loc): unknown type $(repr(record.type)); known types are $(sort(collect(keys(EXECUTOR_REGISTRY))))")
            continue
        end

        example = get(catalogue.examples, record.example, nothing)
        if example === nothing
            push!(errors, "$(loc): references unknown example $(repr(record.example))")
            continue
        end

        try
            executor.validate_parameters(record.name, record.parameters)
        catch e
            push!(errors, "$(loc): $(sprint(showerror, e))")
        end

        try
            expected[record.name] = executor.expected_descending(example, record)
        catch e
            push!(errors, "$(loc): $(sprint(showerror, e))")
        end
    end

    isempty(errors) || throw(CatalogueValidationError(errors))
    return expected
end

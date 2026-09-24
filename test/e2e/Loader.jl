# Loader.jl — resolve the pinned zeta_test_suite artifact and load examples
# (artifact + local, same schema) and TOML execution records into one
# catalogue. Normal `using DeRham` never triggers a download: the artifact
# is lazy and is only resolved when the E2E suite actually runs.

const SCHEMA_VERSION = "3"
const ARTIFACT_NAME = "zeta_test_suite"
const ARTIFACT_CASES_SUBDIR = "zeta_test_suite-3.1.0"

function artifact_examples_dir(e2e_root::AbstractString)
    artifacts_toml = joinpath(dirname(e2e_root), "Artifacts.toml")
    isfile(artifacts_toml) || error("Artifacts.toml not found at $(artifacts_toml)")
    path = ensure_artifact_installed(ARTIFACT_NAME, artifacts_toml)
    casesdir = joinpath(path, ARTIFACT_CASES_SUBDIR, "cases")
    isdir(casesdir) ||
        error("resolved artifact $(ARTIFACT_NAME) is missing expected cases/ at $(casesdir)")
    return casesdir
end

function _find_files(dir::AbstractString, ext::AbstractString)
    out = String[]
    isdir(dir) || return out
    for (root, _, names) in walkdir(dir)
        for name in names
            if endswith(name, ext)
                push!(out, joinpath(root, name))
            end
        end
    end
    return sort(out)
end

# Validates the required top-level shape of one v3 document and returns its
# `cases` array; throws with the offending file path on any malformed
# document.
function _read_example_document(path::AbstractString)
    local data
    try
        data = JSON.parsefile(path)
    catch e
        error("$(path): failed to parse as JSON ($(e))")
    end
    data isa AbstractDict || error("$(path): top level must be a JSON object")
    get(data, "schema_version", nothing) == SCHEMA_VERSION ||
        error(
            "$(path): schema_version must be $(repr(SCHEMA_VERSION)), got $(repr(get(data, "schema_version", nothing)))",
        )
    haskey(data, "cases") || error("$(path): missing required key 'cases'")
    cases = data["cases"]
    cases isa AbstractVector || error("$(path): 'cases' must be an array")
    return cases
end

const _REQUIRED_CASE_KEYS = ("id", "variety", "results", "notes")
const _REQUIRED_VARIETY_KEYS = ("coeff_domain", "dim", "non_middle_factors", "model")

function _parse_example_case(case::AbstractDict, path::AbstractString, index::Int)
    for k in _REQUIRED_CASE_KEYS
        haskey(case, k) || error("$(path): case #$(index) is missing required key '$(k)'")
    end
    id = case["id"]
    (id isa AbstractString && !isempty(id)) ||
        error("$(path): case #$(index) has an empty or non-string id")
    variety = case["variety"]
    variety isa AbstractDict ||
        error("$(path): case $(repr(id)): 'variety' must be an object")
    for k in _REQUIRED_VARIETY_KEYS
        haskey(variety, k) ||
            error("$(path): case $(repr(id)): variety is missing required key '$(k)'")
    end
    results = case["results"]
    results isa AbstractVector && !isempty(results) ||
        error("$(path): case $(repr(id)): 'results' must be a nonempty array")
    notes = get(case, "notes", "")
    return ExampleCase(id, path, variety, results, string(notes))
end

# Loads every example from `path`, tracking source locations for
# diagnostics, and merges them into `into`. Raises on a duplicate id, either
# within this source or against what is already in `into` (covers
# collisions both within and across the artifact / local sources).
function _load_examples_into!(
    into::Dict{String,ExampleCase},
    path::AbstractString,
)
    for (index, case) in enumerate(_read_example_document(path))
        case isa AbstractDict ||
            error("$(path): case #$(index) must be a JSON object")
        example = _parse_example_case(case, path, index)
        if haskey(into, example.id)
            existing = into[example.id]
            error(
                "duplicate example id $(repr(example.id)): defined in both $(existing.source_path) and $(path)",
            )
        end
        into[example.id] = example
    end
    return into
end

function load_example_catalogue(e2e_root::AbstractString)
    examples = Dict{String,ExampleCase}()

    artifact_dir = artifact_examples_dir(e2e_root)
    for path in _find_files(artifact_dir, ".json")
        _load_examples_into!(examples, path)
    end

    local_dir = joinpath(e2e_root, "examples")
    for path in _find_files(local_dir, ".json")
        _load_examples_into!(examples, path)
    end

    return ExampleCatalogue(examples)
end

# ----------------------------------------------------------------------
# TOML execution records
# ----------------------------------------------------------------------

const _REQUIRED_RECORD_KEYS = ("name", "type", "example")

function _parse_execution_record(raw::AbstractDict, path::AbstractString)
    for k in _REQUIRED_RECORD_KEYS
        haskey(raw, k) || error("$(path): execution record is missing required key '$(k)'")
    end
    name = raw["name"]
    (name isa AbstractString && !isempty(name)) ||
        error("$(path): execution record has an empty or non-string name")
    type = String(raw["type"])
    tags = String.(get(raw, "tags", String[]))
    runtime = Float64(get(raw, "runtime", 0.0))
    example = String(raw["example"])
    p = haskey(raw, "p") ? Int(raw["p"]) : nothing
    parameters = Dict{String,Any}(get(raw, "parameters", Dict{String,Any}()))
    return ExecutionRecord(name, type, tags, runtime, example, p, parameters, path)
end

function load_execution_records(e2e_root::AbstractString)
    records = ExecutionRecord[]
    by_name = Dict{String,String}()
    records_dir = joinpath(e2e_root, "records")
    for path in _find_files(records_dir, ".toml")
        local raw
        try
            raw = TOML.parsefile(path)
        catch e
            error("$(path): failed to parse as TOML ($(e))")
        end
        haskey(raw, "test") || continue
        entries = raw["test"]
        entries isa AbstractVector ||
            error("$(path): 'test' must be an array of tables ([[test]])")
        for entry in entries
            record = _parse_execution_record(entry, path)
            if haskey(by_name, record.name)
                error(
                    "duplicate execution record name $(repr(record.name)): defined in both $(by_name[record.name]) and $(path)",
                )
            end
            by_name[record.name] = path
            push!(records, record)
        end
    end
    return records
end

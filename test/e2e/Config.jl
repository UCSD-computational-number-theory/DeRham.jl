# Config.jl — config.toml model: workflow runtime thresholds and the
# allowed-tag set. Tags have no implied relationships; directory placement
# has no execution semantics.

struct CatalogueConfig
    presubmit_seconds::Float64
    ci_seconds::Float64
    allowed_tags::Set{String}
end

function load_config(path::AbstractString)
    isfile(path) || error("E2E config not found: $(path)")
    raw = TOML.parsefile(path)

    haskey(raw, "thresholds") || error("$(path): missing [thresholds] table")
    thresholds = raw["thresholds"]
    haskey(thresholds, "presubmit_seconds") ||
        error("$(path): thresholds.presubmit_seconds is required")
    haskey(thresholds, "ci_seconds") || error("$(path): thresholds.ci_seconds is required")
    presubmit_seconds = Float64(thresholds["presubmit_seconds"])
    ci_seconds = Float64(thresholds["ci_seconds"])
    presubmit_seconds > 0 || error("$(path): thresholds.presubmit_seconds must be > 0")
    ci_seconds > presubmit_seconds ||
        error("$(path): thresholds.ci_seconds must be > thresholds.presubmit_seconds")

    haskey(raw, "tags") || error("$(path): missing [tags] table")
    haskey(raw["tags"], "allowed") || error("$(path): tags.allowed is required")
    allowed_tags = Set{String}(String.(raw["tags"]["allowed"]))

    return CatalogueConfig(presubmit_seconds, ci_seconds, allowed_tags)
end

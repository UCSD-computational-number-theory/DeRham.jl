# Select.jl — test_args selection: exactly one public mode per invocation
# (--workflow=, --name=, or --tag=); bare Pkg.test() defaults to ci.
#
# Workflow membership uses each record's recorded runtime against the
# configured strict thresholds; presubmit ⊆ ci ⊆ full by construction, since
# each is `runtime < threshold` over the same threshold ordering (full has
# no threshold). Exact name and tag selections ignore runtime thresholds.

struct Selection
    mode::Symbol
    value::Any
end

function parse_selection(args::AbstractVector{<:AbstractString})
    modes = Selection[]
    for arg in args
        if startswith(arg, "--workflow=")
            raw = arg[(length("--workflow=") + 1):end]
            wf = Symbol(raw)
            wf in (:presubmit, :ci, :full) ||
                error("unknown workflow $(repr(raw)); expected presubmit, ci, or full")
            push!(modes, Selection(:workflow, wf))
        elseif startswith(arg, "--name=")
            push!(modes, Selection(:name, arg[(length("--name=") + 1):end]))
        elseif startswith(arg, "--tag=")
            push!(modes, Selection(:tag, arg[(length("--tag=") + 1):end]))
        else
            error(
                "unknown E2E test_args selector $(repr(arg)); expected --workflow=, --name=, or --tag=",
            )
        end
    end

    isempty(modes) && return Selection(:workflow, :ci)
    length(modes) == 1 || error(
        "conflicting E2E selectors ($(join(["$(m.mode)=$(m.value)" for m in modes], ", "))); " *
        "pass exactly one of --workflow=, --name=, or --tag=",
    )
    return only(modes)
end

# Returns (broad, records): `broad` selections (workflow, tag) skip cases
# whose capability is unavailable; `name` (exact) selections fail hard
# instead.
function select_records(
    records::Vector{ExecutionRecord},
    config::CatalogueConfig,
    selection::Selection,
)
    if selection.mode == :workflow
        threshold = if selection.value == :presubmit
            config.presubmit_seconds
        elseif selection.value == :ci
            config.ci_seconds
        else
            Inf
        end
        return (broad = true, records = filter(r -> r.runtime < threshold, records))
    elseif selection.mode == :name
        matches = filter(r -> r.name == selection.value, records)
        isempty(matches) && error("no execution record named $(repr(selection.value))")
        return (broad = false, records = matches)
    elseif selection.mode == :tag
        selection.value in config.allowed_tags || error(
            "unknown tag $(repr(selection.value)); allowed tags are $(sort(collect(config.allowed_tags)))",
        )
        matches = filter(r -> selection.value in r.tags, records)
        isempty(matches) && error("no execution record has tag $(repr(selection.value))")
        return (broad = true, records = matches)
    end
    error("unreachable selection mode $(selection.mode)")
end

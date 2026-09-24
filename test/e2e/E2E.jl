# E2E.jl — the declarative E2E catalogue. See test/e2e/README.md for the
# contributor workflow (adding examples, execution records, running by
# name/tag/workflow).

module E2E

using Test
using TOML
using JSON
using Pkg.Artifacts
import Main.DeRham

include("Model.jl")
include("Config.jl")
include("Loader.jl")
include("executors/ZetaCoefficients.jl")

const EXECUTOR_REGISTRY =
    Dict{String,Module}("zeta_coefficients" => ZetaCoefficientsExecutor)

include("Validate.jl")
include("Select.jl")
include("Run.jl")

# Entry point called from test/runtests.jl with Base.ARGS (i.e. Pkg.test's
# test_args). Bare Pkg.test() passes no args, which defaults to --workflow=ci.
function run!(args::AbstractVector{<:AbstractString} = ARGS)
    return run_e2e(@__DIR__, args)
end

end # module E2E

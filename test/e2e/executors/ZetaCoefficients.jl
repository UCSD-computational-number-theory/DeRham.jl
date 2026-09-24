# executors/ZetaCoefficients.jl — the sole entry in the executor registry.
# Owns resolved-example validation, parameter validation, capability checks,
# execution, output normalization, and exact comparison for
# type = "zeta_coefficients".
#
# Supports exactly the subset of the zeta_test_suite v3 schema the migrated
# and locally-authored examples use: coeff_domain.kind == "integer", model
# kind plane_curve or projective_hypersurface, and non_middle_factors
# {kind = "projective_lefschetz", middle_factor_content = "full"} — the
# convention under which DeRham.zeta_coefficients' return value (the full
# middle-cohomology L-factor, descending / leading-coefficient-first) is
# directly comparable to a result's Lpoly.

module ZetaCoefficientsExecutor

using Oscar
using CUDA: CUDA
import ..DeRham
import ..ExampleCase
import ..ExecutionRecord

const SUPPORTED_MODEL_KINDS = Set(["plane_curve", "projective_hypersurface"])

const KNOWN_PARAMETERS = Set([
    "S",
    "verbose",
    "changef",
    "givefrobmat",
    "algorithm",
    "termorder",
    "vars_reversed",
    "fastevaluation",
    "always_use_bigints",
    "use_gpu",
    "use_threads",
])
const KNOWN_ALGORITHMS = Set(["default", "depthfirst", "varbyvar", "pchunk"])

_to_bigint(c::Integer) = BigInt(c)
_to_bigint(c::AbstractString) = parse(BigInt, c)

function validate_parameters(name::AbstractString, params::AbstractDict)
    for k in keys(params)
        k in KNOWN_PARAMETERS ||
            error("execution record $(repr(name)): unknown parameter $(repr(k))")
    end
    if haskey(params, "algorithm")
        params["algorithm"] isa AbstractString && params["algorithm"] in KNOWN_ALGORITHMS ||
            error(
                "execution record $(repr(name)): unknown algorithm $(repr(get(params, "algorithm", nothing)))",
            )
    end
    if haskey(params, "termorder") && !(params["termorder"] isa AbstractString)
        error("execution record $(repr(name)): termorder must be a string")
    end
    if haskey(params, "S") && !(params["S"] isa AbstractVector)
        error("execution record $(repr(name)): S must be an array of integers")
    end
    return nothing
end

function required_capabilities(record::ExecutionRecord)
    caps = Set{Symbol}()
    if "gpu" in record.tags || get(record.parameters, "use_gpu", false) == true
        push!(caps, :gpu)
    end
    return caps
end

capability_available(::Val{:gpu}) = CUDA.functional()
capability_available(cap::Symbol) = capability_available(Val(cap))

function _find_result(example::ExampleCase, p::Integer)
    matches = [r for r in example.results if _to_bigint(r["p"]) == BigInt(p)]
    isempty(matches) && error("example $(repr(example.id)) has no result for p=$(p)")
    length(matches) > 1 &&
        error("example $(repr(example.id)) has more than one result for p=$(p)")
    return matches[1]
end

function resolve_prime(example::ExampleCase, record::ExecutionRecord)
    if record.p !== nothing
        return record.p
    end
    length(example.results) == 1 || error(
        "execution record $(repr(record.name)) references example $(repr(example.id)), which has $(length(example.results)) results; 'p' is required to disambiguate",
    )
    return Int(_to_bigint(example.results[1]["p"]))
end

# Validates the (example, record) reference cheaply — no Oscar objects are
# built and no mathematics is executed — and returns the expected value in
# DeRham's own descending (leading-coefficient-first) convention.
function expected_descending(example::ExampleCase, record::ExecutionRecord)
    variety = example.variety
    kind = variety["model"]["kind"]
    kind in SUPPORTED_MODEL_KINDS || error(
        "example $(repr(example.id)): zeta_coefficients does not support model.kind $(repr(kind))",
    )
    get(variety["coeff_domain"], "kind", nothing) == "integer" || error(
        "example $(repr(example.id)): zeta_coefficients only supports coeff_domain.kind == \"integer\"",
    )
    nmf = variety["non_middle_factors"]
    (
        get(nmf, "kind", nothing) == "projective_lefschetz" &&
        get(nmf, "middle_factor_content", nothing) == "full"
    ) || error(
        "example $(repr(example.id)): zeta_coefficients only supports non_middle_factors " *
        "{kind = \"projective_lefschetz\", middle_factor_content = \"full\"}, got $(nmf)",
    )

    p = resolve_prime(example, record)
    result = _find_result(example, p)
    haskey(result, "Lpoly") || error(
        "example $(repr(example.id)), p=$(p): zeta_coefficients requires a Lpoly result, got keys $(collect(keys(result)))",
    )
    coeffs_asc = _to_bigint.(result["Lpoly"]["coeffs_asc"])
    isempty(coeffs_asc) &&
        error("example $(repr(example.id)), p=$(p): Lpoly.coeffs_asc must be nonempty")
    if coeffs_asc[1] != 1
        if coeffs_asc[end] == 1
            error(
                "example $(repr(example.id)), p=$(p): coefficient-orientation mistake — " *
                "Lpoly.coeffs_asc looks reversed (descending); the schema requires ascending order with coeffs_asc[1] == 1",
            )
        else
            error(
                "example $(repr(example.id)), p=$(p): invalid Lpoly — coeffs_asc[1] must equal 1, got $(coeffs_asc[1])",
            )
        end
    end
    return reverse(coeffs_asc)
end

function _multivariable_poly(monomials, coeffs, R)
    gens_R = Oscar.gens(R)
    nv = length(gens_R)
    out = zero(R)
    for (mono, c) in zip(monomials, coeffs)
        length(mono) == nv ||
            error("monomial length $(length(mono)) does not match $(nv) variables")
        term = R(_to_bigint(c))
        for (j, e) in enumerate(mono)
            term *= gens_R[j]^Int(e)
        end
        out += term
    end
    return out
end

# Builds the Oscar polynomial for `example` over GF(p). This is the one part
# of resolution that constructs real Oscar objects, so it is deferred to
# execution time for selected cases only, never run during full-catalogue
# validation.
function build_input(example::ExampleCase, record::ExecutionRecord)
    p = resolve_prime(example, record)
    model = example.variety["model"]
    R, _ = Oscar.polynomial_ring(Oscar.GF(Int(p)), Symbol.(model["vars"]))
    F = _multivariable_poly(model["monomials"], model["coeffs"], R)
    return F
end

function parse_kwargs(params::AbstractDict)
    kwargs = Dict{Symbol,Any}()
    for (k, v) in params
        sym = Symbol(k)
        kwargs[sym] = if sym in (:algorithm, :termorder)
            Symbol(v)
        elseif sym == :S
            Int.(v)
        else
            v
        end
    end
    return kwargs
end

function execute(F, params::AbstractDict)
    kwargs = parse_kwargs(params)
    return DeRham.zeta_coefficients(F; kwargs...)
end

normalize(raw) = BigInt.(raw)

compare(expected::Vector{BigInt}, actual::Vector{BigInt}) = expected == actual

end # module ZetaCoefficientsExecutor

# E2E test catalogue

A declarative, metadata-driven system for DeRham's end-to-end
`zeta_coefficients` tests. Example data (the mathematical variety plus its
known-correct zeta function) is separate from execution records (which
example to run, with which parameters, tagged and timed). `Pkg.test()` runs
it as part of the conventional test suite; see `test/runtests.jl`.

## Add local JSON example data

Add a file under `test/e2e/examples/` (any subdirectory — placement has no
execution semantics) containing one or more grouped cases in the
[zeta_test_suite v3 schema](https://github.com/jjgarzella/zeta_test_suite)
(`schema.json` there is the normative definition):

```json
{
  "schema_version": "3",
  "cases": [
    {
      "id": "my_new_example",
      "variety": {
        "coeff_domain": { "kind": "integer" },
        "dim": 1,
        "non_middle_factors": {
          "kind": "projective_lefschetz",
          "middle_factor_content": "full"
        },
        "model": {
          "kind": "plane_curve",
          "pretty": "x1^3 + x2^3 + x3^3 = 0",
          "vars": ["x1", "x2", "x3"],
          "monomials": [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
          "coeffs": [1, 1, 1]
        }
      },
      "results": [
        { "p": 7, "Lpoly": { "coeffs_asc": [1, 1, 7] } }
      ],
      "notes": "why this case exists"
    }
  ]
}
```

The `zeta_coefficients` executor supports `model.kind` of `plane_curve` (3
variables) or `projective_hypersurface` (>= 2 variables), `coeff_domain.kind
== "integer"`, and `non_middle_factors = {kind: "projective_lefschetz",
middle_factor_content: "full"}`. `Lpoly.coeffs_asc` is **ascending** (constant
term first); DeRham's own `zeta_coefficients` return value is **descending**
(leading coefficient first) — the executor reverses one to compare against
the other, and rejects an accidentally-reversed `coeffs_asc` (its first entry
must be `1`).

Example ids must be nonempty and globally unique across every local file and
the pinned artifact.

## Add a TOML execution record

Add (or append to) a `[[test]]` array-of-tables entry under
`test/e2e/records/` (any subdirectory):

```toml
[[test]]
name = "my_new_example_default"
type = "zeta_coefficients"
tags = []
runtime = 0.1
example = "my_new_example"
parameters = { algorithm = "depthfirst", fastevaluation = true }
```

- `name` must be nonempty and globally unique across every record file.
- `runtime` is an approximate, bucket-oriented estimate in seconds — not a
  measured average — used only for `--workflow=` membership.
- `example` must match an example `id`. If that example has more than one
  prime result, also set `p = <prime>` to disambiguate.
- `tags` must be a subset of `config.toml`'s `[tags] allowed` list (currently
  `gpu`, `multithreading`, `uint`, `bigint`, `karatsuba`). A `gpu` tag (or
  `parameters.use_gpu = true`) marks the record as requiring a functional
  CUDA device.
- `parameters` are passed straight through as keyword arguments to
  `DeRham.zeta_coefficients`.

## Reference an example from the pinned artifact

Artifact examples (from the `zeta_test_suite` Julia artifact, resolved lazily
— `using DeRham` alone never downloads it) load into the same catalogue as
local examples, keyed by their own `id`. Reference one from a TOML record
exactly as you would a local example; there is no separate syntax.

## Run tests by name, tag, or workflow

`Pkg.test` passes `test_args` straight through to `ARGS`. Exactly one
selector may be given; a bare `Pkg.test()` defaults to `--workflow=ci`.

```julia
Pkg.test("DeRham")                                    # --workflow=ci (default)
Pkg.test("DeRham"; test_args = ["--workflow=presubmit"])
Pkg.test("DeRham"; test_args = ["--workflow=full"])
Pkg.test("DeRham"; test_args = ["--name=my_new_example_default"])
Pkg.test("DeRham"; test_args = ["--tag=bigint"])
```

`--workflow=` selects by each record's `runtime` against the configured
`presubmit_seconds` / `ci_seconds` thresholds (`presubmit ⊆ ci ⊆ full`).
`--name=` and `--tag=` ignore runtime thresholds; a `--tag=` selection (like
a workflow) silently skips cases whose capability (e.g. `gpu`) is
unavailable, while an unavailable `--name=` selection fails.

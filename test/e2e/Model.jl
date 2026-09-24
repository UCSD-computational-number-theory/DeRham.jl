# Model.jl — the catalogue data model.
#
# ExampleCase mirrors one grouped test case of the zeta_test_suite v3 schema
# (schema.json / spec.txt, pinned unmodified — see loader.jl). ExecutionRecord
# is a TOML execution record: a nonempty globally unique name, fixed tags, a
# finite nonnegative runtime, a type, a reference to one example (and,
# if the example carries more than one prime result, which prime), and
# optional parameters.

struct ExampleCase
    id::String
    source_path::String
    variety::Dict{String,Any}
    results::Vector{Any}
    notes::String
end

struct ExampleCatalogue
    examples::Dict{String,ExampleCase}
end

struct ExecutionRecord
    name::String
    type::String
    tags::Vector{String}
    runtime::Float64
    example::String
    p::Union{Nothing,Int}
    parameters::Dict{String,Any}
    source_path::String
end

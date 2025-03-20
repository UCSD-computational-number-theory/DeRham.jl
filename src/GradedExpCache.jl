
abstract type GradedExpCache end

#TODO: convert to SVector, that might help with allocation pressure
#using StaticArrays

struct PolyExpCache <: GradedExpCache
    n::Int
    termorder::Symbol
    exp_vec_pieces::Vector{Union{Vector{Vector{Int64}},Nothing}}
    rev_look_pieces::Vector{Union{Dict{Vector{Int64},Int},Nothing}}

    PolyExpCache(n,termorder) = new(n,termorder,[],[])
end


"""
Gets the exponent vector cache for the degree d,
returning it as a vector.

This throws an error if 
"""
function Base.getindex(c::PolyExpCache,d::Int)
    if exp_vec_pieces[d] == nothing
        errormsg = "This PolyExpCache does not have degree=$d stored. "* 
                   "This probably means that the PolyExpCache hasn't been "* 
                   "fully initialized, or there is a bug. If you want to "* 
                   "initialize this cache to work for degree $d, use the "*
                   "`get_forward` function."
        throw(ArgumentError(errormsg))
    end
    exp_vec_pieces[d]
end

function Base.show(io::IO, c::PolyExpCache)
    msg = """PolyExpCache for the graded ring k[x_1, ..., x_$(c.n)]
    Term order: $(c.termorder)
    Terms cached for degrees: $(cached_degrees(c))
    Reverse lookups cached for degrees: $(cached_reverse_lookups(c))
    """
    print(io, msg)
end

cached_degrees(c::PolyExpCache) = findall(c.exp_vec_pieces .!= nothing)
cached_reverse_lookups(c::PolyExpCache) = findall(c.rev_look_pieces .!= nothing)

function get_forward(c::PolyExpCache,d)
    if c.exp_vec_pieces[d] == nothing
        generate_degree_forward(c,d)
    end
    c.exp_vec_pieces[d]
end

function generate_degree_forward(c::PolyExpCache,d)
    l = length(c.exp_vec_pieces) 
    if l < d
        append!(c.exp_vec_pieces,fill(nothing,d - l))
    end
    c.exp_vec_pieces[d] = gen_exp_vec(c.n,d,c.termorder)
end

function generate_degree_reverse(c::PolyExpCache,d)
    evs = get_forward(c,d)

    l = length(c.rev_look_pieces) 
    if l < d
        append!(c.rev_look_pieces,fill(nothing,d - l))
    end
    c.rev_look_pieces[d] = Dict(evs[i] => i for i in 1:length(evs))
end

function PolyExpCache(n,termorder,degrees_to_prefill,reverse_to_prefill)
    c = PolyExpCache(n,termorder)
    for d in degrees_to_prefill
        generate_degree_forward(c,d)
    end

    for d in reverse_to_prefill
        generate_degree_reverse(c,d)
    end
    c
end

function PolyExpCache(n,termorder,max_degree_prefill)
    PolyExpCache(n,termorder,1:max_degree_prefill,1:max_degree_prefill)
end



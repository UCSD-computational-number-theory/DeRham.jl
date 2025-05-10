
abstract type GradedExpCache end

#TODO: convert to SVector, that might help with allocation pressure
#using StaticArrays

"""
A PolyExpCache is a GradedExpCache for the stadard 
\\Z-grading on the polynomial ring.
"""
struct PolyExpCache <: GradedExpCache
    n::Int
    termorder::Symbol
    vars_reversed::Bool #TODO: remove this and use the termorder instead
    exp_vec_pieces::Vector{Union{Vector{Vector{Int64}},Nothing}}
    rev_look_pieces::Vector{Union{Dict{Vector{Int64},Int},Nothing}}
    constant_term::Vector{Vector{Int}}

    function PolyExpCache(n,termorder;vars_reversed=false) 
        new(n,termorder,vars_reversed,[],[],[zeros(Int,n)])
    end
end


"""
Gets the exponent vector cache for the degree d,
returning it as a vector.

This throws an error if the cached degree doesn't exist.
"""
function Base.getindex(c::PolyExpCache,d::Int)
    d < 0 && return []
    d == 0 && return c.constant_term

    if c.exp_vec_pieces[d] == nothing
        errormsg = "This PolyExpCache does not have degree=$d stored. "* 
                   "This probably means that the PolyExpCache hasn't been "* 
                   "fully initialized, or there is a bug. If you want to "* 
                   "initialize this cache to work for degree $d, use the "*
                   "`get_forward` function."
        throw(ArgumentError(errormsg))
    end
    c.exp_vec_pieces[d]
end

function Base.getindex(c::PolyExpCache,d::Int,kind::Symbol)
    d < 0 && return []
    d == 0 && return c.constant_term

    if kind == :forward
        c[d]
    elseif kind == :reverse

        if c.rev_look_pieces[d] == nothing
            errormsg = "This PolyExpCache does not have degree=$d stored for "*
                       "reverse lookup. "* 
                       "This probably means that the PolyExpCache hasn't been "* 
                       "fully initialized, or there is a bug. If you want to "* 
                       "initialize this cache to work for degree $d, use the "*
                       "`get_forward` and `generate_degree_reverse` functions."
            throw(ArgumentError(errormsg))
        end
        c.rev_look_pieces[d]
    else
        throw(ArgumentError("Invalid option for indexing PolyExpCache"))
    end
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
    d ≤ 0 && return 
    d in cached_degrees(c) && return 

    l = length(c.exp_vec_pieces) 
    if l < d
        append!(c.exp_vec_pieces,fill(nothing,d - l))
    end
    c.exp_vec_pieces[d] = gen_exp_vec(c.n,d,c.termorder,vars_reversed=c.vars_reversed)

    return
end

function generate_degree_reverse(c::PolyExpCache,d)
    d ≤ 0 && return 
    d in cached_reverse_lookups(c) && return 

    evs = get_forward(c,d)

    l = length(c.rev_look_pieces) 
    if l < d
        append!(c.rev_look_pieces,fill(nothing,d - l))
    end
    c.rev_look_pieces[d] = Dict(evs[i] => i for i in 1:length(evs))

    return
end

function PolyExpCache(n,termorder,degrees_to_prefill,reverse_to_prefill;vars_reversed=false)
    c = PolyExpCache(n,termorder)
    for d in degrees_to_prefill
        generate_degree_forward(c,d)
    end

    for d in reverse_to_prefill
        generate_degree_reverse(c,d)
    end
    c
end

function PolyExpCache(n,termorder,max_degree_prefill;vars_reversed=false)
    PolyExpCache(n,termorder,1:max_degree_prefill,1:max_degree_prefill)
end



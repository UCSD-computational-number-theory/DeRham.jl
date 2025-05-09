
"""
PEP stands for *p*artially *e*valuated (linear) *p*olynomial

f(V...,U...)

where f has already been evaluated at V.
f(V...,---) is stored as a vector of coefficients of the variables U,
with key V.

coefficients of the U polynomial are stored in order, with the constant term first,
i.e. 1 + 2U_1 + 3U_3 would be the vector [1,2,0,2]

coefficients of the U polynomials are of type T

The PEP interface
-----------------

PEP types are expected to implement the following methods

Base.getindex(V::Vector{Int})
allpoints()
allcomponents()
"""
abstract type AbstractPEP{T} end

"""
An EagerPEP{T} is a PEP with coefficients of type T.

An eager PEP computes all of its entries when it is created, and 
is does not compute after that.

In the constructor, one may use threads to faster compute all
of the values.
"""
struct EagerPEP{T} <: AbstractPEP{T} 
    Vs::Vector{Vector{Int}}
    # first term is constant, then each variable term
    Ucomponent::Dict{Vector{Int},Vector{T}} 

    function EagerPEP{T}(Vs, compute; usethreads=false) where T
        Ucomponent = Dict{Vector{Int},Vector{T}}()
        if usethreads
            Threads.@threads for V in Vs
                coeffs = compute(V)
                Ucomponent[V] = coeffs #TODO: make thread safe?
            end
        else
            for V in Vs
                coeffs = compute(V)
                Ucomponent[V] = coeffs 
            end
        end

        new{T}(Vs,UComponent)
    end
end

Base.getindex(P::EagerPEP, V::Vector{Int}) = P.Ucomponent[V] # catch error?

allpoints(P::EagerPEP) = P.Vs
allcomponents(P::EagerPEP) = P.Ucomponent


"""
A LazyPEP{T} is a PEP with coefficients of type T.

It computes each entry the first time it gets accessed.
"""
struct LazyPEP{T} <: AbstractPEP{T}
    Vs::Vector{Vector{Int}}
    Ucomponent::Dict{Vector{Int},Vector{T}}
    compute::Function

    function LazyPEP{T}(compute) where T
        Vs = Vector{Vector{Int}}()
        Ucomponent = Dict{Vector{Int},Vector{T}}()

        new{T}(Vs,Ucomponent,compute)
    end
end

allpoints(P::LazyPEP) = P.Vs
allcomponents(P::LazyPEP) = P.Ucomponent

function Base.getindex(P::LazyPEP, V::Vector{Int})
    if haskey(Ucomponent,V)
        return Ucomponent[V]
    end

    coeffs = P.compute(V)

    push!(P.Vs,V)
    P.Ucomponent[V] = coeffs

    return coeffs
end

#"""
#An CachePEP is a PEP with coefficients of type S.
#
#We assume that due to resource constraints we 
#cannot have very many terms of type S.
#However, we can readily have things available in type T,
#and we can convert between the two by taking S(t) for t
#some term of type T.
#
#The main example (for which this was designed) is when
#S is a type that lives on the GPU and T is a corresponding
#type that lives on the CPU.
#Due to memory constraints we might not be able to have everything
#allocated on the GPU at once, so we need to manage swapping over.
#"""
#struct CachePEP{T,S} <: AbstractPEP{S}
#    backing::AbstractPEP{T}
#    # Ucomponent = LRUCache{Vector{Int},Vector{S}}
#
#    function CachePEP{T,S}(backing)
#
#    end
#end
#
#allpoints(P::CachePEP) = allpoints(P.backing)
##cachedpoints(P::LRUCachePEP) = # keys
#allcomponents(P::CachePEP) = allcomponents(P.backing)
#
#function Base.getindex(P::CachePEP, V::Vector{Int})
#    if haskey(P.Ucomponent,V)
#        # cache hit!
#        return P.Ucomponent[V]
#    end
#
#    # cache miss!
#
#    coeffs = backing[V]
#
#    P.Ucomponent[V] = map(x -> S(x), coeffs)
#
#    return P.Ucomponent[V]
#end


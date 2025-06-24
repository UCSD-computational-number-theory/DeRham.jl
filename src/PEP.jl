

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
        #Ucomponent = Dict{Vector{Int},Vector{T}}()
        if usethreads
            Ucomponent = Dict{Vector{Int},Vector{T}}()
            Threads.@threads for V in Vs
                coeffs = compute(V)
                Ucomponent[V] = coeffs #TODO: make thread safe?
            end
        else
            Ucomponent = Dict{Vector{Int},Vector{T}}()
            for V in Vs
                coeffs = compute(V)
                Ucomponent[V] = coeffs 
            end
        end

        new{T}(Vs,Ucomponent)
    end
end

Base.getindex(P::EagerPEP, V::Vector{Int}) = P.Ucomponent[V] # catch error?

allpoints(P::EagerPEP) = P.Vs
allcomponents(P::EagerPEP) = P.Ucomponent


"""
A LazyPEP{T} is a PEP with coefficients of type T.

It computes each entry the first time it gets accessed.

Optionally, you may specify some Vs to be computed at creation time.
"""
struct LazyPEP{T} <: AbstractPEP{T}
    Vs::Vector{Vector{Int}}
    Ucomponent::Dict{Vector{Int},Vector{T}}
    compute::Function

    function LazyPEP{T}(compute;eagerVs=Vector{Vector{Int}}(),usethreads=false) where T
        Ucomponent = Dict{Vector{Int},Vector{T}}()

        if usethreads && 0 < length(eagerVs)
            #OhMyThreads.tforeach(eagerVs) do V 
            Threads.@threads for V in eagerVs
                coeffs = compute(V)
                Ucomponent[V] = coeffs #TODO: make thread safe?
            end
        else
            for V in eagerVs
                coeffs = compute(V)
                Ucomponent[V] = coeffs 
            end
        end

        new{T}(eagerVs,Ucomponent,compute)
    end
end

allpoints(P::LazyPEP) = P.Vs
allcomponents(P::LazyPEP) = P.Ucomponent

function Base.getindex(P::LazyPEP, V::Vector{Int})
    if haskey(P.Ucomponent,V)
        return P.Ucomponent[V]
    end

    coeffs = P.compute(V)

    push!(P.Vs,V)
    P.Ucomponent[V] = coeffs

    return coeffs
end

"""
An CachePEP is a PEP with coefficients of type S.

We assume that due to resource constraints we 
cannot have very many terms of type S.
However, we can readily have things available in type T,
and we can convert between the two by taking S(t) for t
some term of type T.

The main example (for which this was designed) is when
S is a type that lives on the GPU and T is a corresponding
type that lives on the CPU.
Due to memory constraints we might not be able to have everything
allocated on the GPU at once, so we need to manage swapping over.

This cache takes a lazy loading approach, putting elements into
the cache when they are first accessed.

The cache also keeps data in place

Ucomponent - the cache
backing - the AbstractPEP{T} which is the backing
create - create an S for the very first time
convert - function that converts a T to an S
"""
struct CachePEP{T,S} <: AbstractPEP{S}
    Ucomponent::LFUDA{Vector{Int},Vector{S}}
    backing::AbstractPEP{T}
    create::Function
    convert::Function
    temp::Base.RefValue{Union{Vector{S},Nothing}}
    tempV::Base.RefValue{Union{Vector{Int},Nothing}}

    function CachePEP{T,S}(backing,create,convert_entry,maxsize) where {T,S}
        temp = Ref{Union{Vector{S},Nothing}}(nothing)
        tempV = Ref{Union{Vector{Int},Nothing}}(nothing)
        recover = (key, value) -> recover!(tempV,temp,key,value)
        Ucomponent = LFUDA{Vector{Int},Vector{S}}(maxsize=maxsize,finalizer=recover)

        new{T,S}(Ucomponent,backing,create,convert_entry,temp,tempV)
    end
end

allpoints(P::CachePEP) = allpoints(P.backing)
cachedpoints(P::CachePEP) = keys(P.Ucomponent)

allcomponents(P::CachePEP) = allcomponents(P.backing)
cachedcomponents(P::CachePEP) = P.Ucomponent

function recover!(keyref::Ref,valueref::Ref,key::Vector,value::Vector)
    println("Recovering $key")
    keyref[] = key
    valueref[] = value
end

function Base.getindex(P::CachePEP, V::Vector{Int})

    default = () -> begin 
        if V == P.tempV[]
            t = P.temp[]
            P.temp[] = nothing
            P.tempV[] = nothing
            t
        elseif P.tempV[] == nothing
            P.create(P.backing[V])
        else
            # in this case, tempV is already initialized
            P.convert(P.temp[],P.backing[V])
            t = P.temp[]
            P.temp[] = nothing
            P.tempV[] = nothing
            t
        end
    end

    #println("Getting cache entry at $V")
    get!(default,P.Ucomponent,V)

    #if haskey(P.Ucomponent,V)
    #    # cache hit!
    #
    #    return P.Ucomponent[V]
    #end

    ## cache miss!

    #coeffs = backing[V]

    ##TODO: is there a version of this which won't need to allocate?
    ##Does that even matter?
    #P.Ucomponent[V] = map(x -> S(x), coeffs)

    #return P.Ucomponent[V]
end



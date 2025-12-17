
"""
Abstract type that describes a Newton Polygon/Hodge Polygon
like object.

There are a couple of ways to store such data, the idea
is that different implementations can be used to implement 
those different strategies
"""
abstract type AbstractPolygon end

"""
Represents a polygon that starts at
(0,0) and contines for slope slopes[i] 
for slopelengths[i] units.

For example, SlopesPolygon([1,19,1])
is the hodge polygon of 
(the primitive second cohomology of)
a K3 surface.

This is implemented in an "eager" way, i.e.
all relevant constructions are calculated
upon creation of the struct.

`slopes` - a list of slopes of the polygon
`slopelengths` - a list of how long each slope goes for
`values` - the value of the polygon as a function, evaluated starting at x=0
`slopesbefore` - given an x, what is the slope in the coordinate just before x,
                 by definition slopesbefore[1] = 0
"""
struct SlopesPolygon
    slopes::Array{Rational{Int}}
    slopelengths::Array{Int}
    values::Array{Rational{Int}}
    slopesbefore::Array{Rational{Int}}
end

# MARK - constructors

"""
Gives the array of values of the polygon,
starting with zero. Thus,
(i-1,values(sp)[i])
is a point on the polygon.

Right now, we call this in the constructors.
"""
function values(slopes,slopelengths)
    l = sum(slopelengths) + 1
    vals = zeros(Rational{Int},l)
    slopesbefore = zeros(Rational{Int},l)

    vals[1] = 0 # starting point is (0,0)
    slopesbefore[1] = 0

    nSlopes = length(slopelengths)
    j = 1 
    for i = 1:nSlopes
        slope = slopes[i]
        slopeend = j-1 + slopelengths[i]
        while j <= slopeend
            yval = vals[j] + slope 
            vals[j+1] = yval
            slopesbefore[j+1] = slope
            j = j + 1
        end
    end

    (vals,slopesbefore)
end

"""
Given an array of vertices
"""
#function SlopesPolygon(vertices::Array{Tuple{Int,Int}})
function SlopesPolygon(vertices::Vector{Tuple{Int64, Int64}})
    #TODO implement
    
    n = length(vertices)-1
    slopelengths = zeros(Int,n)
    slopes = zeros(Rational{Int},n)
    for i = 2:n+1
        slopevec = vertices[i] .- vertices[i-1]
        #push!(slopes, slopevec[2] // slopevec[1])
        #push!(slopelengths,slopevec[2])
        slopes[i-1] = slopevec[2] // slopevec[1]
        #slopelengths[i-1] = slopevec[2]
        slopelengths[i-1] = slopevec[1]
    end
    
    SlopesPolygon(slopes,slopelengths,values(slopes,slopelengths)...)
end

"""
This constructor assumes integer slopes
starting at zero, and slope i goes for
slopelengths[i] units

E.g. SlopesPolygon([1,19,1]) is the hodge polygon of (the primitive H^2 of) a K3 surface.
"""
function SlopesPolygon(slopelengths::Array{T}) where T<:Integer
    n = length(slopelengths) - 1
    slopes = Rational.(collect(0:n))
    SlopesPolygon(slopes,slopelengths,values(slopes,slopelengths)...)
end

"""
Given an array of points
(xs[i],yvals[i]),
creates a newton polygon by taking the lower
convex hull of the points.

This should work for newton and hodge polygons of varieties,
but I don't know exactly what I'm assuming about the
points here.

This assumes that xs is sorted in ascending order
"""
function SlopesPolygon(xs::Array{T},yvals::Array{T}) where T<:Real
    vertices = [(0,0)]

    points = [(xs[i],yvals[i]) for i in eachindex(yvals)]
    
    x = 0
    while x < xs[end] # i.e. end of the hodge polygon

        startpoint = points[x+1]
        x2 = x + 1
        smallest_slope = typemax(Int128)
        finalx = x

        while x2 ≤ length(yvals) - 1

            thispoint = points[x2 + 1]
            thisslope = (thispoint[2] - startpoint[2]) // (thispoint[1] - startpoint[1])
            if thisslope ≤ smallest_slope # we do \leq even though the slope won't change so we can update finalpoint
              smallest_slope = thisslope
              finalx = x2
            end
            x2 = x2 + 1
        end # postcondition: x2 is the next vertex

        # add x2 to the vertices
        push!(vertices,points[finalx+1])

        x = finalx
    end

    SlopesPolygon(vertices)
end

#"""
#Given a list of coefficients of a polynomial,
#starting with the constant term and ending with the
#highest degree term, gives the newton polygon
#using a very naive algorithm.

#"""
#function newton_polygon(p,coefs)
#  vertices = [(0,0)]

#  points = [(i-1,padic_val(2,coefs[i])) for i in eachindex(coefs)]

#  x = 0
#  while x < length(coefs) - 1 # i.e. the end of the hodge polygon
#    #println("x=$x")
#    startpoint = points[x+1]
#    x2 = x + 1
#    smallest_slope = typemax(Int128)
#    finalx = x
#    while x2 ≤ length(coefs) - 1
#      #println("x2=$x2, y1=$(points[x2 + 1][2])")
#      thispoint = points[x2 + 1]
#      thisslope = (thispoint[2] - startpoint[2]) // (thispoint[1] - startpoint[1])
#      if thisslope ≤ smallest_slope # we do \leq even though the slope won't change so we can update finalpoint
#        smallest_slope = thisslope
#        finalx = x2
#      end
#      x2 = x2 + 1
#    end # postcondition: x2 is the next vertex

#    # add x2 to the vertices
#    push!(vertices,points[finalx+1])

#    x = finalx

#  end

#  vertices
#end

function SlopesPolygon(coefficients,valuation)
    n = length(coefficients)-1
    SlopesPolygon(collect(0:n),reverse!(valuation.(coefficients)))
end


function newton_polygon(p,coeffs)
    #valuation = x -> x == 0 ? typemax(typeof(x)) : padic_val(p,x)
    #padic_val = x -> x == 0 ? typemax(typeof(x)) : valuation(x, p)
    padic_val = x -> x == 0 ? typemax(Int64) : valuation(x, p) # TODO: make this correct

    SlopesPolygon(coeffs,padic_val)
end

# MARK - creating polygons from other polygons

tatetwist(sp::SlopesPolygon,n) = SlopesPolygon(sp.slopes .- n,sp.slopelengths,nothing)

# MARK - methods

slopelengths(sp::SlopesPolygon) = sp.slopelengths

function Base.getindex(sp::SlopesPolygon, i::Integer)
    sp.values[i+1]
end

function Base.:(==)(x::SlopesPolygon,y::SlopesPolygon) 
    (x.slopes == y.slopes) && (x.slopelengths == y.slopelengths)
end

function Base.hash(x::SlopesPolygon, h::UInt64)
    hash((x.slopes, x.slopelengths))
end

function removetypefromstring(arr) 
    s = string(arr)
    i = findfirst(==('['), s)
    s[i:end]
end

function Base.show(io::IO, sp::SlopesPolygon)
    println(io, "SlopesPolygon of length $(length(sp.values) - 1)")


    println(io, "Slopes:        $(removetypefromstring(sp.slopes))")
    println(io, "Slope Lengths: $(sp.slopelengths)")
    println(io, "Values:        $(removetypefromstring(sp.values))")

    plt = lineplot(collect(0:(length(sp.values) - 1)),
                   sp.values,
                   # width=:auto,
                   width=60,
                   height=22,
                   compact=true)
    println(io, plt)
end

endcoord(sp::SlopesPolygon) = length(values(sp)) - 1
    
function vertices(sp::SlopesPolygon)
    error("not implemented")
end

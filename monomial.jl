using Oscar


#Very rough code to generate monomials, will clean up later


#Generates an array of vectors corresponding to the degrees of the variables in the monomial
function genVec(n,d)
    z = Any[]
    if n == 1
        return [[d]]
    end
    if d == 1
        for i in 1:n
            s = zeros(Int64,n)
            s[i] = 1
            push!(z,s)
        end
        return z
    end
    for i in 0:d
        y = genVec(n-1,d-i)
        for j in axes(y,1)
            append!(y[j],i)
        end
        append!(z,y)
    end
    return z
end

#generates the monomials
function genMon(n,d)
    m = []
    z = genVec(n,d)
    for i in axes(z,1)
        B = MPolyBuildCtx(R)
        push_term!(B, QQ(1), z[i])
        p = finish(B)
        push!(m,p)
    end
    return m
end

#=
n = number of variables
d = degree of the monomials
=#

n = 3
d = 3
R, x = polynomial_ring(QQ, "x" => (1:n))

y = genMon(n,d)

for i in axes(y,1)
    println(y[i])
end
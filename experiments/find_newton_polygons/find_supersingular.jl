
"""
pencil - function that takes in an element of GF(p) and returns a 
EXAMPLE:
dwork(a) = x^4 + y^4 + z^4 + w^4 + a*x*y*z*w
p - the prime
"""
function evaluate_pencil(pencil,p)
    for i in 0:p-1
        f = pencil(i)
        np = DeRham.newton_polygon(f,fastevaluation=true,algorithm=:naive)
        println("$i: $np")
    end
end


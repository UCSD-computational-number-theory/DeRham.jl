
using CSV
using DataFrames

function test_filepath(filename)
    prefix = "varieties/"

    if basename(pwd()) != "test"
        prefix = "test/" * prefix
    end

    prefix * filename
end

"""
Saves a test in a csv file with name "filename" 
in the test/varieties folder

"""
function savetest(n,p,f,zeta,filename)

    #TODO: use pwd() to decide if we're in
    #   the test folder are not, if not, 
    #   append the prefix test/ to the filename

    filepath = test_filepath(filename)


    if !isfile(filepath)
        open(filepath,"a") do file
            write(file,"n;p;f;zeta_function\n")
        end
    end

    open(filepath,"a") do file
        write(file,"$n;$p;$f;$(ZZ.(vec(zeta)))\n")
    end
end

"""
    runcsvtest(filename; zeta_function=nothing)

Runs all the tests in the file test/variables/[filename]

optionally, you may specify a particular zeta function function
by using the optional zeta_function argument.

Examples:

runtest("ellipticcurves.csv")

zf = f -> DeRham.zeta_function(f,fastevaluation=true)
runtest("ellipticcurves.csv"; zeta_function=zf)
"""
function runcsvtest(filename; zeta_function=nothing)

    filepath = test_filepath(filename)
    
    tests = CSV.read(filepath,DataFrame,delim=";")

    if zeta_function == nothing
        zeta_function = DeRham.zeta_function
    end

    for test in eachrow(tests)
        n = test.n
        varstring = prod(["x$i," for i in 1:n])[1:end-1]

        # this isn't secure!
        codestring = """begin
            using Oscar
            R, ($varstring) = polynomial_ring(GF($(test.p)),$n);
            $(test.f)
        end
        """

        f = eval(Meta.parse(codestring))
        zeta_result = eval(Meta.parse("begin using Oscar; $(test.zeta_function) end"))

        zeta_evaluated = zeta_function(f)
        println()
        println("Testing: p=$(test.p), $n variables, f = $f")
        println()
        @test zeta_result == zeta_evaluated
    end

end


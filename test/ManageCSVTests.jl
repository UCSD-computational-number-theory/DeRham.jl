
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

runcsvtest("ellipticcurves.csv")

zf = f -> DeRham.zeta_function(f,fastevaluation=true)
runcsvtest("ellipticcurves.csv"; zeta_function=zf)
"""
function runcsvtest(filename; zeta_function=nothing)

    filepath = test_filepath(filename)
    
    tests = CSV.read(filepath,DataFrame,delim=";")

    if zeta_function == nothing
        zeta_function = DeRham.zeta_function
    end

    for test in eachrow(tests)
        testzetafunction(test.n,test.p,test.f,test.zeta_function,zeta_function)
    end

end

function teststring(ts, zeta_function=nothing)
    if zeta_function == nothing
        zeta_function = DeRham.zeta_function
    end

    testparts = split(ts,(";"))

    n = Meta.parse(testparts[1])
    p = Meta.parse(testparts[2])
    testzetafunction(n,p,testparts[3],testparts[4],zeta_function)
end

function testzetafunction(n,p,fstring,correctzeta,zeta_function)
    varstring = prod(["x$i," for i in 1:n])[1:end-1]

    # this isn't secure!
    codestring = """begin
        using Oscar
        R, ($varstring) = polynomial_ring(GF($(p)),$n);
        $(fstring)
    end
    """

    f = eval(Meta.parse(codestring))
    zeta_result = eval(Meta.parse("begin using Oscar; $(correctzeta) end"))

    println("Testing: p=$(p), $n variables, f = $f")
    println()
    zeta_evaluated = zeta_function(f)
    println()
    @test zeta_result == zeta_evaluated
end

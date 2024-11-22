
function hodge_polygon(p,n,d)
    if n == 2 # curve
        g = divexact((d-1)*(d-2),2)
        DeRham.SlopesPolygon([g,g])
    elseif n == 3 && d == 4
        DeRham.SlopesPolygon([1,19,1]) #normally [1,20,1]
    elseif n == 3 && d == 5
        DeRham.SlopesPolygon([4,44,4]) #normally [4,45,4]
    elseif n == 3 && d == 6 
        DeRham.SlopesPolygon([10,85,10]) #normally [10,86,10]
    else
        error("uniplemented newton polygon formula")
    end
end

# The following code is dual licensed with the GPL,
# since it originally came from being copied/pasted
# from controlledrecution [insert link here]
#
# MARK - Begin GPL dual license

# TODO: copy of gpl license

function test_series_precision_example(p,n,d)
    r_m = DeRham.calculate_relative_precision(hodge_polygon(p,n,d),n-1,p)
    N = DeRham.calculate_series_precision(p,n,r_m)

    if ( n == 2 && d == 3 ) 
        if ( p < 5) 
            @test N == [2 3]
        elseif ( 5 <= p && p < 17 ) 
            @test N == [2,2]
        elseif ( 17 <= p  ) 
            @test N == [0,1]
        end
    elseif ( n == 2 && d == 4 ) 
        if ( p < 5) 
            @test N == [4,4]
        elseif ( 5 <= p && p < 17 ) 
            @test N == [3,3]
        elseif ( 17 <= p  ) 
            @test N == [2,2]
        end
    elseif ( n == 2 && d == 5 ) 
        if ( p < 5) 
            @test N == [6,6]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [4,5]
        elseif ( 7 <= p  ) 
            @test N == [4,4]
        end
    elseif ( n == 2 && d == 6 ) 
        if ( p < 5) 
            @test N == [8,9]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [7,7]
        elseif ( 7 <= p && p < 11 ) 
            @test N == [6,7]
        elseif ( 11 <= p  ) 
            @test N == [6,6]
        end
    elseif ( n == 2 && d == 7 ) 
        if ( p < 5) 
            @test N == [11,11]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [10,10]
        elseif ( 7 <= p && p < 11 ) 
            @test N == [10,10]
        elseif ( 11 <= p && p < 17 ) 
            @test N == [9,9]
        elseif ( 17 <= p  ) 
            @test N == [8,8]
        end
    elseif ( n == 2 && d == 8 ) 
        if ( p < 5) 
            @test N == [14,14]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [13,13]
        elseif ( 7 <= p && p < 13 ) 
            @test N == [13,13]
        elseif ( 13 <= p && p < 17 ) 
            @test N == [12,13]
        elseif ( 17 <= p  ) 
            @test N == [11,11]
        end
    elseif ( n == 2 && d == 9 ) 
        if ( p < 5) 
            @test N == [18,18]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [16,16]
        elseif ( 7 <= p && p < 11 ) 
            @test N == [16,16]
        elseif ( 11 <= p && p < 17 ) 
            @test N == [16,16]
        elseif ( 17 <= p  ) 
            @test N == [15,15]
        end
    elseif ( n == 3 && d == 4 ) 
        if ( p < 5) 
            @test N == [7,7,8]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [4,5,5]
        elseif ( 7 <= p && p < 23 ) 
            @test N == [4,4,3]
        elseif ( 23 <= p && p < 43 ) 
            @test N == [3,3,3]
        elseif ( 43 <= p  ) 
            @test N == [3,3,2]
        end
    elseif ( n == 3 && d == 5 ) 
        if ( p < 5) 
            @test N == [11,11,10]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [8,8,9]
        elseif ( 7 <= p && p < 11 ) 
            @test N == [8,8,7]
        elseif ( 11 <= p && p < 23 ) 
            @test N == [7,7,6]
        elseif ( 23 <= p && p < 29 ) 
            @test N == [6,6,6]
        elseif ( 29 <= p  ) 
            @test N == [6,6,5]
        end
    elseif ( n == 3 && d == 6 ) 
        if ( p < 5) 
            @test N == [17,18,17]
        elseif ( 5 <= p && p < 7 ) 
            @test N == [15,15,14]
        elseif ( 7 <= p && p < 11 ) 
            @test N == [15,15,14]
        elseif ( 11 <= p && p < 17 ) 
            @test N == [14,14,13]
        elseif ( 17 <= p && p < 23 ) 
            @test N == [13,13,12]
        elseif ( 23 <= p  ) 
            @test N == [12,12,11]
        end
    end
    
end


function test_algorithm_precision_example(p,n,d)

    r_m = DeRham.calculate_relative_precision(hodge_polygon(p,n,d),n-1,p)
    N_m = DeRham.series_precision(p,n,d, r_m)

    precision = DeRham.algorithm_precision(p,n,d,r_m,N_m)

    if ( n == 2 && d == 3 ) 
        if ( p < 5) 
            @test precision == 5;
        elseif ( 5 <= p && p < 17 ) 
            @test precision == 3;
        elseif ( 17 <= p  ) 
            @test precision == 1;
        end
    elseif ( n == 2 && d == 4 ) 
        if ( p < 5) 
            @test precision == 7;
        elseif ( 5 <= p && p < 17 ) 
            @test precision == 5;
        elseif ( 17 <= p  ) 
            @test precision == 3;
        end
    elseif ( n == 2 && d == 5 ) 
        if ( p < 5) 
            @test precision == 12;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 9;
        elseif ( 7 <= p  ) 
            @test precision == 7;
        end
    elseif ( n == 2 && d == 6 ) 
        if ( p < 5) 
            @test precision == 19;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 13;
        elseif ( 7 <= p && p < 11 ) 
            @test precision == 13;
        elseif ( 11 <= p  ) 
            @test precision == 11;
        end
    elseif ( n == 2 && d == 7 ) 
        if ( p < 5) 
            @test precision == 23;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 20;
        elseif ( 7 <= p && p < 11 ) 
            @test precision == 19;
        elseif ( 11 <= p && p < 17 ) 
            @test precision == 17;
        elseif ( 17 <= p  ) 
            @test precision == 15;
        end
    elseif ( n == 2 && d == 8 ) 
        if ( p < 5) 
            @test precision == 30;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 26;
        elseif ( 7 <= p && p < 13 ) 
            @test precision == 25;
        elseif ( 13 <= p && p < 17 ) 
            @test precision == 25;
        elseif ( 17 <= p  ) 
            @test precision == 21;
        end
    elseif ( n == 2 && d == 9 ) 
        if ( p < 5) 
            @test precision == 41;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 33;
        elseif ( 7 <= p && p < 11 ) 
            @test precision == 32;
        elseif ( 11 <= p && p < 17 ) 
            @test precision == 31;
        elseif ( 17 <= p  ) 
            @test precision == 29;
        end
    elseif ( n == 3 && d == 4 ) 
        if ( p < 5) 
            @test precision == 16;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 9;
        elseif ( 7 <= p && p < 23 ) 
            @test precision == 6;
        elseif ( 23 <= p && p < 43 ) 
            @test precision == 5;
        elseif ( 43 <= p  ) 
            @test precision == 4;
        end
    elseif ( n == 3 && d == 5 ) 
        if ( p < 5) 
            @test precision == 21;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 17;
        elseif ( 7 <= p && p < 11 ) 
            @test precision == 14;
        elseif ( 11 <= p && p < 23 ) 
            @test precision == 12;
        elseif ( 23 <= p && p < 29 ) 
            @test precision == 11;
        elseif ( 29 <= p  ) 
            @test precision == 10;
        end
    elseif ( n == 3 && d == 6 ) 
        if ( p < 5) 
            @test precision == 38;
        elseif ( 5 <= p && p < 7 ) 
            @test precision == 29;
        elseif ( 7 <= p && p < 11 ) 
            @test precision == 28;
        elseif ( 11 <= p && p < 17 ) 
            @test precision == 26;
        elseif ( 17 <= p && p < 23 ) 
            @test precision == 24;
        elseif ( 23 <= p  ) 
            @test precision == 22;
        end
    else
        println("Example not supported by `controlledreduction`")
        @test false
    end
end

# MARK - End GPL dual licesnse

function test_algorithm_precision()
    #TODO: figure out precision for small primes
    #primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    curveprimes = [7,11,13]

    for p in curveprimes

        test_algorithm_precision_example(p,2,3)
        test_algorithm_precision_example(p,2,4)
        test_algorithm_precision_example(p,2,5)
        test_algorithm_precision_example(p,2,6)
        test_algorithm_precision_example(p,2,7)
        test_algorithm_precision_example(p,2,8)
        test_algorithm_precision_example(p,2,9)

    end

    primes = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    for p in primes 

        #TODO: consolidate preivous for loop here 
        #      once curves are working w.r.t. precision
        #
        #test_algorithm_precision_example(p,2,3)
        #test_algorithm_precision_example(p,2,4)
        #test_algorithm_precision_example(p,2,5)
        #test_algorithm_precision_example(p,2,6)
        #test_algorithm_precision_example(p,2,7)
        #test_algorithm_precision_example(p,2,8)
        #test_algorithm_precision_example(p,2,9)

        test_algorithm_precision_example(p,3,4)
        test_algorithm_precision_example(p,3,5)
    end
end

function test_series_precision()
    #TODO: add the bound and fix for small p

    primes = [3,5,7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    for p in primes
        n = 2
        for d in 3:9
            println("Testing precision estimate for (p,n,d) = ($p,$n,$d)")
            test_series_precision_example(p,n,d)
        end

        n = 3
        for d in 4:6
            println("Testing precision estimate for (p,n,d) = ($p,$n,$d)")
            test_series_precision_example(p,n,d)
        end
    end

    #test_series_precision_example(5,2,3)
    #test_series_precision_example(5,2,4)

    #test_series_precision_example(7,2,3)
    #test_series_precision_example(7,2,4)
    #
    #test_series_precision_example(7,3,4)
    #

    ## p = 11
    #test_series_precision_example(11,2,3)
    #test_series_precision_example(11,2,4)

    ## p = 13
    #test_series_precision_example(13,2,3)
    #test_series_precision_example(13,2,4)
    #test_series_precision_example(13,2,5)

    ## k3 surfaces
    #test_series_precision_example(11,3,4)
    #test_series_precision_example(13,3,4)

    #test_series_precision_example(13,3,5)

    #test_series_precision_example(17,2,5)
    #test_series_precision_example(17,2,6)
end


function test_hodge_polygon_values()
    k3 = DeRham.SlopesPolygon([1,20,1])

    @test k3[1] == 0
    @test k3[2] == 1
    @test k3[10] == 9
    @test k3[21] == 20
    @test k3[22] == 22

    genus5= DeRham.SlopesPolygon([5,5])

    @test genus5[0] == 0
    @test genus5[2] == 0
    @test genus5[5] == 0
    @test genus5[6] == 1
    @test genus5[8] == 3
    @test genus5[10] == 5
end

function test_hodge_polygon_examples()
    R, (x,y,z,w) = polynomial_ring(GF(7),4)

    f = x^4 + y^4 + z^4 + w^4

    @test DeRham.hodgepolygon(f) == DeRham.SlopesPolygon([1,19,1])
end


function test_ibm_crank()

end

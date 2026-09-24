
function test_hodge_polygon_values()
    k3 = DeRham.SlopesPolygon([1, 20, 1])

    @test k3[1] == 0
    @test k3[2] == 1
    @test k3[10] == 9
    @test k3[21] == 20
    @test k3[22] == 22

    genus5 = DeRham.SlopesPolygon([5, 5])

    @test genus5[0] == 0
    @test genus5[2] == 0
    @test genus5[5] == 0
    @test genus5[6] == 1
    @test genus5[8] == 3
    @test genus5[10] == 5
end

function test_hodge_polygon_examples()
    R, (x, y, z, w) = polynomial_ring(GF(7), 4)

    f = x^4 + y^4 + z^4 + w^4

    @test DeRham.hodgepolygon(f) == DeRham.SlopesPolygon([1, 19, 1])
end

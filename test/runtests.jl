using Exact
using Test

@testset "Exact model values" begin
    @test getgse(TFIC(1.0, 1.0))[1] ≈ -0.5317722049866825 atol=1e-12
    @test getfe(TFIC(1.0, 1.0), 1.0)[1] ≈ -0.8398140716988789 atol=1e-12
    @test getie(TFIC(1.0, 1.0), 1.0)[1] ≈ -0.2761925884297794 atol=1e-12
    @test getgse(TFIC(1.0, 0.0))[1] == -0.25
    @test isfinite(getfe(TFIC(1.0, 0.0), 1.0)[1])

    @test getgse(TFXYC(1.0, 0.8, 0.2))[1] ≈ -0.29896218401601693 atol=1e-12

    xxz_energy = getgse(XXZC(1.0, 0.5))[1]
    @test xxz_energy isa Real
    @test xxz_energy ≈ -1.8750000000001146 atol=1e-10

    @test getfe(IsingSquare(), 2.5)[1] ≈ -2.198409551937308 atol=1e-9
    @test getfe(IsingKagome(), 2.5)[1] ≈ -6.086454949666 atol=1e-7
    @test getfe(IsingTriangular(), 2.5)[1] ≈ -2.126184176901414 atol=1e-8
end

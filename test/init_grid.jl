@testset "Grid" begin

    nx, ny = 8, 8
    grid = Grid( 1e6, nx, 1e6, ny )

    tx = SurfaceQuasiGeostrophic.fct_unity_approx5(nx)
    ty = SurfaceQuasiGeostrophic.fct_unity_approx5(ny)
    
    @show tx
    @show ty

    @test true

end

@testset "Poisson" begin

    xmin = 0.0
    xmax = 2π
    ymin = 0.0
    ymax = 2π

    nx = 128
    ny = 128

    ex = zeros(Float64, nx, ny)
    ey = zeros(Float64, nx, ny)

    rho = zeros(Float64, nx, ny)
    phi = zeros(Float64, nx, ny)

    ex_exact = zeros(nx, ny)
    ey_exact = zeros(nx, ny)
    phi_exact = zeros(nx, ny)

    grid = Grid(xmax-xmin, nx, ymax - ymin, ny)

    poisson = Poisson(grid)

    mode = 2
    for i = 1:nx, j = 1:ny
        x = xmin + (i - 1) * (xmax - xmin) / nx
        y = ymin + (j - 1) * (ymax - ymin) / ny
        phi_exact[i, j] = -mode * sin(mode * x) * cos(mode * y)
        ex_exact[i, j] = mode^2 * cos(mode * x) * cos(mode * y)
        ey_exact[i, j] = -mode^2 * sin(mode * x) * sin(mode * y)
        rho[i, j] = -2 * mode^3 * sin(mode * x) * cos(mode * y)
    end

    rhs = copy(rho)
    solve!(phi, poisson, rhs)
    @test phi_exact ≈ phi

    rhs = copy(rho)
    solve!(ex, ey, poisson, rhs)
    @test ex_exact ≈ ex
    @test ey_exact ≈ ey

end

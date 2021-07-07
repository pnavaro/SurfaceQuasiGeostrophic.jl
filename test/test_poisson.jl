@testset "Poisson" begin

    xmin = 0.0
    xmax = 2π
    ymin = 0.0
    ymax = 2π

    nx = 128
    ny = 128

    ux = zeros(Float64, nx, ny)
    uy = zeros(Float64, nx, ny)

    omega = zeros(Float64, nx, ny)
    psi = zeros(Float64, nx, ny)

    ux_exact = zeros(nx, ny)
    uy_exact = zeros(nx, ny)
    psi_exact = zeros(nx, ny)

    grid = Grid(xmax-xmin, nx, ymax - ymin, ny)

    poisson = Poisson(grid)

    mode = 2
    for i = 1:nx, j = 1:ny
        x = xmin + (i - 1) * (xmax - xmin) / nx
        y = ymin + (j - 1) * (ymax - ymin) / ny
        psi_exact[i, j] = mode * sin(mode * x) * cos(mode * y)
        ux_exact[i, j] = mode^2 * sin(mode * x) * sin(mode * y)
        uy_exact[i, j] = mode^2 * cos(mode * x) * cos(mode * y)
        omega[i, j] = 2 * mode^3 * sin(mode * x) * cos(mode * y)
    end

    rhs = copy(omega)
    compute_psi!(psi, poisson, rhs)
    @test psi_exact ≈ psi

    rhs = copy(omega)
    compute_velocities!(ux, uy, poisson, rhs)
    @test ux_exact ≈ ux
    @test uy_exact ≈ uy

end

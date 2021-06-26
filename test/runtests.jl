using SurfaceQuasiGeostrophic
using Test

@testset "Grid" begin

    nx, ny = 8, 8
    grid = Grid( 1e6, nx, 1e6, ny )

    tx = SurfaceQuasiGeostrophic.fct_unity_approx5(nx)
    ty = SurfaceQuasiGeostrophic.fct_unity_approx5(ny)
    
    @show tx
    @show ty

    @test true

end

@testset "Init Buoyancy" begin

# Duration of the simulation [in seconds]
advection_duration = 3600*24*20; # 20 days

# Resolution [even integer]
nx, ny  = 64, 64
angle_grid = π/4 
omega = 2π/(24*60*60)
f0 = 2 * omega * sin( angle_grid )

grid = Grid( 1e6, nx, 1e6, ny )

# Gather parameters in the structure model
model = ( 
    rho=1e3, # Background density
    g = 9.81, # Gravity
    buoyancy_freq_N = 3 * f0, # Background stratification
    odg_b = 1e-3 # Amplitude of the buoyancy
)

# Initial condition for the buoyancy
fft_buoy = fct_buoyancy_init(model, grid)

# Advection

# fct_fft_advection_sto!(model, fft_buoy)

@test true

end

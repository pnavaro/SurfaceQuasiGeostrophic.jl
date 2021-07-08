using SurfaceQuasiGeostrophic
using Test
using FFTW
using Plots
using ProgressMeter

ENV["GKSwstype"] = "100"

include("init_grid.jl")
include("test_poisson.jl")
include("init_buoyancy.jl")

@testset "Advection" begin

# Duration of the simulation [in seconds]
advection_duration = 3600*24*20; # 20 days

# Resolution [even integer]
nx, ny  = 64, 64
angle_grid = π/4 
omega = 2π/(24*60*60)
f0 = 2 * omega * sin( angle_grid )

## Grid [spatial & Fourier]
grid = Grid( 1e6, nx, 1e6, ny )

# Gather parameters in the SQG model
sqg = SQG( grid, f0 ) 

# Initial condition for the buoyancy
init_buoyancy!(sqg)

@show sqg.buoyancy_freq_n

## Initialisation of the spatial fields

update_velocities!( sqg)
contourf(sqg.u_x)
savefig("u_x.png")
contourf(sqg.u_y)
savefig("u_y.png")
contourf(irfft(sqg.ψ̂, nx))
savefig("stream_function.png")

update_advection_term!( sqg )

contourf(irfft(sqg.â, nx))
savefig("advection_term.png")

## Hyperviscosity [order & coefficient]
@show sqg.hv_order = 8
@show sqg.hv_val = compute_hypercoef(sqg)

## Choice of time step : CFL
@show dt = compute_dt(sqg)
println("Time step: $dt seconds")

nstep = advection_duration ÷ dt
@show nstep

# Printing some information
println("Final time: $(nstep*dt÷(3600*24)) days ")

rk4 = TimeSolver(nx, ny)

@showprogress 1000 for i = 1:nstep
    step!( sqg, rk4, dt )
end

contourf(irfft(sqg.b̂, nx))
savefig("buoyancy.png")

@test true

end

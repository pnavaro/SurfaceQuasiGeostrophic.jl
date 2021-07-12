using SurfaceQuasiGeostrophic
using Test
using FFTW
using ProgressMeter
using TimerOutputs

const to = TimerOutput()

@testset "Advection" begin

# Duration of the simulation [in seconds]
advection_duration = 3600*24; # 1 day

# Resolution [even integer]
nx, ny  = 128, 128
lx, ly = 1e6, 1e6
angle_grid = π/4 
omega = 2π/(24*60*60)
f0 = 2 * omega * sin( angle_grid )

## Grid [spatial & Fourier]
grid = Grid( lx, nx, ly, ny )

# Gather parameters in the SQG model
sqg = SQG( grid, f0 ) 

# Initial condition for the buoyancy
init_buoyancy!(sqg)

## Initialisation of the spatial fields

update_velocities!( sqg)

## Hyperviscosity [order & coefficient]
@show sqg.hv_order = 8
@show sqg.hv_val = compute_hypercoef(sqg)

## Choice of time step : CFL
@show dt = compute_dt(sqg)
println("Time step: $dt seconds")

nstep = floor(Int, advection_duration / dt)
@show nstep

# Printing some information
println("Final time: $(nstep*dt÷(3600*24)) days ")

rk4 = TimeSolver(nx, ny)

pbar = Progress(nstep)

for i = 1:nstep

    rk4.b̂ .= sqg.b̂

    @timeit to "advection" update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt / 2
    rk4.db̂ .= sqg.â ./ 2

    @timeit to "advection" update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt / 2
    rk4.db̂ .+= sqg.â

    @timeit to "advection" update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt 
    rk4.db̂ .+= sqg.â

    @timeit to "advection" update_advection_term!(sqg)

    rk4.db̂ .+= sqg.â ./ 2

    sqg.b̂ .= rk4.b̂ .+ dt / 3 .* rk4.db̂


    next!(pbar)

end 

show(to)

@test true

end

include("init_grid.jl")
include("test_poisson.jl")
include("init_buoyancy.jl")

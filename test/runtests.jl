using SurfaceQuasiGeostrophic
using Test
using FFTW
using Plots
using ProgressMeter
using TimerOutputs

const to = TimerOutput()

ENV["GKSwstype"] = "100"


@testset "Advection" begin

# Duration of the simulation [in seconds]
advection_duration = 3600*24*20; # 20 days

# Resolution [even integer]
nx, ny  = 64, 64
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

contourf(grid.x, grid.y, irfft(sqg.b̂, nx), aspect_ratio = :equal)
savefig("initial_buoyancy.png")

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

pbar = Progress(nstep)

# anim = @animate for i = 1:nstep
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

    # contourf!(grid.x, grid.y, irfft(sqg.b̂, nx), aspect_ratio=:equal, 
    #           axis=([], false), colorbar=false, clims=(-sqg.odg_b,sqg.odg_b))

    next!(pbar)

end 

# end every 10

# gif(anim, joinpath(@__DIR__, "buoyancy.gif"), fps = 10)

show(to)

contourf(grid.x, grid.y, irfft(sqg.b̂, nx), aspect_ratio = :equal)
savefig("final_buoyancy.png")

@test true


end

include("init_grid.jl")
include("test_poisson.jl")
include("init_buoyancy.jl")

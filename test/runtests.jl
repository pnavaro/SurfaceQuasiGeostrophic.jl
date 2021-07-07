using SurfaceQuasiGeostrophic
using Test
using FFTW
using Plots

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


function update_advection_term!(sqg)

    update_velocities!( sqg )

    # Grid of wave vectors

    px = sqg.grid.nx ÷ 2 # frequency of aliasing
    py = sqg.grid.ny ÷ 2 # frequency of aliasing

    # Advection term

    sqg.adv .=  sqg.u_x .* ifft(-1im .* sqg.grid.kx  .* sqg.b̂) 
    sqg.adv .+= sqg.u_y .* ifft(-1im .* sqg.grid.ky' .* sqg.b̂)

    fft!(sqg.adv)
    sqg.adv[px+1, :] .= 0 # Remove aliasing
    sqg.adv[:, py+1] .= 0 # Remove aliasing

    # Summing Hyperviscosity term

    sqg.adv .+= - sqg.hv_val .* (sqg.grid.k.^sqg.hv_order) .* sqg.b̂

end

nstep = 10

for istep = 1:nstep

    b̂_old = sqg.b̂
    db̂ = zero(sqg.b̂)

    update_advection_term!(sqg)

    sqg.b̂ .= b̂_old .+ sqg.adv .* dt / 2
    db̂ .+= sqg.adv ./ 2

    update_advection_term!(sqg)

    sqg.b̂ .= b̂_old .+ sqg.adv * dt / 2
    db̂ .+= sqg.adv

    update_advection_term!(sqg)

    sqg.b̂ .= b̂_old .+ sqg.adv * dt 
    db̂ .+= sqg.adv

    update_advection_term!(sqg)

    db̂ .+= sqg.adv ./ 2

    sqg.b̂ .= b̂_old .+ dt / 3 .* db̂

end


@test true

end

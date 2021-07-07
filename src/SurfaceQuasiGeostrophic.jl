module SurfaceQuasiGeostrophic

using FFTW

include("grid.jl")

export SQG

"""
    SQG

Data structure to store arrays and parameters for the
Surface-Quasi-Geostrophic model

- rho : Background density
- g : Gravity
- buoyancy_freq_n : Background stratification
- odg_b : Amplitude of the buoyancy
"""
mutable struct SQG

    grid :: Grid
    rho :: Float64
    g :: Float64
    buoyancy_freq_n :: Float64
    odg_b :: Float64
    hv_order :: Int
    hv_val :: Float64

    b :: Array{Float64, 2}
    ψ :: Array{Float64, 2}
    b̂ :: Array{ComplexF64, 2}
    ψ̂ :: Array{ComplexF64, 2}

    u_x :: Array{Float64, 2}
    u_y :: Array{Float64, 2}
    û_x :: Array{ComplexF64, 2}
    û_y :: Array{ComplexF64, 2}

    adv :: Array{ComplexF64, 2}

    function SQG( grid, f0 )

        rho=1e3
        g = 9.81
        buoyancy_freq_n = 3 * f0
        odg_b = 1e-3 
        hv_order = 8
        hv_val = 0

        b = zeros(grid.nx, grid.ny)
        ψ = zeros(grid.nx, grid.ny)
        b̂ = complex(b)
        ψ̂ = complex(ψ)

        u_x = similar(ψ)
        u_y = similar(ψ)
        û_x = complex(u_x)
        û_y = complex(u_y)

        adv = zeros(ComplexF64, grid.nx, grid.ny)

        new( grid, rho, g, buoyancy_freq_n, odg_b, hv_order, hv_val,
             b, ψ, b̂, ψ̂, u_x, u_y, û_x, û_y, adv)

    end

end

export update_velocities!

"""
    update_velocities!(u_x, u_y, model, b̂)

Compute the streamfunction and the velocity from the
buoyancy according to the Surface-Quasi-Geostrophic model.


"""
function update_velocities!( sqg :: SQG )

    sqg.b  = real(ifft(sqg.b̂))
    sqg.ψ̂ .= sqg.b̂ ./ sqg.grid.k ./ sqg.buoyancy_freq_n
    
    sqg.û_x .= - 1im .* sqg.grid.ky' .* sqg.ψ̂
    sqg.û_y	.=  1im .* sqg.grid.kx .* sqg.ψ̂
    
    sqg.u_x .= real(ifft(sqg.û_x))
    sqg.u_y .= real(ifft(sqg.û_y))

end

include("poisson.jl")
include("init_buoyancy.jl")
include("rk4_step.jl")
include("compute_dt.jl")
include("compute_hypercoef.jl")

end

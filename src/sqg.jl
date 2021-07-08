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

    a :: Array{Float64, 2}
    â :: Array{ComplexF64, 2}

    function SQG( grid, f0 )

        rho=1e3
        g = 9.81
        buoyancy_freq_n = 3 * f0
        odg_b = 1e-3 
        hv_order = 8
        hv_val = 0

        nx, ny = grid.nx, grid.ny

        b = zeros(nx, ny)
        ψ = zeros(nx, ny)
        b̂ = zeros(ComplexF64, nx÷2+1, ny)
        ψ̂ = zeros(ComplexF64, nx÷2+1, ny)

        u_x = zeros(nx, ny)
        u_y = zeros(nx, ny)
        û_x = zeros(ComplexF64, nx÷2+1, ny)
        û_y = zeros(ComplexF64, nx÷2+1, ny)

        a = zeros(Float64, nx, ny)
        â = zeros(ComplexF64, nx÷2+1, ny)

        new( grid, rho, g, buoyancy_freq_n, odg_b, hv_order, hv_val,
             b, ψ, b̂, ψ̂, u_x, u_y, û_x, û_y, a, â)

    end

end

export update_velocities!

"""
    update_velocities!(u_x, u_y, model, b̂)

Compute the streamfunction and the velocity from the
buoyancy according to the Surface-Quasi-Geostrophic model.


"""
function update_velocities!( sqg :: SQG )

    nx, ny = sqg.grid.nx, sqg.grid.ny
    sqg.b  = irfft(sqg.b̂, nx)
    sqg.ψ̂ .= sqg.b̂ .* sqg.grid.on_k ./ sqg.buoyancy_freq_n
    
    sqg.û_x .= - 1im .* sqg.grid.ky .* sqg.ψ̂
    sqg.û_y	.=  1im .* sqg.grid.kx .* sqg.ψ̂
    
    sqg.u_x .= irfft(sqg.û_x, nx)
    sqg.u_y .= irfft(sqg.û_y, nx)

end

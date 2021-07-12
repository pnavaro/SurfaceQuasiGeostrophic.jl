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

        b̂ = zeros(ComplexF64, nx÷2+1, ny)
        ψ̂ = zeros(ComplexF64, nx÷2+1, ny)

        u_x = zeros(nx, ny)
        u_y = zeros(nx, ny)
        û_x = zeros(ComplexF64, nx÷2+1, ny)
        û_y = zeros(ComplexF64, nx÷2+1, ny)

        a = zeros(Float64, nx, ny)
        â = zeros(ComplexF64, nx÷2+1, ny)

        new( grid, rho, g, buoyancy_freq_n, odg_b, hv_order, hv_val,
             b̂, ψ̂, u_x, u_y, û_x, û_y, a, â)

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
    sqg.ψ̂ .= sqg.b̂ .* sqg.grid.on_k ./ sqg.buoyancy_freq_n
    
    sqg.û_x .= - 1im .* sqg.grid.ky .* sqg.ψ̂
    sqg.û_y	.=  1im .* sqg.grid.kx .* sqg.ψ̂
    
    sqg.u_x .= irfft(sqg.û_x, nx)
    sqg.u_y .= irfft(sqg.û_y, nx)

end


export update_advection_term!

"""
    update_advection_term(model)

Compute the Fourier transform of the partial derivative
of the advection term plus the hyperviscosity term

"""
function update_advection_term!(sqg)

    nx = sqg.grid.nx

    update_velocities!( sqg )

    # Advection term

    sqg.a .=  sqg.u_x .* irfft(-1im .* sqg.grid.kx .* sqg.b̂, nx)
    sqg.a .+= sqg.u_y .* irfft(-1im .* sqg.grid.ky .* sqg.b̂, nx)

    sqg.â .= rfft(sqg.a)

    # Summing Hyperviscosity term

    sqg.â .-= sqg.hv_val .* (sqg.grid.k.^sqg.hv_order) .* sqg.b̂

end

export TimeSolver

struct TimeSolver

    b̂ :: Array{ComplexF64, 2}
    db̂ :: Array{ComplexF64, 2}

    function TimeSolver( nx, ny)

        b̂ = zeros(ComplexF64, nx÷2+1, ny)
        db̂ = zeros(ComplexF64, nx÷2+1, ny)

        new( b̂, db̂ )

    end

end

export step!

"""
    step!(model, time_solver, dt)

Time integration of the advection equation of b 
by 4th order Runge-Kutta method
"""
function step!( sqg :: SQG, rk4 :: TimeSolver, dt )

    rk4.b̂ .= sqg.b̂

    update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt / 2
    rk4.db̂ .= sqg.â ./ 2

    update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt / 2
    rk4.db̂ .+= sqg.â

    update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt 
    rk4.db̂ .+= sqg.â

    update_advection_term!(sqg)

    rk4.db̂ .+= sqg.â ./ 2

    sqg.b̂ .= rk4.b̂ .+ dt / 3 .* rk4.db̂

end

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

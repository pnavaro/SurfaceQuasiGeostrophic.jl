"""
    deriv_fft_advection(model, b̂, w)

Compute the Fourier transform of the partial derivative
of the advection term plus the hyperviscosity term

"""
function deriv_fft_advection(model, grid, b̂, w)


    # Grid of wave vectors
    px = grid.nx ÷ 2 # frequency of aliasing
    py = grid.ny ÷ 2 # frequency of aliasing
    kx = grid.kx
    ky = grid.ky
    k2 = grid.k2

    # Advection term


    adv = zeros(ComplexF64, grid.nx, grid.ny)
    adv .= w .* (ifft(-1im .* kx  .* b̂) .+ ifft(-1im .* ky' .* b̂))

    fft!(adv)
    adv[px+1, :] .= 0 # Remove aliasing
    adv[:, py+1] .= 0 # Remove aliasing

    # Summing Hyperviscosity term

    adv .+= - model.hv_val .* k2 .^ (model.hv_order÷2) .* b̂

    return adv

end


"""
    rk4(model, fft_b, w)

Time integration of the advection equation of b 
by 4th order Runge-Kutta method with the speed w
"""
function rk4_step!(b̂, model, ux, uy)

    dt = model.dt

    k1 = deriv_fft_advection(model, grid, fft_b, w)
    k2 = deriv_fft_advection(model, grid, fft_b + k1 * dt / 2, w)
    k3 = deriv_fft_advection(model, grid, fft_b + k2 * dt / 2, w)
    k4 = deriv_fft_advection(model, grid, fft_b + k3 * dt, w)

    fft_b .+= dt / 3 .* (k1 ./ 2 .+ k2 .+ k3 .+ k4 ./ 2)

end

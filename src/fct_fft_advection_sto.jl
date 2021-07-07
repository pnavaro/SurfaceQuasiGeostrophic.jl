"""
    compute_dt( model, ogird, w)
CFL of the (large-scale) advection
"""
function compute_dt(model, grid, w)

    dx = grid.dx
    dy = grid.dy

    bound1 = sum(abs(w) * pi / dx, 3)
    bound1 = maximum(bound1)
    bound1 = 1 / bound1 / 4

    # CFL of the hyperviscosity

    dx2 = ([dx, dy] ./ pi) .^ 2
    bound2 = 1 / model.hv_val * (prod(dx2) / sum(dX2))^(model.hv_order / 2)

    # Minimum of the CFL
    minimum([bound1; bound2])

end

"""
    compute_hyper_coef( model, grid, w )

compute hyperviscosity coefficient 
"""
function compute_hyper_coef(model, grid, w)

    # compute gradient of w using 'gradient' matlab function

    # dxUx,dyUx = gradient(w(:,:,1)',grid.x_ref,grid.y_ref)
    # dxUy,dyUy = gradient(w(:,:,2)',grid.x_ref,grid.y_ref)

    dxUx = real(ifft(1im .* kx .* fft(w[1])))
    dyUx = real(ifft(1im .* ky .* fft(w[1])))
    dxUy = real(ifft(1im .* kx .* fft(w[2])))
    dyUy = real(ifft(1im .* ky .* fft(w[2])))


    s = (dxUx .+ dyUy) ./ 2
    d = (dyUx .+ dxUy) ./ 2
    lambda = sqrt.(s .^ 2 + d .^ 2)
    lambda_rms = sqrt(mean(lambda .^ 2))

    # Hyperviscosity coefficient
    coef = 40 * model / lambda_rms * (0.5 * (grid.dx + grid.dy) / pi)^model.hv_order

    return coef

end


function fct_fft_advection_sto(model, grid, fft_buoy_part)

    ## Grid [spatial & Fourier]
    model = init_grid_k(model)

    ## Initialisation of the spatial fields
    fft_w = sqg_large_uq(model, fft_buoy_part)
    w = real(ifft(fft_w))

    ## Hyperviscosity [order & coefficient]
    model.hv_order = 8
    model.hv_val = compute_hyper_coef(model, w)

    ## Choice of time step : CFL
    model.dt_adv = compute_dt(model, w)

    N_t = model.advection_duration / model.dt_adv

    # Printing some information
    println("Time step: $(model.dt_adv) seconds")
    println("Final time: $(N_t*model.dt_adv/3600/24) days ")

    ## Loop on time
    for t = 1:N_t
        ## solve "Poisson" to get w from buoy
        fft_w = sqg_large_uq(model, grid, fft_buoy_part)
        w = real(ifft(fft_w))
        ## Runge-Kutta 4 scheme
        fft_buoy_part = rk4_fft_advection(model, fft_buoy_part, w)

    end

    tt = floor(Int, Nt * model.dt_adv / (3600 * 24)) # Number of days
    println("$(t*model.dt_adv/(24*3600)) days of advection ")
    day = trunc(Int64, t * model.dt_adv / 24 / 3600)

    return fft_buoy_part

end

function fft_advection(model, grid, fft_buoy_part)

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

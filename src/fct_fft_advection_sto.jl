function fct_fft_advection_sto(model, grid, fft_buoy_part)

## Grid [spatial & Fourier]
model = init_grid_k(model)

## Initialisation of the spatial fields
fft_w = sqg_large_uq(model, fft_buoy_part)
w = real(ifft2(fft_w))

## Hyperviscosity [order & coefficient]
model.advection.HV.order=8
model.advection.HV.val=compute_hyper_coef(model, w)

## Choice of time step : CFL
model.advection.dt_adv  = compute_dt(model, w)

tt_last = -inf

N_t = model.advection.advection_duration/model.advection.dt_adv); #nb iterations

# Printing some information
println("Time step: $(model.dt_adv) seconds")
println("Final time: $(N_t*model.dt_adv/3600/24) days ")

## Loop on time
for t=1:N_t
    ## solve "Poisson" to get w from buoy
    fft_w = sqg_large_uq(model, grid, fft_buoy_part)
    w  = real(ifft2(fft_w))
    ## Runge-Kutta 4 scheme
    fft_buoy_part = rk4_fft_advection(model,fft_buoy_part, w)
    ## Plots & save()
    tt = floor(Int, t*model.dt_adv/(3600*24)) # Number of days
    if tt .> tt_last
        tt_last = tt
        println("$(t*model.dt_adv/(24*3600)) days of advection ")
        day = trunc( Int64, t * model.dt_adv/24/3600)
        # Plots
        #fct_plot(model,fft_buoy_part,day)
    end
end

fft_buoy_part

end

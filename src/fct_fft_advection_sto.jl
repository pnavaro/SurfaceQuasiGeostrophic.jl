"""
    RK4_fft_advection(model, fft_b, w)

Time integration of the advection equation of b 
by 4th order Runge-Kutta method with the speed w
"""
function fft_b_adv = RK4_fft_advection(model, fft_b, w)

dt = model.advection.dt_adv
k1 = deriv_fft_advection(model, fft_b, w)
k2 = deriv_fft_advection(model, fft_b + k1*dt/2, w)
k3 = deriv_fft_advection(model, fft_b + k2*dt/2, w)
k4 = deriv_fft_advection(model, fft_b + k3*dt, w)

fft_b_adv = fft_b + (dt/3)*(k1/2 + k2 + k3 + k4/2);

end

function fct_fft_advection_sto(model, grid, fft_buoy_part)

## Grid [spatial & Fourier]
model = init_grid_k [model]

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
fprintf(["Time step: ' num2str(model.advection.dt_adv) ' seconds \n"])
fprintf(["Final time: ' num2str(N_t*model.advection.dt_adv/3600/24) ' days \n"])

## Loop on time
for t=1:N_t
    ## solve "Poisson" to get w from buoy
    fft_w = sqg_large_uq(model, fft_buoy_part)
    w=real(ifft2(fft_w))
    ## Runge-Kutta 4 scheme
    fft_buoy_part = RK4_fft_advection(model,fft_buoy_part, w)
    ## Plots & save()
    tt = floor(t*model.advection.dt_adv/(3600*24)); # Number of days
    if tt .> tt_last
        tt_last = tt
        println("$(t*model.advection.dt_adv/(24*3600)) days of advection ")
        day = trunc( Int64, t * model.advection.dt_adv/24/3600)
        # Plots
        #fct_plot(model,fft_buoy_part,day)
    end
end

fft_buoy_part

end

"""
# Compute the streamfunction sqrt[Delta] psi = b
#and the velocity u_x = ?_y psi; u_y = ?_x psi()
"""
function sqg_large_uq(model, fft_b)


## Fourier grid()

kx=model.grid.k.kx
ky=model.grid.k.ky
on_k=model.grid.k.on_k

## Streamfunction bsxfun(@times,on_k/model.physical_constant.buoyancy_freq_N,fft_b)

## Velocity
fft_u[:,:,1,:] = -1im .* ky .*  fft_psi
fft_u[:,:,2,:] = 1im .* kx .*  fft_psi  

end

function fct_fft_advection_sto(model,  fft_buoy_part)

## Grid [spatial & Fourier]
model.grid.x = model.grid.dX[1]*(0:model.grid.MX[1]-1)
model.grid.y = model.grid.dX[2]*(0:model.grid.MX[2]-1)
model = init_grid_k [model]

## Initialisation of the spatial fields
fft_w = SQG_large_UQ[model, fft_buoy_part]
w=real(ifft2(fft_w))

## Hyperviscosity [order & coefficient]
model.advection.HV.order=8
model.advection.HV.val=compute_hyper_coef[model, w]

## Choice of time step : CFL
model.advection.dt_adv  = compute_dt[model, w]
tt_last = -inf
N_t = ceil(model.advection.advection_duration/model.advection.dt_adv); #nb iterations

# Printing some information
fprintf(["Time step: ' num2str(model.advection.dt_adv) ' seconds \n"])
fprintf(["Final time: ' num2str(N_t*model.advection.dt_adv/3600/24) ' days \n"])

## Loop on time
for t=1:N_t
    ## solve "Poisson" to get w from buoy
    fft_w = SQG_large_UQ[model, fft_buoy_part];w=real(ifft2(fft_w))
    ## Runge-Kutta 4 scheme
    fft_buoy_part = RK4_fft_advection[model,fft_buoy_part, w]
    ## Plots & save()
    tt = floor(t*model.advection.dt_adv/(3600*24)); # Number of days
    if tt .> tt_last
        tt_last = tt
        fprintf([ num2str(t*model.advection.dt_adv/(24*3600)) " days of advection \n"])
        day = num2str(floor(t*model.advection.dt_adv/24/3600))
        # Plots
        fct_plot[model,fft_buoy_part,day]
    end
end

clear fft_w s d dxUx dxUy dyUx dyUy lambda
return [fft_buoy_part, model]
end

## Initial condition for the buoyancy
[fft_buoy,model] = fct_buoyancy_init[model,resolution]

## Advection
[fft_buoy_final, model] = fct_fft_advection_sto[model, fft_buoy]
end

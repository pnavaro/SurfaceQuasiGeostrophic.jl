function compute_dt(model, w)

# CFL of the (large-scale) advection

dX=model.grid.dX;
bound1=sum(abs(w)*pi/dX(1),3);
bound1=max(bound1(:));
bound1=1/bound1/4;
# CFL of the hyperviscosity

dX2=(model.grid.dX /pi).^2;
bound2=1/model.advection.HV.val*(prod(dX2)/sum(dX2))^(model.advection.HV.order/2);

# Minimum of the CFL
dt = min([bound1 bound2]); 
clear  bound1 bound2 dX dX2

end


"""
compute hyperviscosity coefficient 
"""
function compute_hyper_coef( model, grid, w )

# compute gradient of w using 'gradient' matlab function
dxUx,dyUx = gradient(w(:,:,1)',grid.x_ref,grid.y_ref)
dxUy,dyUy = gradient(w(:,:,2)',grid.x_ref,grid.y_ref)
s=(dxUx+dyUy)/2
d=(dyUx+dxUy)/2
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));

# Hyperviscosity coefficient
coef = 40 * sqrt(mean(lambda(:).^2)) * ...
    (mean(model.grid.dX)/pi)^model.advection.HV.order;
    return coef

end


function fct_fft_advection_sto(model, grid, fft_buoy_part)

## Grid [spatial & Fourier]
model = init_grid_k(model)

## Initialisation of the spatial fields
fft_w = sqg_large_uq(model, fft_buoy_part)
w = real(ifft2(fft_w))

## Hyperviscosity [order & coefficient]
model.hv_order = 8
model.hv_val = compute_hyper_coef(model, w)

## Choice of time step : CFL
model.dt_adv  = compute_dt(model, w)

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

end

tt = floor(Int, Nt*model.dt_adv/(3600*24)) # Number of days
println("$(t*model.dt_adv/(24*3600)) days of advection ")
day = trunc( Int64, t * model.dt_adv/24/3600)

fft_buoy_part

end

"""
    rk4(model, fft_b, w)

Time integration of the advection equation of b 
by 4th order Runge-Kutta method with the speed w
"""
function fft_b_adv = rk4_fft_advection(model, fft_b, w)

dt = model.advection.dt_adv
k1 = deriv_fft_advection(model, grid, fft_b, w)
k2 = deriv_fft_advection(model, grid, fft_b + k1*dt/2, w)
k3 = deriv_fft_advection(model, grid, fft_b + k2*dt/2, w)
k4 = deriv_fft_advection(model, grid, fft_b + k3*dt, w)

fft_b_adv = fft_b .+ (dt/3)*(k1/2 + k2 + k3 + k4/2);

end


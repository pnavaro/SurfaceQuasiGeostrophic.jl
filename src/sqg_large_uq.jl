"""
    sqg_large_uq(model, fft_b)

Compute the streamfunction √Δ ψ = b
and the velocity 
   u_x = - ∂ψ / ∂y 
   u_y = ∂ψ / ∂x

"""
function sqg_large_uq(fft_b, model, grid)


## Fourier grid()

kx = grid.kx
ky = grid.ky
on_k = 1 / grid.k

## Streamfunction on_k / buoyancy_freq_N .* fft_b

## Velocity
fft_u[:,:,1,:] = -1im .* ky .*  fft_psi
fft_u[:,:,2,:] = 1im .* kx .*  fft_psi  

end

using Statistics

export compute_hypercoef

"""
    compute_hypercoef( model )

compute hyperviscosity coefficient 
"""
function compute_hypercoef(model)

    # compute gradient of w using 'gradient' matlab function

    # dxUx,dyUx = gradient(w(:,:,1)',grid.x_ref,grid.y_ref)
    # dxUy,dyUy = gradient(w(:,:,2)',grid.x_ref,grid.y_ref)

    nx, ny = model.grid.nx, model.grid.ny

    update_velocities!(model)

    dxUx = irfft(1im .* model.grid.kx .* model.û_x, nx)
    dyUx = irfft(1im .* model.grid.ky .* model.û_x, nx)
    dxUy = irfft(1im .* model.grid.kx .* model.û_y, nx)
    dyUy = irfft(1im .* model.grid.ky .* model.û_y, nx)

    s = (dxUx .+ dyUy) ./ 2
    d = (dyUx .+ dxUy) ./ 2

    @show lambda_rms = sqrt(mean(s .^ 2 .+ d .^ 2))

    dx, dy = model.grid.dx, model.grid.dy

    # Hyperviscosity coefficient
    coef = 40 * lambda_rms * ((dx+dy)/2π)^model.hv_order

    return coef

end

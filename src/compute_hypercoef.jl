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

    dxUx = real(ifft(1im .* model.grid.kx .* model.û_x))
    dyUx = real(ifft(1im .* model.grid.ky' .* model.û_x))
    dxUy = real(ifft(1im .* model.grid.kx .* model.û_y))
    dyUy = real(ifft(1im .* model.grid.ky' .* model.û_y))

    s = (dxUx .+ dyUy) ./ 2
    d = (dyUx .+ dxUy) ./ 2
    lambda = sqrt.(s .^ 2 + d .^ 2)
    lambda_rms = sqrt(mean(lambda .^ 2))

    dx, dy = model.grid.dx, model.grid.dy

    # Hyperviscosity coefficient
    coef = 40 * sqrt(mean(lambda.^2)) * ((dx+dy)/2π)^model.hv_order

    return coef

end

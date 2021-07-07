export init_buoyancy


"""
    init_buoyancy(model, resolution)

Initial condition for the buoyancy
"""
function init_buoyancy(model, grid)

    nx, ny = grid.nx, grid.ny
    x, y = grid.x, grid.y
    ee = 4
    σ = 2 * grid.lx / 15 # Length scale close to the Rossby radius

    ## Spatial buoyancy field
    # Warm anticyclones

    center1x = grid.x[nx÷4+1]
    center1y = grid.y[ny÷4+1]

    buoyancy = zeros(ComplexF64, nx, ny)

    @. buoyancy = exp(-1 / (2σ^2) * (ee * (x - center1x)^2 + (y' - center1y)^2))

    center2x = x[3nx÷4+1]
    center2y = y[ny÷4+1]

    @. buoyancy += exp(-1 / (2σ^2) .* (ee * (x .- center2x) .^ 2 + (y' .- center2y)^2))

    # Cold cyclones

    center1x = x[nx÷4+1]
    center1y = y[3ny÷4+1]

    @. buoyancy -= exp(-1 / (2σ^2) * (ee * (x - center1x)^2 + (y' - center1y)^2))

    center2x = x[3nx÷4+1]
    center2y = y[3ny÷4+1]

    @. buoyancy -= exp(-1 / (2σ^2) * (ee * (x - center2x)^2 + (y' - center2y)^2))

    # Specify the amplitude of the buoyancy
    @. buoyancy *= model.odg_b

    return buoyancy

end

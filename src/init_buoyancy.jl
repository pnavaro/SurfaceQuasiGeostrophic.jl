export init_buoyancy!


"""
    init_buoyancy!(model)

Initial condition for the buoyancy
"""
function init_buoyancy!(model :: SQG)

    nx, ny = model.grid.nx, model.grid.ny
    x, y = model.grid.x, model.grid.y
    ee = 4
    σ = 2 * model.grid.lx / 15 # Length scale close to the Rossby radius

    ## Spatial buoyancy field
    # Warm anticyclones

    center1x = model.grid.x[nx÷4+1]
    center1y = model.grid.y[ny÷4+1]

    @. model.b = exp(-1 / (2σ^2) * (ee * (x - center1x)^2 + (y' - center1y)^2))

    center2x = x[3nx÷4+1]
    center2y = y[ny÷4+1]

    @. model.b += exp(-1 / (2σ^2) .* (ee * (x .- center2x) .^ 2 + (y' .- center2y)^2))

    # Cold cyclones

    center1x = x[nx÷4+1]
    center1y = y[3ny÷4+1]

    @. model.b -= exp(-1 / (2σ^2) * (ee * (x - center1x)^2 + (y' - center1y)^2))

    center2x = x[3nx÷4+1]
    center2y = y[3ny÷4+1]

    @. model.b -= exp(-1 / (2σ^2) * (ee * (x - center2x)^2 + (y' - center2y)^2))

    # Specify the amplitude of the buoyancy
    @. model.b *= model.odg_b

    # Initialize the Fourier transform of the buyoancy
    model.b̂ .= fft(model.b)


end

export init_buoyancy!


"""
    init_buoyancy!(model)

Initial condition for the buoyancy
"""
function init_buoyancy!(model :: SQG)

    nx, ny = model.grid.nx, model.grid.ny
    lx, ly = model.grid.lx, model.grid.ly
    x, y = model.grid.x, model.grid.y
    ee = 4
    σ = 2 * model.grid.lx / 15 # Length scale close to the Rossby radius

    ## Spatial buoyancy field

    anticyclone(cx, cy) = (exp.(- ee .* (x .- cx).^2 ./ 2σ^2) 
              .* transpose(exp.(- (y .- cy).^2 ./ 2σ^2)))

    # Warm anticyclones

    b = anticyclone(lx/4, ly/4)
    b .+= anticyclone(3lx/4, ly/4)

    # Cold cyclones

    b .-= anticyclone(lx/4, 3ly/4)
    b .-= anticyclone(3lx/4, 3ly/4)

    # Specify the amplitude of the buoyancy
    b .*= model.odg_b

    # Initialize the Fourier transform of the buyoancy
    model.b̂ .= rfft(b)

end

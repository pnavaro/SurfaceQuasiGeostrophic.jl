export fct_buoyancy_init


"""
    fct_buoyancy_init(model, resolution)

Initial condition for the buoyancy
"""
function fct_buoyancy_init(model, grid)

nx, ny = grid.nx, grid.ny
x, y = grid.x, grid.y
ee = 4
σ = 2 * grid.lx/15 # Length scale close to the Rossby radius

## Spatial buoyancy field
# Warm anticyclones

center1x = grid.x[nx÷4+1]
center1y = grid.y[ny÷4+1]

b_S = zeros(nx, ny)

@. b_S = exp( -1/(2σ^2) * (ee*(x - center1x)^2 + ( y' - center1y)^2) )

center2x = x[3nx÷4+1]
center2y = y[ny÷4+1]

@. b_S += exp( -1/(2σ^2) .* (ee*(x .- center2x).^2 + (y' .- center2y)^2) )

# Cold cyclones

center1x = x[nx÷4+1]
center1y = y[3ny÷4+1]

@. b_S -= exp( -1/(2σ^2) * (ee * ( x - center1x)^2 + (y' - center1y)^2) )

center2x = x[3nx÷4+1]
center2y = y[3ny÷4+1]

@. b_S -= exp( -1/(2σ^2) * (ee * ( x - center2x)^2 + (y' - center2y)^2) )

# Specify the amplitude of the buoyancy
@. b_S *= model.odg_b

fft(b_S)

end

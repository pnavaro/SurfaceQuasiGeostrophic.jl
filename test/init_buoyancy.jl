@testset "Init Buoyancy" begin

# Resolution [even integer]
nx, ny  = 64, 64
angle_grid = π/4 
omega = 2π/(24*60*60)
f0 = 2 * omega * sin( angle_grid )

grid = Grid( 1e6, nx, 1e6, ny )

# Gather parameters in the structure model
model = SQG( grid, f0 ) 

# Initial condition for the buoyancy
init_buoyancy!(model)

@test true

end

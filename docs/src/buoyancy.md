# Buoyancy

```@example 1
using SurfaceQuasiGeostrophic, Plots

nx, ny  = 64, 64
angle_grid = π/4 
omega = 2π/(24*60*60)
f0 = 2 * omega * sin( angle_grid )

grid = Grid( 1e6, nx, 1e6, ny )

# Gather parameters in the structure model
model = ( 
    rho=1e3, # Background density
    g = 9.81, # Gravity
    buoyancy_freq_N = 3 * f0, # Background stratification
    odg_b = 1e-3 # Amplitude of the buoyancy
)

# Initial condition for the buoyancy
buoyancy = init_buoyancy(model, grid)

contourf(real(buoyancy))
```

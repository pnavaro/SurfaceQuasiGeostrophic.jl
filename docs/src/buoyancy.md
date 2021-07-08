# Buoyancy


## Intialize the buoyancy with 4 vortices

```@example 1
using SurfaceQuasiGeostrophic, Plots

nx, ny  = 64, 64
angle_grid = π/4 
omega = 2π/(24*60*60)
f0 = 2 * omega * sin( angle_grid )

grid = Grid( 1e6, nx, 1e6, ny )

sqg = SQG( grid, f0 ) 

init_buoyancy!(sqg)

contourf(sqg.b)
```

## Compute velocities by using the Surface Quasi-Geostrophic (SQG) model


```@example 1

update_velocities!( sqg )

p = plot(layout=(1,2))

contourf!(p[1,1], sqg.u_x)
contourf!(p[1,2], sqg.u_y)
display(p)

```

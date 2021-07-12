# Buoyancy

- Initialize the buoyancy with 4 vortices

```@example 1
using SurfaceQuasiGeostrophic, Plots, FFTW

nx, ny  = 64, 64
angle_grid = π/4 
omega = 2π/(24*60*60)
f0 = 2 * omega * sin( angle_grid )

grid = Grid( 1e6, nx, 1e6, ny )

sqg = SQG( grid, f0 ) 

init_buoyancy!(sqg)

contourf(irfft(sqg.b̂, nx), aspect_ratio=:equal, axis=([], false), colorbar=false)
    
```

- Compute velocities by using the Surface Quasi-Geostrophic (SQG) model


```@example 1

update_velocities!( sqg )

p = plot(layout=(1,2))

contourf!(p[1,1], sqg.u_x, aspect_ratio=:equal, axis=([], false), colorbar=false)
contourf!(p[1,2], sqg.u_y, aspect_ratio=:equal, axis=([], false), colorbar=false)

```

- Hyperviscosity [order & coefficient]
- Choice of time step : CFL

```@example 1
sqg.hv_order = 8
sqg.hv_val = compute_hypercoef(sqg)

dt = compute_dt(sqg)
println("Time step: $dt seconds")
```

- Compute the advection term `` ( u \cdot \nabla ) b `` 

```@example 1

update_advection_term!( sqg )

contourf( irfft(sqg.â, nx), aspect_ratio=:equal, axis=([], false), colorbar=false)

```

```@example 1

advection_duration = 3600*24*30

nstep = floor(Int, advection_duration / dt)

rk4 = TimeSolver(nx, ny)

anim = @animate for i = 1:nstep

    rk4.b̂ .= sqg.b̂

    update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt / 2
    rk4.db̂ .= sqg.â ./ 2

    update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt / 2
    rk4.db̂ .+= sqg.â

    update_advection_term!(sqg)

    sqg.b̂ .= rk4.b̂ .+ sqg.â .* dt 
    rk4.db̂ .+= sqg.â

    update_advection_term!(sqg)

    rk4.db̂ .+= sqg.â ./ 2

    sqg.b̂ .= rk4.b̂ .+ dt / 3 .* rk4.db̂

    contourf(irfft(sqg.b̂, nx), linewidth=0, aspect_ratio=:equal, axis=([], false), 
        colorbar=false, clims=(-sqg.odg_b,sqg.odg_b) )

end every 50 

gif(anim, "assets/buoyancy.gif", fps = 10) # hide
nothing # hide

```

![buoyancy](assets/buoyancy.gif)

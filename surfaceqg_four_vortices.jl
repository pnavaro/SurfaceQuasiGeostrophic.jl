# -*- coding: utf-8 -*-
using GeophysicalFlows, Plots, Printf, Random

 dev = CPU() 

n = 256
stepper = "FilteredETDRK4"
dt = 0.03
tf =60
nsteps = Int(tf / dt)
nsubs = round(Int, nsteps/100)

L = 2π
ν = 1e-19
nν = 4

prob = SurfaceQG.Problem(dev; nx=n, Lx=L, dt=dt, 
    stepper=stepper, ν = ν , nν = nν )

sol, clock, vars, params, grid = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
x, y = grid.x, grid.y

X, Y = gridpoints(grid)
b = @. exp(-((X+L/4)^2 + 4*(Y+L/4)^2))
b .+= @. exp(-((X+L/4)^2 + 4*(Y-L/4)^2))
b .-= @. exp(-((X-L/4)^2 + 4*(Y+L/4)^2))
b .-= @. exp(-((X-L/4)^2 + 4*(Y-L/4)^2))

SurfaceQG.set_b!(prob, b )

heatmap(x, y, vars.b', aspectratio = 1,
               c = :deep,
            clim = (-1, 1),
           xlims = (-grid.Lx/2, grid.Lx/2),
           ylims = (-grid.Ly/2, grid.Ly/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
ylabel = "y",
title = "buoyancy b ",
      framestyle = :box)

B  = Diagnostic(SurfaceQG.buoyancy_variance, prob; nsteps=nsteps)
KE = Diagnostic(SurfaceQG.kinetic_energy, prob; nsteps=nsteps)
Dᵇ = Diagnostic(SurfaceQG.buoyancy_dissipation, prob; nsteps=nsteps)
diags = [B, KE, Dᵇ] # A list of Diagnostics types passed to `stepforward!`. Diagnostics are updated every timestep.

# +
base_filename = string("four_vortices_", n)

datapath = "./"
plotpath = "./"

dataname = joinpath(datapath, base_filename) 
plotname = joinpath(plotpath, base_filename)

if !isdir(plotpath); mkdir(plotpath); end 
if !isdir(datapath); mkdir(datapath); end
# -

get_sol(prob) = sol 
get_u(prob) = irfft(im * grid.l .* sqrt.(grid.invKrsq) .* sol, grid.nx) 
out = Output(prob, dataname, (:sol, get_sol), (:u, get_u))

function plot_output(prob)
  b = prob.vars.b

  pb = heatmap(x, y, Array(b'),
       aspectratio = 1,
                 c = :deep,
              clim = (-1, 1),
             xlims = (-grid.Lx/2, grid.Lx/2),
             ylims = (-grid.Ly/2, grid.Ly/2),
            xticks = -3:3,
            yticks = -3:3,
            xlabel = "x",
            ylabel = "y",
             title = "buoyancy bₛ",
        framestyle = :box)

  pKE = plot(1,
             label = "kinetic energy ∫½(uₛ²+vₛ²)dxdy/L²",
         linewidth = 2,
            legend = :bottomright,
             alpha = 0.7,
             xlims = (0, tf),
             ylims = (0, 1e-2),
            xlabel = "t")

  pb² = plot(1,
             label = "buoyancy variance ∫bₛ²dxdy/L²",
         linecolor = :red,
            legend = :bottomright,
         linewidth = 2,
             alpha = 0.7,
             xlims = (0, tf),
             ylims = (0, 2e-2),
            xlabel = "t")

  layout = @layout [a{0.5w} Plots.grid(2, 1)]
  p = plot(pb, pKE, pb², layout=layout, size = (900, 500))

  return p
end

startwalltime = time()

p = plot_output(prob)

anim = @animate for j = 0:round(Int, nsteps/nsubs)
  if j % (500 / nsubs) == 0
    cfl = clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])

    log1 = @sprintf("step: %04d, t: %.1f, cfl: %.3f, walltime: %.2f min",
          clock.step, clock.t, cfl, (time()-startwalltime)/60)

    log2 = @sprintf("buoyancy variance: %.2e, buoyancy variance dissipation: %.2e",
              B.data[B.i], Dᵇ.data[Dᵇ.i])

    println(log1)

    println(log2)
  end

  p[1][1][:z] = vars.b
  p[1][:title] = "buoyancy, t=" * @sprintf("%.2f", clock.t)
  push!(p[2][1], KE.t[KE.i], KE.data[KE.i])
  push!(p[3][1], B.t[B.i], B.data[B.i])

  stepforward!(prob, diags, nsubs)
  SurfaceQG.updatevars!(prob)
end

mp4(anim, "four_vortices.mp4", fps=14)

pu = heatmap(x, y, Array(vars.u'),
     aspectratio = 1,
               c = :balance,
            clim = (-maximum(abs.(vars.u)), maximum(abs.(vars.u))),
           xlims = (-L/2, L/2),
           ylims = (-L/2, L/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
          ylabel = "y",
           title = "uₛ(x, y, t=" * @sprintf("%.2f", clock.t) * ")",
      framestyle = :box)

pv = heatmap(x, y, Array(vars.v'),
     aspectratio = 1,
               c = :balance,
            clim = (-maximum(abs.(vars.v)), maximum(abs.(vars.v))),
           xlims = (-L/2, L/2),
           ylims = (-L/2, L/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
          ylabel = "y",
           title = "vₛ(x, y, t=" * @sprintf("%.2f", clock.t) * ")",
      framestyle = :box)

pb = heatmap(x, y, Array(vars.b'),
     aspectratio = 1,
               c = :deep,
            clim = (-1, 1),
           xlims = (-L/2, L/2),
           ylims = (-L/2, L/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
          ylabel = "y",
           title = "bₛ(x, y, t=" * @sprintf("%.2f", clock.t) * ")",
      framestyle = :box)

layout = @layout [a{0.5h}; b{0.5w} c{0.5w}]

plot_final = plot(pb, pu, pv, layout=layout, size = (800, 800))



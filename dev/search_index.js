var documenterSearchIndex = {"docs":
[{"location":"buoyancy/#Buoyancy","page":"Buoyancy","title":"Buoyancy","text":"","category":"section"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"Initialize the buoyancy with 4 vortices","category":"page"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"using SurfaceQuasiGeostrophic, Plots, FFTW\n\nnx, ny  = 64, 64\nangle_grid = π/4 \nomega = 2π/(24*60*60)\nf0 = 2 * omega * sin( angle_grid )\n\ngrid = Grid( 1e6, nx, 1e6, ny )\n\nsqg = SQG( grid, f0 ) \n\ninit_buoyancy!(sqg)\n\ncontourf(sqg.b, aspect_ratio=:equal, axis=([], false), colorbar=false)","category":"page"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"Compute velocities by using the Surface Quasi-Geostrophic (SQG) model","category":"page"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"\nupdate_velocities!( sqg )\n\np = plot(layout=(1,2))\n\ncontourf!(p[1,1], sqg.u_x, aspect_ratio=:equal, axis=([], false), colorbar=false)\ncontourf!(p[1,2], sqg.u_y, aspect_ratio=:equal, axis=([], false), colorbar=false)\n","category":"page"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"Hyperviscosity [order & coefficient]\nChoice of time step : CFL","category":"page"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"sqg.hv_order = 8\nsqg.hv_val = compute_hypercoef(sqg)\n\ndt = compute_dt(sqg)\nprintln(\"Time step: $dt seconds\")","category":"page"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"Compute the advection term ( u cdot nabla ) b ","category":"page"},{"location":"buoyancy/","page":"Buoyancy","title":"Buoyancy","text":"\nupdate_advection_term!( sqg )\n\ncontourf( irfft(sqg.â, nx), aspect_ratio=:equal, axis=([], false), colorbar=false)\n","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SurfaceQuasiGeostrophic","category":"page"},{"location":"#SurfaceQuasiGeostrophic","page":"Home","title":"SurfaceQuasiGeostrophic","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SurfaceQuasiGeostrophic.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The unknown of the model is the buoyancy b(txy) and the stream function psi(txy)  which are related through a Poisson like equation  i nabla  psi = b, which becomes in Fourier k hatpsi = hatb","category":"page"},{"location":"","page":"Home","title":"Home","text":"The dynamical model is the following","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracpartial bpartial t + (nabla times psi) cdot  nabla b = 0 ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The solver uses Fourier techniques so that","category":"page"},{"location":"","page":"Home","title":"Home","text":"nabla times psi cdot nabla b = ifft(ik) times hatpsi cdot ifft(ik)hatb","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SurfaceQuasiGeostrophic]","category":"page"},{"location":"#SurfaceQuasiGeostrophic.Grid","page":"Home","title":"SurfaceQuasiGeostrophic.Grid","text":"init_grid_k(model)\n\nCreate a grid in the Fourier space\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceQuasiGeostrophic.Poisson","page":"Home","title":"SurfaceQuasiGeostrophic.Poisson","text":"Poisson( grid )\n\nDerived type to solve the Poisson equation on 2d regular  cartesian mesh with periodic boundary conditions on both sides\n\nkx   : wave number in x\nky   : wave number in y\nk2   : k_x^2 + k_y^2\nnc_x : cells number in x\nnc_y : cells number in y\ndx   : x step size\ndy   : y step size\nrht  : fft(rho)\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceQuasiGeostrophic.SQG","page":"Home","title":"SurfaceQuasiGeostrophic.SQG","text":"SQG\n\nData structure to store arrays and parameters for the Surface-Quasi-Geostrophic model\n\nrho : Background density\ng : Gravity\nbuoyancyfreqn : Background stratification\nodg_b : Amplitude of the buoyancy\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceQuasiGeostrophic.compute_dt-Tuple{Any}","page":"Home","title":"SurfaceQuasiGeostrophic.compute_dt","text":"compute_dt( model, ogird, w)\n\nCFL of the (large-scale) advection\n\n%% Choice of time step : CFL\n\n% CFL of the diffusion (or CFL of the white noise advection) dX2=(model.grid.dX /pi).^2; bound1=2/model.sigma.a0*prod(dX2)/sum(dX2);\n\n% CFL of the (large-scale) advection dX=permute(model.grid.dX,[1 3 2]); bound2=sum(bsxfun(@times,abs(w),pi./dX),3); bound2=max(bound2(:)); bound2=1/bound2/4;\n\n% CFL of the hyperviscosity bound3=1/model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...     (model.advection.HV.order/2);\n\ndt = min([bound1 bound2 bound3]);\n\n\n\n\n\n","category":"method"},{"location":"#SurfaceQuasiGeostrophic.compute_hypercoef-Tuple{Any}","page":"Home","title":"SurfaceQuasiGeostrophic.compute_hypercoef","text":"compute_hypercoef( model )\n\ncompute hyperviscosity coefficient \n\n\n\n\n\n","category":"method"},{"location":"#SurfaceQuasiGeostrophic.compute_psi!-Tuple{Any, Poisson, Any}","page":"Home","title":"SurfaceQuasiGeostrophic.compute_psi!","text":"compute_psi!( poisson, psi, omega )\n\ncomputes psi from omega\n\n-Delta psi(xy) = omega(xy)\n\n\n\n\n\n","category":"method"},{"location":"#SurfaceQuasiGeostrophic.compute_velocities!-Tuple{Any, Any, Poisson, Any}","page":"Home","title":"SurfaceQuasiGeostrophic.compute_velocities!","text":"compute_velocities!( ux, uy, poisson, omega )\n\nsolves Poisson equation to compute velocity fields\n\nDelta psi(xy) = - omega(xy) \nu(xy) = nabla times psi(xy) \n\n\n\n\n\n","category":"method"},{"location":"#SurfaceQuasiGeostrophic.init_buoyancy!-Tuple{SQG}","page":"Home","title":"SurfaceQuasiGeostrophic.init_buoyancy!","text":"init_buoyancy!(model)\n\nInitial condition for the buoyancy\n\n\n\n\n\n","category":"method"},{"location":"#SurfaceQuasiGeostrophic.step!-Tuple{SQG, TimeSolver, Any}","page":"Home","title":"SurfaceQuasiGeostrophic.step!","text":"step!(model, time_solver, dt)\n\nTime integration of the advection equation of b  by 4th order Runge-Kutta method\n\n\n\n\n\n","category":"method"},{"location":"#SurfaceQuasiGeostrophic.update_advection_term!-Tuple{Any}","page":"Home","title":"SurfaceQuasiGeostrophic.update_advection_term!","text":"update_advection_term(model)\n\nCompute the Fourier transform of the partial derivative of the advection term plus the hyperviscosity term\n\n\n\n\n\n","category":"method"},{"location":"#SurfaceQuasiGeostrophic.update_velocities!-Tuple{SQG}","page":"Home","title":"SurfaceQuasiGeostrophic.update_velocities!","text":"update_velocities!(u_x, u_y, model, b̂)\n\nCompute the streamfunction and the velocity from the buoyancy according to the Surface-Quasi-Geostrophic model.\n\n\n\n\n\n","category":"method"}]
}

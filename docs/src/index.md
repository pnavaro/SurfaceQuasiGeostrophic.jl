```@meta
CurrentModule = SurfaceQuasiGeostrophic
```

# SurfaceQuasiGeostrophic

Documentation for [SurfaceQuasiGeostrophic](https://github.com/pnavaro/SurfaceQuasiGeostrophic.jl).

The unknown of the model is the buoyancy ``b(t,x,y)`` and the stream function ``\psi(t,x,y)`` 
which are related through a Poisson like equation
``| −i \nabla | \psi = b``, which becomes in Fourier ``|k| \hat{\psi} = \hat{b}``

The dynamical model is the following

```math
\frac{\partial b}{\partial t} + (\nabla \times \psi) \cdot  \nabla b = 0. 
```

The solver uses Fourier techniques so that

```math
\nabla \times \psi \cdot \nabla b = ifft[(ik) \times \hat{\psi}] \cdot ifft[(ik)\hat{b}],
```

```@index
```

```@autodocs
Modules = [SurfaceQuasiGeostrophic]
```

using FFTW

export Poisson

"""
    Poisson( grid )

Derived type to solve the Poisson equation on 2d regular 
cartesian mesh with periodic boundary conditions on both sides
- kx   : wave number in x
- ky   : wave number in y
- k2   : ``k_x^2 + k_y^2``
- nc_x : cells number in x
- nc_y : cells number in y
- dx   : x step size
- dy   : y step size
- rht  : fft(rho)
"""
struct Poisson

    grid::Grid
    kx::Array{Float64,2}
    ky::Array{Float64,2}
    k2::Array{Float64,2}
    rht::Array{ComplexF64,2}

    function Poisson(grid::Grid)

        nx = grid.nx
        ny = grid.ny

        rht = zeros(ComplexF64, (div(nx, 2) + 1, ny))

        kx = zeros(nx ÷ 2 + 1, ny)
        ky = zeros(nx ÷ 2 + 1, ny)
        k2 = zeros(nx ÷ 2 + 1, ny)

        kx0 = 2π / grid.lx
        ky0 = 2π / grid.ly

        for ik = 1:nx÷2+1
            kx1 = (ik - 1) * kx0
            for jk = 1:ny÷2
                kx[ik, jk] = kx1
                ky[ik, jk] = (jk - 1) * ky0
            end
            for jk = ny÷2+1:ny
                kx[ik, jk] = kx1
                ky[ik, jk] = (jk - 1 - ny) * ky0
            end
        end

        kx[1, 1] = 1.0
        k2 .= kx .* kx .+ ky .* ky
        kx .= kx ./ k2
        ky .= ky ./ k2

        new(grid, kx, ky, k2, rht)

    end

end


export compute_psi!

"""
    compute_psi!( poisson, psi, omega )

computes ``\\psi`` from ``\\omega``

```math
-\\Delta \\psi(x,y) = \\omega(x,y)
```

"""
function compute_psi!(psi, poisson::Poisson, omega)

    poisson.rht .= rfft(omega)
    poisson.rht ./= poisson.k2
    psi .= irfft(poisson.rht, poisson.grid.nx)

end

export compute_velocities!

"""
    compute_velocities!( ux, uy, poisson, omega )

solves Poisson equation to compute velocity fields

```math
\\Delta \\psi(x,y) = - \\omega(x,y) \\\\
u(x,y) = \\nabla \\times \\psi(x,y) 
```

"""
function compute_velocities!(ux, uy, poisson::Poisson, omega)

    poisson.rht .= rfft(omega)

    poisson.rht[1, 1] = 0.0
    poisson.rht .*= -1im .* poisson.ky

    ux .= irfft(poisson.rht, poisson.grid.nx)

    poisson.rht .= rfft(omega)

    poisson.rht[1, 1] = 0.0
    poisson.rht .*= 1im .* poisson.kx

    uy .= irfft(poisson.rht, poisson.grid.nx)

end

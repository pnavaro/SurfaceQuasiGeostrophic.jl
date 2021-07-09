function fct_unity_approx5(n)

    @assert n % 2 == 0

    slop_size_ratio = 6
    t = ones(n)
    p = n ÷ 2
    s = ceil(Int, n / slop_size_ratio)

    t[(p-s+1):p] .= (-tanh.(-3 .+ 6 / (s - 1) .* collect(0:(s-1))) .+ 1) ./ 2
    t[(p+2):(p+1+s)] .= (tanh.(-3 .+ 6 / (s - 1) .* collect(0:(s-1))) .+ 1) ./ 2

    t[p+1] = 0

    return t

end

export Grid

"""
    init_grid_k(model)

Create a grid in the Fourier space

"""
struct Grid

    nx::Int
    ny::Int
    lx::Float64
    ly::Float64
    dx::Float64
    dy::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    k::Array{Float64,2}
    on_k::Array{Float64,2}

    function Grid(lx, nx, ly, ny; meth_anti_alias = :deriv_LowPass)

        dx = lx / nx
        dy = ly / ny

        kx = zeros(nx ÷ 2 + 1, ny)
        ky = zeros(nx ÷ 2 + 1, ny)
        k2 = zeros(nx ÷ 2 + 1, ny)
        k  = zeros(nx ÷ 2 + 1, ny)
        on_k = zeros(nx ÷ 2 + 1, ny)

        kx0 = 2π / lx
        ky0 = 2π / ly

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

        k2 .= kx .* kx .+ ky .* ky
        k .= sqrt.(k2)

        on_k[ k .!= 0] .= 1 ./ k[ k .!= 0] 
        
        x = LinRange(0, lx, nx + 1)[1:end-1]
        y = LinRange(0, ly, ny + 1)[1:end-1]

        new(nx, ny, lx, ly, dx, dy, x, y, k, on_k)

    end

end

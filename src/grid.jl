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
    kx::Vector{Float64}
    ky::Vector{Float64}
    k::Array{Float64,2}

    function Grid(lx, nx, ly, ny; meth_anti_alias = :deriv_LowPass)

        px = nx ÷ 2
        py = ny ÷ 2

        dx = lx / nx
        dy = ly / ny

        kx = Float64[0:px-1; -px:-1]
        ky = Float64[0:py-1; -py:-1]

        if meth_anti_alias == :deriv_LowPass
            kx .*= fct_unity_approx5(nx)
            ky .*= fct_unity_approx5(ny)
        end

        kx .*= 2π / lx
        ky .*= 2π / ly

        k2 = kx .^ 2 .+ transpose(ky) .^ 2

        k2[ k2 .== 0] .= 1

        x = LinRange(0, lx, nx + 1)[1:end-1]
        y = LinRange(0, ly, ny + 1)[1:end-1]

        new(nx, ny, lx, ly, dx, dy, x, y, kx, ky, sqrt.(k2))

    end

end

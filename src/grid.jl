function fct_unity_approx5(n)

    @assert n % 2 == 0

    slop_size_ratio = 6
    t = ones(n)
    p = n ÷ 2
    s = ceil(Int, n/slop_size_ratio)

    t[(p-s+1):p] .= (-tanh.(-3 .+ 6/(s-1) .* collect(0:(s-1))) .+1)./2
    t[(p+2):(p+1+s)] .= (tanh.(-3 .+ 6/(s-1) .* collect(0:(s-1))).+1)./2

    t[p+1] = 0

    return t

end

export Grid

"""
    init_grid_k(model)

Create a grid in the Fourier space

"""
struct Grid

    nx :: Int
    ny :: Int
    lx :: Float64
    ly :: Float64
    dx :: Float64
    dy :: Float64
    x :: Vector{Float64}
    y :: Vector{Float64}
    kx :: Vector{Float64}
    ky :: Vector{Float64}
    k2 :: Array{Float64, 2}

    function Grid( lx, nx, ly, ny; meth_anti_alias = :deriv_LowPass)

        px = nx ÷ 2
        py = ny ÷ 2

        dx = lx / nx
        dy = ly / ny

        kx = zeros(nx ÷ 2 + 1, ny)
        ky = zeros(ny ÷ 2 + 1, ny)
        k2 = zeros(nx ÷ 2 + 1, ny)

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

        kx[1, 1] = 1.0
        k2 .= kx .* kx .+ ky .* ky
        kx .= kx ./ k2
        ky .= ky ./ k2


        kx = Float64[ 0:px-1; -px:-1] 
        ky = Float64[ 0:py-1; -py:-1]

        if meth_anti_alias == :deriv_LowPass
            kx .*= fct_unity_approx5(nx)
            ky .*= fct_unity_approx5(ny)
        end
        
        kx .*= 2π/lx
        ky .*= 2π/ly
        
        k2 = kx.^2 .+ transpose(ky).^2
        k2[px+1,:] .= 0
        k2[:,py+1] .= 0

        x = LinRange(0, lx, nx+1)[1:end-1]
        y = LinRange(0, ly, ny+1)[1:end-1]

        new( nx, ny, lx, ly, dx, dy, x, y, kx, ky, k2 )

    end
        
end


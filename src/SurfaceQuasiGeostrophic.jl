module SurfaceQuasiGeostrophic

using FFTW

include("grid.jl")
include("poisson.jl")
include("fct_buoyancy_init.jl")
include("sqg_large_uq.jl")
include("deriv_fft_advection.jl")
include("rk4_fft_advection.jl")
include("fct_fft_advection_sto.jl")

end

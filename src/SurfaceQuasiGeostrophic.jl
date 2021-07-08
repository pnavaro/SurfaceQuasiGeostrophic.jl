module SurfaceQuasiGeostrophic

using FFTW

include("grid.jl")
include("sqg.jl")
include("poisson.jl")
include("init_buoyancy.jl")
include("rk4_step.jl")
include("compute_dt.jl")
include("compute_hypercoef.jl")

end

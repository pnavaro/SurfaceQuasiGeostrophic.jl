export Model

"""
    Model

Data structure to store some parameters about the simulation

- rho : Background density
- g : Gravity
- buoyancy_freq_n : Background stratification
- odg_b : Amplitude of the buoyancy
"""
mutable struct Model

    rho :: Float64
    g :: Float64
    buoyancy_freq_n :: Float64
    odg_b :: Float64
    hv_order :: Int
    hv_val :: Float64

    function Model( f0 )

    rho=1e3
    g = 9.81
    buoyancy_freq_n = 3 * f0
    odg_b = 1e-3 
    hv_order = 8
    hv_val = 0

    new( rho, g, buoyancy_freq_n, odg_b, hv_order, hv_val)

    end



end

```@meta
CurrentModule = SurfaceQuasiGeostrophic
```

# SurfaceQuasiGeostrophic

Documentation for [SurfaceQuasiGeostrophic](https://github.com/pnavaro/SurfaceQuasiGeostrophic.jl).


dans le main 
paramètres divers  
initialisation (fct_buoyancy_init.m) 
appel au solver fct_fft_advection_sto  

dans fct_fft_advection_sto 
définition des maillages et modes
initialisation du Poisson (avec sqrt !)
paramètres (pour le terme hyper-visqueux, dt)
boucle en temps : appel à Poisson puis RK4 

dans RK4_fft_advection
calcul des étapes de RK4 avec appel à derive_fft_advection 

dans derive_fft_advection 
calcul des termes d'advection u_x ∂_x b et u_y ∂_y b en pseudo-spectral (fft) 
calcul du terme hyper-visqueux C (∂_x^8 + ∂_y^8) b en Fourier (--> C |k|^8 F(u)) 

```@index
```

```@autodocs
Modules = [SurfaceQuasiGeostrophic]
```

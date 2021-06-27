"""
    deriv_fft_advection(model, fft_b, w)

Compute the Fourier transform of the partial derivative
of the advection term (adv1) and hyperviscosity (adv2) term

"""
function deriv_fft_advection(model, grid, fft_b, w)


# Grid of wave vectors
px = grid.nx ÷ 2 # frequency of aliasing
py = grid.ny ÷ 2 # frequency of aliasing
kx = grid.kx
ky = grid.ky
k2 = grid.k2

# Advection term
gradb = [similar(fft_b), similar(fft_b)]
adv1 = similar(fft_b)
adv2 = similar(fft_b)

gradb[1] .= real(ifft(- 1im * kx .* fft_b))
gradb[2] .= real(ifft(- 1im * ky .* fft_b))

wgradT = w .* ( gradb[1] .+ gradb[2] )

adv1 = fft(wgradT)
adv1[px+1,:] .= 0 # Remove aliasing
adv1[:,py+1] .= 0 # Remove aliasing

# Hyperviscosity

# PN # adv2 = - model.advection.HV.val * k2 .^ (model.advection.HV.order/2) .* fft_b

# Summing terms
return adv1 .+ adv2

end



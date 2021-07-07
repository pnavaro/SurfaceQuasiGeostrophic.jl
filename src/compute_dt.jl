export compute_dt

"""
    compute_dt( model, ogird, w)

CFL of the (large-scale) advection
"""
function compute_dt(model)

    dx = model.grid.dx
    dy = model.grid.dy

    bound1 = sum(abs.(model.u_x) * pi / dx .+ abs.(model.u_y) * pi / dy )
    bound1 = maximum(bound1)
    bound1 = 1 / bound1 / 4

    # CFL of the hyperviscosity

    dx2 = ((dx, dy) ./ pi) .^ 2
    bound2 = 1 / model.hv_val * (prod(dx2) / sum(dx2))^(model.hv_order ÷ 2)

    # Minimum of the CFL
    min(minimum(bound1), bound2)

end


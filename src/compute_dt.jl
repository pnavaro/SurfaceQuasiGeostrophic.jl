export compute_dt

"""
    compute_dt( sqg, ogird, w)

CFL of the (large-scale) advection


"""
function compute_dt(sqg)

    update_velocities!(sqg)
    dx = sqg.grid.dx
    dy = sqg.grid.dy

    bound1 = maximum(abs.(sqg.u_x .* π ./ dx) .+ abs.(sqg.u_y .* π ./ dy))
    bound1 = 1 / bound1 / 4

    # CFL of the hyperviscosity

    dx2 = ((dx, dy) ./ pi) .^ 2
    bound2 = 1 / sqg.hv_val * (prod(dx2) / sum(dx2))^(sqg.hv_order ÷ 2)

    # Minimum of the CFL
    min(minimum(bound1), bound2)

end


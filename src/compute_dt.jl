export compute_dt

"""
    compute_dt( model, ogird, w)

CFL of the (large-scale) advection


%% Choice of time step : CFL

% CFL of the diffusion (or CFL of the white noise advection)
dX2=(model.grid.dX /pi).^2;
bound1=2/model.sigma.a0*prod(dX2)/sum(dX2);

% CFL of the (large-scale) advection
dX=permute(model.grid.dX,[1 3 2]);
bound2=sum(bsxfun(@times,abs(w),pi./dX),3);
bound2=max(bound2(:));
bound2=1/bound2/4;

% CFL of the hyperviscosity
bound3=1/model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...
    (model.advection.HV.order/2);

dt = min([bound1 bound2 bound3]);

"""
function compute_dt(model)

    dx = model.grid.dx
    dy = model.grid.dy

    bound1 = sum(abs.(model.u_x) .* pi ./ dx .+ abs.(model.u_y) .* pi ./ dy )
    bound1 = maximum(bound1)
    @show bound1 = 1 / bound1 / 4

    # CFL of the hyperviscosity

    dx2 = ((dx, dy) ./ pi) .^ 2
    @show bound2 = 1 / model.hv_val * (prod(dx2) / sum(dx2))^(model.hv_order ÷ 2)

    # Minimum of the CFL
    min(minimum(bound1), bound2)

end


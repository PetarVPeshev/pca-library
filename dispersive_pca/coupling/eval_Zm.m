function Zm = eval_Zm(f, slot, x, args)
%EVAL_ZM Summary of this function goes here
%   Detailed explanation goes here
    arguments
        f       (1,:) double             {mustBePositive}
        slot          SlotInDielectrics
        x       (1,1) double             {mustBeReal}
        args.dx (1,1) double             {mustBePositiveOrNan} = NaN
    end

    dx = args.dx;
    if isnan(dx)
        dx = slot.d_gap;
    end
    
    c0 = get_phys_const('LightSpeed');

    Nf     = length(f);
    k0     = 2 * pi * f / c0;
    kx_max = slot.get_integration_domain(f);

    Dfs    = @(kx, kfs) slot.compute_D(kx, kfs, 'Top', 'Top');
    feed_1 = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    feed_2 = @(kx) sinc(kx * dx / (2 * pi));
    
    Zm = NaN(1, Nf);
    parfor idx = 1 : Nf
        start_pt = - (kx_max(idx) * k0(idx) + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        D         = @(kx) Dfs(kx, k0(idx));
        integrand = @(kx) (feed_1(kx) .* feed_2(kx)) .* exp(- 1j * kx * x) ./ D(kx);

        Zm(idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts);
    end

    Zm = Zm / (2 * pi);
end


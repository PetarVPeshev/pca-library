function Zm = eval_Zmp(f, slot, x, args)
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

    Nf  = length(f);
    k0  = 2 * pi * f / c0;
    kxp = slot.find_kxp(f);

    Dpfs = @(kxp, kfs) slot.compute_Dp(kxp, kfs, 'Top', 'Bottom');
    
    Dp   = NaN(1, Nf);
    parfor idx = 1 : Nf
        Dp(idx) = Dpfs(kxp(idx), k0(idx));
    end
    
    feed_1 = sinc(kxp * slot.d_gap / (2 * pi));
    feed_2 = sinc(kxp * dx / (2 * pi));

    Zm = - 1j * feed_1 .* feed_2 .* exp(- 1j * kxp * abs(x)) ./ Dp;
end


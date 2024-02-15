close all;
clear;
clc;

env;
sim_params;

COLORS = ["#0072BD", "#A2142F", "#EDB120"];

%% PARAMETERS
f = linspace(eps, 2, 4001) * 1e12;
% FEED
d_gap = 4.5 * 1e-6;
d     = (100 : 100 : 300) * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% NUMBER OF POINTS
Nf = length(f);
Nd = length(d);

%% SLOT, PHOTO-CONDUCTOR, LASER, AND TIME-STEP
slot = create_slot(d_gap, ws, params_slot);

%% MUTUAL IMPEDANCE
Zm = NaN(Nd, Nf);
for idx = 1 : Nd
    Zm(idx, :) = eval_Zm(f, slot, d(idx));
end

% REAL PART PLOT
figure('Position', [250 250 750 450]);

for idx = 1 : Nd
    plot(f * 1e-12, real(Zm(idx, :)), 'LineWidth', 1.5, 'Color', COLORS(idx), ...
         'DisplayName', ['d = ' num2str(d(idx) * 1e6) ' \mum']);
    hold on;
end
hold off;

grid on;
legend('location', 'bestoutside');
    
xlabel('f [THz]');
ylabel('\Re\{ Z_{m} \} [\Omega]');
title(['\Re @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum']);

% IMAG PART PLOT
figure('Position', [250 250 750 450]);

for idx = 1 : Nd
    plot(f * 1e-12, imag(Zm(idx, :)), 'LineWidth', 1.5, 'Color', COLORS(idx), ...
         'DisplayName', ['d = ' num2str(d(idx) * 1e6) ' \mum']);
    hold on;
end
hold off;

grid on;
legend('location', 'bestoutside');
    
xlabel('f [THz]');
ylabel('\Im\{ Z_{m} \} [\Omega]');
title(['\Im @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum']);

%% FUNCTIONS
function slot = create_slot(d_gap, ws, params_slot)
    slot = SlotInDielectrics(d_gap, ws, params_slot.er_up, params_slot.er_dn);
end

function Zm = eval_Zm(f, slot, d)
    c0 = get_phys_const('LightSpeed');

    Nf     = length(f);
    k0     = 2 * pi * f / c0;
    kx_max = slot.get_integration_domain(f);

    Dfs  = @(kx, kfs) slot.compute_D(kx, kfs, 'Top', 'Top');
    feed = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    
    Zm = NaN(1, Nf);
    parfor idx = 1 : Nf
        start_pt = - (kx_max(idx) * k0(idx) + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        D         = @(kx) Dfs(kx, k0(idx));
        integrand = @(kx) (feed(kx) .^ 2) .* exp(- 1j * kx * d) ./ D(kx);

        Zm(idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts);
    end

    Zm = Zm / (2 * pi);
end

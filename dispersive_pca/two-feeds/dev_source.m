close all;
clear;
clc;

env;
sim_params;

COLORS = ["#0072BD", "#A2142F"];

%% PARAMETERS
f = linspace(eps, 2e12, 1001);
% FEED
d_gap = 4.5 * 1e-6;
d     = 200 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;

%% SLOT
slot = create_slot(d_gap, ws, params_slot);

%% INPUT IMPEDANCE
tic;
Zin = slot.compute_zin(f);
toc

% PLOT - PLOT IS CORRECT
% figure('Position', [250 250 750 450]);
% 
% plot(f * 1e-12, real(Zin), 'LineWidth', 1.5, 'DisplayName', '\Re');
% hold on;
% plot(f * 1e-12, imag(Zin), '--', 'LineWidth', 1.5, 'DisplayName', '\Im');
% 
% grid on;
% legend('location', 'bestoutside');
% ylim([0 80]);
% 
% xlabel('f [THz]');
% ylabel('Z_{in} [\Omega]');
% title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum']);

%% MUTUAL IMPEDANCE
tic;
Znm = eval_Znm(f, slot, d);
toc

% PLOT - PLOT LOOKS CORRECT
% fig = figure('Position', [250 250 750 450]);
% 
% AX1 = axes(fig);
% L = gobjects(1, 2);
% for idx = 1 : 2
%     L(idx) = plot(f * 1e-12, real(Znm{1, idx}), 'LineWidth', 1.5, 'Color', COLORS(idx), 'Parent', AX1);
%     hold on;
%              plot(f * 1e-12, imag(Znm{1, idx}), '--', 'LineWidth', 1.5, 'Color', COLORS(idx), 'Parent', AX1);
%     hold on;
% end
% 
% LGND1 = legend(AX1, L, {'Z_{11}', 'Z_{12}'}, 'location', 'bestoutside');
% hold on;
% grid on;
%     
% AX2 = copyobj(AX1, gcf);
% delete( get(AX2, 'Children') );
% H1 = plot(NaN, 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
% H2 = plot(NaN, '--', 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
% 
% set(AX2, 'Color', 'none', 'XTick', [], 'Visible', 'off');
% LGND2 = legend([H1 H2], {'\Re', '\Im'}, 'location', 'bestoutside');
% LGND2.Position(1) = LGND1.Position(1);
% LGND2.Position(2) = LGND1.Position(2) - LGND2.Position(4) - 0.02;
% LGND2.Box = 'off';
% LGND2.Color = 'none';
% hold off;
%     
% xlabel('f [THz]');
% ylabel('Z [\Omega]');
% title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
%        num2str(d * 1e6) ' \mum']);

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

function Znm = eval_Znm(f, slot, d)
    Z11 = slot.compute_zin(f);
    Z12 = eval_Zm(f, slot, d);

    Znm{1,1} = Z11;
    Znm{1,2} = Z12;
    Znm{2,1} = Z12;
    Znm{2,2} = Z11;
end

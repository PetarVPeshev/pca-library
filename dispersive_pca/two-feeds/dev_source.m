close all;
clear;
clc;

env;
sim_params;

% c = physconst('LightSpeed');

LineLegend = {'\Re', '\Im'};
Location   = 'bestoutside';
LineWidth  = 1.5;
Color      = ["#0072BD", "#EDB120", "#77AC30", "#A2142F"];
Position   = [680 558 700 420];

%% PARAMETERS
f = linspace(0.05, 2, 4001) * 1e12;
% SLOT
d_gap = 4.5 * 1e-6;
ws    = (10 : - 2.5 : 2.5) * 1e-6;
% NUMBER OF POINTS
Nf  = length(f);
Nws = length(ws);

%% INPUT IMPEDANCE
Zin         = NaN(Nws, Nf);
ColorLegend = cell(1, 4);
tic;
for idx = 1 : Nws
    slot = create_slot(d_gap, ws(idx), params_slot);

    Zin(idx, :) = eval_Zin(slot, f);

    ColorLegend(idx) = cellstr( ['w_{s} = ' num2str(ws(idx) * 1e6) ' \mum'] );
end
toc

% PLOT MULTIPLE WIDTH
XLabel = 'f [THz]';
YLabel = 'z_{in} [\Omega]';
Title = ['@ \Delta = ' num2str(d_gap * 1e6) ' \mum'];
YLim = [0 80];

plot_two(f * 1e-12, real(Zin), imag(Zin), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'YLim', YLim, 'Position', Position);
grid on;
yticks(0 : 20 : 80);

% PLOT SINGLE WIDTH
% figure('Position', [250 250 750 450]);
% 
% plot(f * 1e-12, real(Zin(1, :)), 'LineWidth', 1.5, 'DisplayName', '\Re');
% hold on;
% plot(f * 1e-12, imag(Zin(1, :)), '--', 'LineWidth', 1.5, 'DisplayName', '\Im');
% 
% grid on;
% legend('location', 'bestoutside');
% ylim([0 80]);
% 
% xlabel('f [THz]');
% ylabel('Z_{in} [\Omega]');
% title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum']);

%% FUNCTIONS
function slot = create_slot(ws, d_gap, params_slot)
    slot = SlotInDielectrics(d_gap, ws, params_slot.er_up, params_slot.er_dn);
end

function a = get_integration_domain(d_gap, f, N)
    arguments
        d_gap (1,1) double {mustBePositive}
        f     (1,1) double {mustBePositive}
        N     (1,1) double {mustBePositive} = 15
    end

    c = physconst('LightSpeed');
    wlen = c / f;

    a = N * wlen / d_gap;
end

function D = eval_D(kx, k0, ws, er_dn, er_up)
    eta0 = get_phys_const('VacuumImpedance');
    
    k1 = k0 * sqrt(er_dn);
    k2 = k0 * sqrt(er_up);

    % Top sheet
    K1 = - 1j * sqrt(- k1 ^ 2 + kx .^ 2);
    K2 = - 1j * sqrt(- k2 ^ 2 + kx .^ 2);

    D = (K1 .^ 2) .* compute_J0H02(K1 * ws / 4) + (K2 .^ 2) .* compute_J0H02(K2 * ws / 4);
    D = D / (2 * k0 * eta0);
end

function Zin = eval_Zin(slot, f)
    c0 = physconst('LightSpeed');

    NFT = length(f);
    d_gap = slot.d_gap;
    ws    = slot.ws;
    er_dn = slot.er_dn;
    er_up = slot.er_up;
    
    Zin = NaN(1, NFT);
    parfor idx = 1 : NFT
        k0 = 2 * pi * f(idx) / c0;

        % INTEGRATION DOMAIN
        kx_max = get_integration_domain(d_gap, f(idx));

        start_pt = - (kx_max * k0 + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        % INPUT IMPEDANCE
        integrand = @(kx) (sinc(kx * d_gap / (2 * pi)) .^ 2) ./ eval_D(kx, k0, ws, er_dn, er_up);
        Zin(idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
    end
end

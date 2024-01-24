close all;
clear;
clc;

addpath('..\slots');
addpath('..\utils');

%% PARAMETERS
% f = linspace(eps, 2, 4001) * 1e12;
f = (0.05 : 0.005 : 2) * 1e12;
% t = (-2 : 0.001 : 10) * 1e-12;
t = linspace(-2, 10, 2001) * 1e-12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1
% SLOT POSITIONS 
dx = d_gap;
x = (0 : 100 : 300) * 1e-6;

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% IMPEDANCE
Zin = NaN(1, length(f));
Z = NaN(length(x), length(f));
for f_idx = 1 : 1 : length(f)
    Zin(f_idx) = slot.compute_zin(f(f_idx));
    Z(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx)';
end

Zplt = NaN(length(x) + 1, length(f));
Zplt(1, :) = Zin;
Zplt(2 : length(x) + 1, :) = Z;

ColorLegend = cell(1, 5);
ColorLegend(1) = {'Z_{in}'};
for x_idx = 2 : 1 : length(x) + 1
    ColorLegend(x_idx) = cellstr(['Z_{x}, x = ' num2str(x(x_idx - 1) * 1e6) ' \mum']);
end
LineLegend = {'\Re', '\Im'};
Location = 'bestoutside';
LineWidth = 1.5;
Color = ["#0072BD", "#EDB120", "#7E2F8E", "#77AC30", "#A2142F"];
Position = [680 558 700 420];
XLabel = 'f [THz]';
YLabel = 'Z [\Omega]';
Title = ['@ \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' num2str(dx * 1e6) ' \mum'];
YLim = [-20 80];

plot_two(f * 1e-12, real(Zplt), imag(Zplt), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'YLim', YLim, 'Position', Position);
grid on;

%% ADMITTANCE
g = NaN(length(x) + 1, length(t));
for x_idx = 1 : 1 : length(x)
    if x_idx == 1
        G = 1 ./ Zplt(x_idx, :);
    else
        G = Zplt(x_idx, :) ./ (Zplt(2, :) .^ 2);
    end
    g(x_idx, :) = 2 * real(eval_IFT(t, f, G));
end

figure('Position', [680 558 700 420]);

plot(t * 1e12, g(1, :), 'LineWidth', 1.5, 'DisplayName', 'g_{in}');
hold on;
for x_idx = 3 : 1 : length(x) + 1
    plot(t * 1e12, g(x_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', ['g_{x} / w_{x}, x = ' num2str(x(x_idx - 1) * 1e6) ' \mum']);
    hold on;
end
hold off;

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('g [S]');
title(['w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum']);

%% WEIGHT
w = NaN(2, length(t));
W = 1 ./ (Zplt(1, :) .^ 2);
w(1, :) = 2 * real(eval_IFT(t, f, W));
Wm = 1 ./ (Zplt(2, :) .^ 2);
w(2, :) = 2 * real(eval_IFT(t, f, Wm));

figure('Position', [680 558 700 420]);

plot(t * 1e12, w(1, :), 'Color', "#0072BD", 'LineWidth', 1.5, 'DisplayName', 'w');
hold on;
plot(t * 1e12, w(2, :), '--', 'Color', "#EDB120", 'LineWidth', 1.5, 'DisplayName', 'w_{m}');

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('w [S^{2}]');
title(['w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum']);

%% FUNCTIONS
function Zx = evaluate_Zx(x, f, slot, dx)
    c0 = get_phys_const('LightSpeed');
    k0 = 2 * pi * f / c0;

    D = @(kx) slot.compute_D(kx, f, 'TopSheet');
    F = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx = @(kx) sinc(- kx * dx / (2 * pi));

    Zx = NaN(1, length(x));
    for x_idx = 1 : 1 : length(x)
        if x(x_idx) == 0
            start_pt = (1000 * k0 + 1j * 0.01);
            end_pt = - start_pt;
            waypts = [-(1 + 1j) (1 + 1j)] * 0.01;

            integrand = @(kx) F(kx) .* Fx(kx) ./ D(kx);
        else
            start_pt = - (0.01 + 1j * 50 * k0);
            end_pt = k0 + 0.01 - 1j * 50 * k0;
            waypts = [-(0.01 + 1j * 0.01) (0.01 + 1j * 0.01) (k0 + 1j * 0.01) (k0 + 0.01)];

            integrand = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D(kx);
        end
        
        Zx(x_idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
    end
end

function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

function Zx = evaluate_Zx_analyt(x, f, slot, dx)
    kxp = slot.find_kxp(f);

    Dp = slot.compute_Dp(kxp, f, 'BottomSheet');
    F = sinc(kxp * slot.d_gap / (2 * pi));
    Fx = sinc(- kxp * dx / (2 * pi));

    Zx = - 1j * F * Fx * exp(- 1j * kxp * abs(x)) / Dp;
end

% OLD FUNCTIONS
function D = compute_D(kx, f, ws, er_up, er_dn, sheet)
    %COMPUTE_D Summary of this method goes here
    %   Detailed explanation goes here
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    k0 = 2 * pi * f / c0;
    k1 = k0 * sqrt(er_dn);
    k2 = k0 * sqrt(er_up);

    % Solution in top Riemann sheet
    K1 = - 1j * sqrt(- k1 ^ 2 + kx .^ 2);
    K2 = - 1j * sqrt(- k2 ^ 2 + kx .^ 2);

    if strcmp(sheet, 'BottomSheet')
        K2 = - K2;
    end

    D = (K1 .^ 2) .* compute_J0H02(K1 * ws / 4) + (K2 .^ 2) .* compute_J0H02(K2 * ws / 4);
    D = D / (2 * k0 * eta0);
end

close all;
clear;
clc;

addpath('..\slots', '..\utils');

Position   = [680 558 700 420];
LineLegend = {'\Re', '\Im'};
Location   = 'bestoutside';
LineWidth  = 1.5;
Color      = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30"];

%% PARAMETERS
f = linspace(eps, 15, 4001) * 1e12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1
% SLOT POSITIONS
dx = 0.5 * 1e-6;
x = [0.5 100] * 1e-6;
% NUMBER OF POINTS
Nf = length(f);
Nx = length(x);

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% FREQUENCY DOMAIN
Zx = NaN(Nx, Nf);
for f_idx = 1 : 1 : Nf
    Zx(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx);
end

ColorLegend = cell(1, Nx);
for x_idx = 1 : 1 : Nx
    ColorLegend(x_idx) = cellstr(['x = ' num2str(x(x_idx) * 1e6) ' \mum']);
end
Title = ['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum'];
XLabel = 'f [THz]';

YLabel = 'Z_{x} [\Omega]';
YLim   = [-200 500];

plot_two(f * 1e-12, real(Zx), imag(Zx), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'Position', Position, 'YLim', YLim);
grid on;

figure('Position', Position);
plot(f * 1e-12, abs(Zx), 'LineWidth', 1.5);
hold on;
xline(7.5, 'k--', 'LineWidth', 1.5, 'DisplayName', 'MODEL FAILS');

grid on;
legend(ColorLegend, 'location', 'bestoutside');

xlabel(XLabel);
ylabel('|Z_{x}| [\Omega]');
title(Title);

%% FUNCTIONS
function Zx = evaluate_Zx(x, f, slot, dx)
    c0 = get_phys_const('LightSpeed');
    k0 = 2 * pi * f / c0;

    Nx = length(x);

    D  = @(kx) slot.compute_D(kx, k0, 'Top', 'Top');
    F  = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx = @(kx) sinc(- kx * dx / (2 * pi));

    Zx = NaN(1, Nx);
    for x_idx = 1 : 1 : Nx
        start_pt = - (100 * k0 + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        integrand = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D(kx);
        Zx(x_idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
    end
end

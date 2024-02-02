close all;
clear;
clc;

addpath('..\slots', '..\utils');

Position   = [680 558 700 420];
LineLegend = {'Z_{x}', 'Res\{Z_{x}\}'};
Location   = 'bestoutside';
LineWidth  = 1.5;
Color      = ["#0072BD", "#EDB120", "#7E2F8E", "#77AC30", "#A2142F"];

%% PARAMETERS
f = linspace(eps, 2, 4001) * 1e12;
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
x = [4.5 100 200 300] * 1e-6;
% NUMBER OF POINTS
Nf = length(f);
Nt = length(t);
Nx = length(x);

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% FREQUENCY DOMAIN
Zx = NaN(Nx, Nf);
Zx_res = Zx;
for f_idx = 1 : 1 : Nf
    Zx(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx);
    Zx_res(:, f_idx) = evaluate_Zx_analyt(x, f(f_idx), slot, dx);
end

ColorLegend = cell(1, Nx);
for x_idx = 1 : 1 : Nx
    ColorLegend(x_idx) = cellstr(['x = ' num2str(x(x_idx) * 1e6) ' \mum']);
end
Title = ['\Re @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum'];
YLabel = '\Re\{Z_{x}\} [\Omega]';
XLabel = 'f [THz]';

plot_two(f * 1e-12, real(Zx), real(Zx_res), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'Position', Position);
grid on;

Title = ['\Im @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum'];
YLabel = '\Im\{Z_{x}\} [\Omega]';

plot_two(f * 1e-12, imag(Zx), imag(Zx_res), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'Position', Position);
grid on;

%% TIME DOMAIN
zx = NaN(Nx - 1, Nt);
zx_res = zx;
for x_idx = 1 : 1 : Nx - 1
    zx(x_idx, :) = 2 * real(eval_IFT(t, f, Zx(x_idx + 1, :)));
    zx_res(x_idx, :) = 2 * real(eval_IFT(t, f, Zx_res(x_idx + 1, :)));
end

Title = ['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum'];
YLabel = 'z_{x} [\Omega]';
XLabel = 't [ps]';

plot_two(t * 1e12, zx, zx_res, 'ColorLegend', ColorLegend(2 : Nx), 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color(2 : Nx), 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'Position', Position);
grid on;

%% FUNCTIONS
function Zx = evaluate_Zx(x, f, slot, dx)
    c0 = get_phys_const('LightSpeed');
    k0 = 2 * pi * f / c0;

    Nx = length(x);

    D  = @(kx) slot.compute_D(kx, f, 'TopSheet');
    F  = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx = @(kx) sinc(- kx * dx / (2 * pi));

    Zx = NaN(1, Nx);
    for x_idx = 1 : 1 : Nx
        start_pt = - (1000 * k0 + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        integrand = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D(kx);
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

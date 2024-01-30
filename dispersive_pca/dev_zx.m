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

    D_tt = @(kx) slot.compute_D(kx, f, 'TopSheet');
    D_tb = @(kx) slot.compute_D(kx, f, 'BottomSheet');
    D_bb = @(kx) compute_Dbb(slot, kx, f);
    F = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx = @(kx) sinc(- kx * dx / (2 * pi));

    Zx = NaN(1, Nx);
    for x_idx = 1 : 1 : Nx
        if x(x_idx) == 0
            start_pt = - (1000 * k0 + 1j * 0.01);
            end_pt = - start_pt;
            waypts = [-(1 + 1j) (1 + 1j)] * 0.01;

            integrand = @(kx) F(kx) .* Fx(kx) ./ D_tt(kx);
            Zx(x_idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        else
            k2 = k0 * sqrt(slot.er_up);
            kxp = slot.find_kxp(f);

            integrand_tt = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D_tt(kx);
            integrand_tb = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D_tb(kx);
            integrand_bb = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D_bb(kx);

            %% integration path 1
%             % Top-Top sheet integral 1
%             start_pt = - (0.01 + 1j * 50 * k0);
%             end_pt = k0 + 0.01;
%             waypts = [-(0.01 + 1j * 0.01) (0.01 + 1j * 0.01) (k0 + 1j * 0.01)];
%             int_p1 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
% 
%             % Top-Bottom sheet
%             start_pt = k0 + 0.01;
%             end_pt = k2 - 0.01;
%             waypts = [(k0 - 1j * 0.01) (0.01 - 1j * 0.01) (0.01 - 1j * 50 * k0) (k2 - 0.01 - 1j * 50 * k0) ...
%                       (k2 - 0.01 + 1j * (imag(kxp) - 1)) (real(kxp) - 1 + 1j * (imag(kxp) - 1)) ...
%                       (real(kxp) - 1 + 1j * (imag(kxp) + 1)) (k2 - 0.01 + 1j * (imag(kxp) + 1))];
%             int_p2 = integral(integrand_tb, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
% 
%             % Top-Top sheet integral 2
%             start_pt = k2 - 0.01;
%             end_pt = k2 + 0.01 - 1j * 50 * k0;
%             waypts = [(k2 - 0.01 + 1j * 0.01) (k2 + 0.01 + 1j * 0.01)];
%             int_p3 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
% 
%             Zx(x_idx) = int_p1 + int_p2 + int_p3;

            %% integration path 2
%             % Top-Top sheet integral 1
%             start_pt = - (0.01 + 1j * 50 * k0);
%             end_pt = k0 + 0.01;
%             waypts = [-(0.01 + 1j * 0.01) (0.01 + 1j * 0.01) (k0 + 1j * 0.01)];
%             int_p1 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
% 
%             % Top-Bottom sheet
%             start_pt = k0 + 0.01;
%             end_pt = k2 - 0.01;
%             waypts = [(k0 + 0.01 - 1j * 50 * k0) (k2 - 0.01 - 1j * 50 * k0) ...
%                       (k2 - 0.01 + 1j * (imag(kxp) - 1)) (real(kxp) - 1 + 1j * (imag(kxp) - 1)) ...
%                       (real(kxp) - 1 + 1j * (imag(kxp) + 1)) (k2 - 0.01 + 1j * (imag(kxp) + 1))];
%             int_p2 = integral(integrand_tb, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
% 
%             % Top-Top sheet integral 2
%             start_pt = k2 - 0.01;
%             end_pt = k2 + 0.01 - 1j * 50 * k0;
%             waypts = [(k2 - 0.01 + 1j * 0.01) (k2 + 0.01 + 1j * 0.01)];
%             int_p3 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
% 
%             Zx(x_idx) = int_p1 + int_p2 + int_p3;

            %% integration path 3
            % Bottom-Bottom sheet integral
            start_pt = k0 - 0.01 - 1j * 50 * k0;
            end_pt = k0 - 0.01;
            int_p1 = integral(integrand_bb, start_pt, end_pt) / (2 * pi);

            % Top-Top sheet integral 1
            start_pt = k0 - 0.01;
            end_pt = k0 + 0.01;
            waypts = [(k0 - 0.01 + 1j * 0.01) (k0 + 0.01 + 1j * 0.01)];
            int_p2 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);

            % Top-Bottom sheet integral
            start_pt = k0 + 0.01;
            end_pt = k2 - 0.01;
            waypts = [(k0 + 0.01 - 1j * 50 * k0) (k2 - 0.01 - 1j * 50 * k0) ...
                      (k2 - 0.01 + 1j * (imag(kxp) - 1)) (real(kxp) - 1 + 1j * (imag(kxp) - 1)) ...
                      (real(kxp) - 1 + 1j * (imag(kxp) + 1)) (k2 - 0.01 + 1j * (imag(kxp) + 1))];
            int_p3 = integral(integrand_tb, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);

            % Top-Top sheet integral 2
            start_pt = k2 - 0.01;
            end_pt = k2 + 0.01 - 1j * 50 * k0;
            waypts = [(k2 - 0.01 + 1j * 0.01) (k2 + 0.01 + 1j * 0.01)];
            int_p4 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);

            Zx(x_idx) = int_p1 + int_p2 + int_p3 + int_p4;
        end
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

function D = compute_Dbb(slot, kx, f)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    k0 = 2 * pi * f / c0;
    k1 = k0 * sqrt(slot.er_dn);
    k2 = k0 * sqrt(slot.er_up);

    % Solution in bottom-bottom Riemann sheet
    K1 = 1j * sqrt(- k1 ^ 2 + kx .^ 2);
    K2 = 1j * sqrt(- k2 ^ 2 + kx .^ 2);

    D = (K1 .^ 2) .* compute_J0H02(K1 * slot.ws / 4) + (K2 .^ 2) .* compute_J0H02(K2 * slot.ws / 4);
    D = D / (2 * k0 * eta0);
end

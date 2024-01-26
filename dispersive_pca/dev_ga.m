close all;
clear;
clc;

addpath('..\slots');
addpath('..\utils');

%% PARAMETERS
t = linspace(-2, 10, 2001) * 1e-12;
f = (0.05 : 0.005 : 2) * 1e12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1
% SLOT POSITIONS 
dx = d_gap;
x = (100 : 100 : 300) * 1e-6;
% NUMBER OF POINTS
Nf = length(f);
Nt = length(t);
Nx = length(x);

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% WEIGHT
% frequency domain
W = NaN(1, Nf);
Wm = W;
for f_idx = 1 : 1 : Nf
    W(f_idx) = slot.compute_zin(f(f_idx));
    Wm(f_idx) = evaluate_Zx(0, f(f_idx), slot, dx);
end
W = (1 ./ W) .^ 2;
Wm = (1 ./ Wm) .^ 2;

% time domain
w = 2 * real(eval_IFT(t, f, W));
wx = 2 * real(eval_IFT(t, f, Wm));

figure('Position', [680 558 700 420]);

plot(t * 1e12, w, 'Color', "#0072BD", 'LineWidth', 1.5, 'DisplayName', 'w');
hold on;
plot(t * 1e12, wx, '--', 'Color', "#EDB120", 'LineWidth', 1.5, 'DisplayName', 'w_{x}');

grid on;
ylim([min(w) max(w)]);
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('w [S^{2}]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum']);

%% IMPEDANCE
Zin = NaN(1, Nf);
Zm = NaN(Nx, Nf);
for f_idx = 1 : 1 : Nf
    Zin(f_idx) = slot.compute_zin(f(f_idx));
    Zm(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx);
end

%% ADMITTANCE
Hin = 1 ./ Zin;
Hm = Zm .* Wm;

% time domain
hin = 2 * real(eval_IFT(t, f, Hin));
hm = NaN(Nx, Nt);
for x_idx = 1 : 1 : Nx
    hm(x_idx, :) = 2 * real(eval_IFT(t, f, Hm(x_idx, :)));
end

figure('Position', [680 558 700 420]);

plot(t * 1e12, hin, 'LineWidth', 1.5, 'DisplayName', 'g_{in}');
hold on;
for x_idx = 1 : 1 : Nx
    plot(t * 1e12, hm(x_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', ['h_{m}, x = ' num2str(x(x_idx) * 1e6) ' \mum']);
    hold on;
end
hold off;

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('g [S]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum']);

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
close all;
clear;
clc;

addpath('../utils');

%% PARAMETERS
% DIELECTRICS
er_up = [11.7 4 2.54];  % er2
er_dn = 1;              % er1
% FREQUENCY
f = 10 * 1e9;
wlen = physconst('LightSpeed') / f;
k0 = 2 * pi / wlen;
% SLOT WIDTH
ws = (0.001 : 0.001 : 0.1) * wlen;

%% LEAKY-WAVE PROPAGATION CONSTANT
kx_lw = NaN(length(er_up), length(ws));
kx_lw_analyt = NaN(length(er_up), length(ws));
for er_idx = 1 : 1 : length(er_up)
    for ws_idx = 1 : 1 : length(ws)
        kx_lw(er_idx, ws_idx) = find_kxp(f, ws(ws_idx), er_up(er_idx), er_dn);
        kx_lw_analyt(er_idx, ws_idx) = evaluate_kxp_analyt(f, ws(ws_idx), er_up(er_idx), er_dn);
    end
end

figure('Name', 'Leaky-Wave Propagation Constant', 'Position', [150 150 800 500]);

subplot(2, 1, 1);
plot(ws / wlen, real(kx_lw(1, :) / k0), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws / wlen, real(kx_lw_analyt(1, :) / k0), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws / wlen, real(kx_lw(2, :) / k0), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws / wlen, real(kx_lw_analyt(2, :) / k0), '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws / wlen, real(kx_lw(3, :) / k0), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
hold on;
plot(ws / wlen, real(kx_lw_analyt(3, :) / k0), '--', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
grid on;
legend('location', 'bestoutside');
ylabel('\Re\{k_{x}^{LW} / k_{0}\}');
ylim([1 3.5]);
xticks(0 : 0.02 : 0.1);

subplot(2, 1, 2);
plot(ws / wlen, imag(kx_lw(1, :) / k0), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws / wlen, imag(kx_lw_analyt(1, :) / k0), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws / wlen, imag(kx_lw(2, :) / k0), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws / wlen, imag(kx_lw_analyt(2, :) / k0), '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws / wlen, imag(kx_lw(3, :) / k0), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
hold on;
plot(ws / wlen, imag(kx_lw_analyt(3, :) / k0), '--', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
grid on;
legend('location', 'bestoutside');
ylabel('\Im\{k_{x}^{LW} / k_{0}\}');
xlabel('w_{s} / \lambda_{0}');
ylim([-0.5 0]);
xticks(0 : 0.02 : 0.1);

sgtitle('Leaky-Wave Propagation Constant', 'FontSize', 11, 'FontWeight', 'bold');

%% LEAKY-WAVE MAGNETIC CURRENT
kxp = find_kxp(f, 0.03 * wlen, 11.7, 1);
x = (-1 : 0.001 : 1) * wlen;
vlw = compute_v_lw(x, kxp, f, 0.03 * wlen, 11.7, 1, 0.06 * wlen);

figure('Name', 'Normalized Leaky-Wave Magnetic Current', 'Position', [150 150 800 350]);
plot(x / wlen, real(vlw), 'LineWidth', 1.5, 'DisplayName', '\Re\{v_{LW}\}');
hold on;
plot(x / wlen, imag(vlw), 'LineWidth', 1.5, 'DisplayName', '\Im\{v_{LW}\}');
grid on;
legend('location', 'bestoutside');
xlabel('x / \lambda_{0}');
ylabel('v_{LW} / V');
title(['Normalized Leaky-Wave Magnetic Current @ \Deltaw = 0.06\lambda_{0}, w_{s} = 0.03\lambda_{0}, ' ...
    '\epsilon_{r,2} = 11.7, \epsilon_{r,1} = 1']);

%% FRINGE MAGNETIC CURRENT
v = NaN(1, length(x));
for x_idx = 1 : 1 : length(x)
    v(x_idx) = compute_v(x(x_idx), f, 0.03 * wlen, 11.7, 1, 0.06 * wlen);
end

figure('Name', 'Normalized Fringe Magnetic Current', 'Position', [150 150 800 350]);
plot(x / wlen, real(v - vlw), 'LineWidth', 1.5, 'DisplayName', '\Re\{v_{fc}\}');
hold on;
plot(x / wlen, imag(v - vlw), 'LineWidth', 1.5, 'DisplayName', '\Im\{v_{fc}\}');
grid on;
legend('location', 'bestoutside');
xlabel('x / \lambda_{0}');
ylabel('v_{fc} / V');
title(['Normalized Fringe Magnetic Current @ \Deltaw = 0.06\lambda_{0}, w_{s} = 0.03\lambda_{0}, ' ...
    '\epsilon_{r,2} = 11.7, \epsilon_{r,1} = 1']);

%% FUNCTIONS
function v = compute_v(x, f0, ws, er1, er2, dw)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    const = 1 / (2 * k0 * eta0);
    delta = 0.01;

    D1 = @(kx) (k1 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4);
    D2 = @(kx) (k2 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4);
    D = @(kx) const * ( D1(kx) + D2(kx) );
    integrand = @(kx) exp(- 1j * kx * x) .* sinc(kx * dw / (2 * pi)) ./ D(kx);
    
    start_int = - 50 * k0 - 1j * delta;
    end_int = 50 * k0 + 1j * delta;

    % Try Dunja's waypoints
    waypts = [-(1 + 1j) (1 + 1j)] * delta;

    v = integral(integrand, start_int, end_int, 'Waypoints', waypts);
    v = v / (2 * pi);
end

function v = compute_v_lw(x, kxp, f0, ws, er1, er2, dw)
    Dp = compute_Dp(kxp, f0, ws, er1, er2);
    v = - 1j * sinc(kxp * dw / (2 * pi)) * exp(- 1j * kxp * abs(x)) / Dp;
end

function kxp = evaluate_kxp_analyt(f0, ws, er1, er2)
    c0 = get_phys_const('LightSpeed');
    
    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    beta = sqrt( (k1 ^ 2 + k2 ^ 2) / 2 );
    kd = sqrt( (k1 ^ 2 - k2 ^ 2) / 2 );

    kxp = beta + ( (kd ^ 2) / beta ) * (compute_J0H2(- ws * 1j * kd / 4) - compute_J0H2(ws * kd / 4)) ...
        / (compute_A(- ws * 1j * kd / 4) + compute_A(ws * kd / 4));
end

function A = compute_A(z)
    A = - 2 * compute_J0H2(z) + z * besselj(1, z) * besselh(0, 2, z) + z * besselh(1, 2, z) * besselj(0, z);
end

function kxp = find_kxp(f0, ws, er1, er2)
    c0   = get_phys_const('LightSpeed');

    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);

    kxp = sqrt( (k1 ^ 2 + k2 ^ 2) / 2 );
    err = 1;

    while err > 1e-5
        kxp_prev = kxp;

        D = compute_D(kxp_prev, f0, ws, er1, er2);
        Dp = compute_Dp(kxp_prev, f0, ws, er1, er2);
        kxp = kxp_prev - D / Dp;

        err = abs(kxp - kxp_prev);
    end
end

function Dp = compute_Dp(kx, f0, ws, er1, er2)
    c0   = get_phys_const('LightSpeed');

    k0 = 2 * pi * f0 / c0;

    delta = 500;
    delta_k = k0 / delta;

    D1 = compute_D(kx + delta_k / 2, f0, ws, er1, er2);
    D2 = compute_D(kx - delta_k / 2, f0, ws, er1, er2);

    Dp = (D1 - D2) / delta_k;
end

function D = compute_D(kx, f0, ws, er1, er2)
    c0   = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');

    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);

    const = 1 / (2 * k0 * eta0);

    % Valid solution is when imag(K1) < 0 and imag(K2) > 0
    K1 = 1j * sqrt(- k1 ^ 2 + kx ^ 2);
    K2 = - 1j * sqrt(- k2 ^ 2 + kx ^ 2);

    D1 = (K1 ^ 2) * compute_J0H2(K1 * ws / 4);
    D2 = (K2 ^ 2) * compute_J0H2(K2 * ws / 4);

    D = const * (D1 + D2);
end

function J0H2 = compute_J0H2(z)
    J0H2 = besselj(0, z) .* besselh(0, 2, z);
end

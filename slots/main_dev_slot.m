close all;
clear;
clc;

addpath('../utils');

%% PARAMETERS
% DIELECTRICS
er_up = 11.7;
er_dn = 1;
% FREQUENCY
f = 10 * 1e9;
wlen = physconst('LightSpeed') / f;
k0 = 2 * pi / wlen;
% SLOT WIDTH
ws = 0.03 * wlen;
% d_gap = 0.06 * wlen;
% SOURCE WIDTH
dw = 0.06 * wlen;

%% INFINITE SLOT BETWEEN DIELECTRICS
% slot = SlotBetweenDielectrics(er_up, er_dn, f, ws);

%% TEST DISPERSION EQUATION FUNCTION
% D = slot.compute_D(5);
% PROPAGATION CONSTANTS IN MEDIUMS
k1 = 2 * pi * f * sqrt(er_up) / get_phys_const('LightSpeed');
k2 = 2 * pi * f * sqrt(er_dn) / get_phys_const('LightSpeed');
% DISPERSION EQUATION
kx_g = sqrt( (k1 ^ 2 + k2 ^ 2) / 2 );
D = compute_D(kx_g, f, ws, er_up, er_dn);
Dp = compute_D_deriv(kx_g, f, ws, er_up, er_dn);
% VOLTAGE AS FUNCTION OF X
x = (-1 : 0.001 : 1) * wlen;
v = NaN(1, length(x));
for x_idx = 1 : 1 : length(x)
    v(x_idx) = compute_v(x(x_idx), f, ws, er_up, er_dn, dw);
end

figure('Name', 'Normalized Magnetic Current', 'Position', [150 150 800 350]);
plot(x / wlen, real(v), 'LineWidth', 1.5, 'DisplayName', '\Re\{v\}');
hold on;
plot(x / wlen, imag(v), 'LineWidth', 1.5, 'DisplayName', '\Im\{v\}');
grid on;
legend('location', 'bestoutside');
xlabel('x / \lambda_{0}');
ylabel('v / V');
title(['Normalized Magnetic Current @ \Deltaw = 0.06\lambda_{0}, w_{s} = 0.03\lambda_{0}, ' ...
    '\epsilon_{r,2} = 11.7, \epsilon_{r,1} = 1']);

% LEAKY-WAVE PROPAGATION CONSTANT
ws_vec = (0.001 : 0.001 : 0.1) * wlen;
er_vec = [11.7 4 2.54];
kx_lw = NaN(length(er_vec), length(ws_vec));
kx_lw_analyt = NaN(length(er_vec), length(ws_vec));
for er_idx = 1 : 1 : length(er_vec)
    for ws_idx = 1 : 1 : length(ws_vec)
        kx_lw(er_idx, ws_idx) = find_kxp(f, ws_vec(ws_idx), er_vec(er_idx), er_dn);
        kx_lw_analyt(er_idx, ws_idx) = evaluate_kxp_analyt(f, ws_vec(ws_idx), er_vec(er_idx), er_dn);
    end
end

figure('Name', 'Leaky-Wave Propagation Constant', 'Position', [150 150 800 500]);

subplot(2, 1, 1);
plot(ws_vec / wlen, real(kx_lw(1, :) / k0), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws_vec / wlen, real(kx_lw_analyt(1, :) / k0), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws_vec / wlen, real(kx_lw(2, :) / k0), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws_vec / wlen, real(kx_lw_analyt(2, :) / k0), '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws_vec / wlen, real(kx_lw(3, :) / k0), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
hold on;
plot(ws_vec / wlen, real(kx_lw_analyt(3, :) / k0), '--', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
grid on;
legend('location', 'bestoutside');
ylabel('\Re\{k_{x}^{LW} / k_{0}\}');
ylim([1 3.5]);
xticks(0 : 0.02 : 0.1);

subplot(2, 1, 2);
plot(ws_vec / wlen, imag(kx_lw(1, :) / k0), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws_vec / wlen, imag(kx_lw_analyt(1, :) / k0), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 11.7');
hold on;
plot(ws_vec / wlen, imag(kx_lw(2, :) / k0), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws_vec / wlen, imag(kx_lw_analyt(2, :) / k0), '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 4');
hold on;
plot(ws_vec / wlen, imag(kx_lw(3, :) / k0), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
hold on;
plot(ws_vec / wlen, imag(kx_lw_analyt(3, :) / k0), '--', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'DisplayName', '\epsilon_{r,2} = 2.54');
grid on;
legend('location', 'bestoutside');
ylabel('\Im\{k_{x}^{LW} / k_{0}\}');
xlabel('w_{s} / \lambda_{0}');
ylim([-0.5 0]);
xticks(0 : 0.02 : 0.1);

sgtitle('Leaky-Wave Propagation Constant', 'FontSize', 11, 'FontWeight', 'bold');

% DISPERSION EQUATION PEAKS
% kx_vec = (1 : 0.001 : 5) * k0;
% D_vec = NaN(1, length(kx_vec));
% for kx_idx = 1 : 1 : length(kx_vec)
%     D_vec(kx_idx) = compute_D(kx_vec(kx_idx), f, ws, er_up, er_dn);
% end
% 
% figure()
% findpeaks( imag(1 ./ D_vec) )

%% RIEMANN'S SHEET
% [A, theta] = meshgrid((1 : 1 : 200) * k0, (-6 : 0.1 : 6) * pi);
% kx = A .* exp(1j * theta);
% D_mesh = NaN(size(A, 1, 2));
% for row = 1 : 1 : size(A, 1)
%     for col = 1 : 1 : size(A, 2)
%         D_mesh(row, col) = compute_D(kx(row, col), f, ws, er_up, er_dn);
%     end
% end
% 
% figure();
% surf(real(kx), imag(kx), imag(D_mesh), 'LineStyle', 'none');
% colormap jet;
% colorbar;
% zlim([0.0050 0.0125]);

% r = (0:1:30)';                           % create a matrix of complex inputs
% theta = pi*(-2:0.05:2);
% z = r*exp(1i*theta);
% w = z.^(1/2);                          % calculate the complex outputs
% w = sqrt(r)*exp(1i*theta/2);
% w = log(r) + 1j * theta;
% w = sqrt(r) * exp(1i*theta/2) .* ( log(r) + 1j * theta );
% zp = r + 0.01j * exp(- 0.5 * (r) .^ 2);

% figure('Name','Graphique complexe','units','normalized','outerposition',[ 0.08 0.1 0.8 0.55]);
% subplot(1, 2, 1)
% 
% surf(real(z),imag(z),real(w),'LineStyle','-')    % visualize the complex function using surf
% % hold on;
% % surf(real(z),imag(z),real(zp),'LineStyle','-')    % visualize the complex function using surf
% xlabel('Real(z)')
% ylabel('Imag(z)')
% zlabel('Real(u)')
% cb = colorbar;
% colormap jet;                            % gradient from blue to red
% cb.Label.String = 'Imag(v)';
% 
% subplot(1, 2, 2)
% surf(real(z),imag(z),imag(w),'LineStyle','-')    % visualize the complex function using surf
% xlabel('Real(z)')
% ylabel('Imag(z)')
% zlabel('Imag(v)')
% cb = colorbar;
% colormap jet;                            % gradient from blue to red
% cb.Label.String = 'Real(u)';

%% FUNCTIONS
function kxp = find_kxp(f0, ws, er1, er2)
    c0 = get_phys_const('LightSpeed');
    
    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    kxp_prev = 0;
%     kxp_prev = evaluate_kxp_analyt(f0, ws, er1, er2);
    kxp = sqrt((k1 ^ 2 + k2 ^ 2) / 2);

    while abs(kxp - kxp_prev) > 1e-5
        kxp_prev = kxp;
        
        D = compute_D(kxp_prev, f0, ws, er1, er2);
        Dp = compute_D_deriv(kxp_prev, f0, ws, er1, er2);
%         Dp = evaluate_D_deriv_analyt(kxp_prev, f0, ws, er1, er2);

        kxp = kxp_prev - D / Dp;
    end
end

function kxp = evaluate_kxp_analyt(f0, ws, er1, er2)
    c0 = get_phys_const('LightSpeed');
    
    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    beta = sqrt( (k1 ^ 2 + k2 ^ 2) / 2 );
    kd = sqrt( (k1 ^ 2 - k2 ^ 2) / 2 );

    kxp = beta + ( (kd ^ 2) / beta ) * (evaluate_fun(- ws * 1j * kd / 4, 'B') - evaluate_fun(ws * kd / 4, 'B')) / (evaluate_fun(- ws * 1j * kd / 4, 'A') + evaluate_fun(ws * kd / 4, 'A'));
end

function Dp = compute_D_deriv(kx, f0, ws, er1, er2, options)
    arguments
        kx (1,1) double
        f0 (1,1) double
        ws (1,1) double
        er1 (1,1) double
        er2 (1,1) double
        options.delta (1,1) double = 500
    end

    c0 = get_phys_const('LightSpeed');
    
    k0 = 2 * pi * f0 / c0;
    delta_k = k0 / options.delta;
    
    D1 = compute_D(kx + delta_k / 2, f0, ws, er1, er2);
    D2 = compute_D(kx - delta_k / 2, f0, ws, er1, er2);
    Dp = (D1 - D2) / delta_k;
end

function Dp = evaluate_D_deriv_analyt(kx, f0, ws, er1, er2)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);

    K1 = - 1j * sqrt(- k1 ^ 2 + kx ^ 2);
    K2 = - 1j * sqrt(- k2 ^ 2 + kx ^ 2);
    
    const = 1 / (2 * k0 * eta0);

    Dp = const * kx * ( evaluate_fun(ws * K1 / 4, 'A') + evaluate_fun(ws * K2 / 4, 'A') );
end

function D = compute_D(kx, f0, ws, er1, er2)
    c0   = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');

    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);

    const = 1 / (2 * k0 * eta0);

    % Need to force evaluation in top Riemann sheet
    D1 = (k1 ^ 2 - kx ^ 2) * besselj(0, 1j * ws * sqrt(- k1 ^ 2 + kx ^ 2) / 4) * besselh(0, 2, 1j * ws * sqrt(- k1 ^ 2 + kx ^ 2) / 4);
    D2 = (k2 ^ 2 - kx ^ 2) * besselj(0, - 1j * ws * sqrt(- k2 ^ 2 + kx ^ 2) / 4) * besselh(0, 2, - 1j * ws * sqrt(- k2 ^ 2 + kx ^ 2) / 4);

    D = const * (D1 + D2);
end

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

function J0H2 = compute_J0H2(z)
    J0H2 = besselj(0, z) .* besselh(0, 2, z);
end

function fun = evaluate_fun(z, fun_type)
    if strcmp(fun_type, 'A')
        B = evaluate_fun(z, 'B');
        fun = - 2 * B + z * besselj(1, z) * besselh(0, 2, z) + z * besselh(1, 2, z) * besselj(0, z);
    else 
        fun = besselj(0, z) .* besselh(0, 2, z);
    end
end

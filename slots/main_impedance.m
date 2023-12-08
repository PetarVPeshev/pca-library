close all;
clear;
clc;

addpath('../utils');

%% PARAMETERS
% DIELECTRICS
er_up = 11.7;   % er1
er_dn = 1;      % er2
% FREQUENCY
f = 10 * 1e9;
wlen = physconst('LightSpeed') / f;
k0 = 2 * pi / wlen;
% SLOT WIDTH
ws = [1/30 1/100 1/1000] * wlen;
% FEED WIDTH
dw = (0 : 0.001 : 3) * wlen / 100;

%% SLOT IMPEDANCE
Zin = NaN(length(ws), length(dw));
Zdyn = NaN(length(ws), length(dw));
for ws_idx = 1 : 1 : length(ws)
    for dw_idx = 1 : 1 : length(dw)
        Zin(ws_idx, dw_idx) = compute_Zin(f, ws(ws_idx), er_up, er_dn, dw(dw_idx));
        Zdyn(ws_idx, dw_idx) = compute_Zdyn(f, ws(ws_idx), er_up, er_dn, dw(dw_idx));
    end
end

figure('Name', 'Infinite Slot Impedance', 'Position', [150 150 800 500]);

plot(100 * dw / wlen, real(Zin(1, :)), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{Z_{in}\}, w_{s}=\lambda_{0}/30');
hold on;
plot(100 * dw / wlen, imag(Zin(1, :)), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{Z_{in}\}, w_{s}=\lambda_{0}/30');
hold on;
plot(100 * dw / wlen, real(Zdyn(1, :)), ':', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{Z_{dyn}\}, w_{s}=\lambda_{0}/30');
hold on;
plot(100 * dw / wlen, imag(Zdyn(1, :)), '-.', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{Z_{dyn}\}, w_{s}=\lambda_{0}/30');
hold on;

plot(100 * dw / wlen, real(Zin(2, :)), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{Z_{in}\}, w_{s}=\lambda_{0}/100');
hold on;
plot(100 * dw / wlen, imag(Zin(2, :)), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{Z_{in}\}, w_{s}=\lambda_{0}/100');
hold on;
plot(100 * dw / wlen, real(Zdyn(2, :)), ':', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{Z_{dyn}\}, w_{s}=\lambda_{0}/100');
hold on;
plot(100 * dw / wlen, imag(Zdyn(2, :)), '-.', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{Z_{dyn}\}, w_{s}=\lambda_{0}/100');
hold on;

plot(100 * dw / wlen, real(Zin(3, :)), '--', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Re\{Z_{in}\}, w_{s}=\lambda_{0}/1000');
hold on;
plot(100 * dw / wlen, imag(Zin(3, :)), 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Im\{Z_{in}\}, w_{s}=\lambda_{0}/1000');

grid on;
ylim([-20 100]);
legend('location', 'bestoutside');
xlabel('100 \Delta_{gap} / \lambda_{0}');
ylabel('Z / \Omega');
title('Infinite Slot Impedance @ \epsilon_{r,1} = 1, \epsilon_{r,2} = 11.7');

%% DYNAMIC AND QUASI-STATIC IMPEDANCE
% FUNCTION OF DELTA GAP
dw_vec = linspace(1e-6, 0.1, 101) * wlen;
Zdyn = NaN(1, length(dw_vec));
Zin = NaN(1, length(dw_vec));
for dw_idx = 1 : 1 : length(dw_vec)
    Zdyn(dw_idx) = compute_Zdyn(f, 0.01 * wlen, 11.7, 1, dw_vec(dw_idx));
    Zin(dw_idx) = compute_Zin(f, 0.01 * wlen, 11.7, 1, dw_vec(dw_idx));
end
Zqs = compute_Zqs(f, 0.01 * wlen, 11.7, 1, dw_vec);
Zin_sum = Zdyn + Zqs;

figure('Name', 'Infinite Slot Impedance Components', 'Position', [150 150 800 500]);

plot(dw_vec / wlen, real(Zin), '-+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(dw_vec), 'MarkerSize', 5, 'DisplayName', '\Re\{Z_{in}\}, num');
hold on;
plot(dw_vec / wlen, imag(Zin), '--+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(dw_vec), 'MarkerSize', 5, 'DisplayName', '\Im\{Z_{in}\}, num');
hold on;

plot(dw_vec / wlen, real(Zin_sum), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{Z_{in}\}');
hold on;
plot(dw_vec / wlen, imag(Zin_sum), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{Z_{in}\}');
hold on;
plot(dw_vec / wlen, imag(Zqs), '-.', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Im\{Z_{qs}\}');
hold on;
plot(dw_vec / wlen, imag(Zdyn), ':', 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560], ...
    'DisplayName', '\Im\{Z_{dyn}\}');
% plot(dw_vec / wlen, imag(Zin - Zqs), ':', 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560], ...
%     'DisplayName', '\Im\{Z_{dyn}\}');

grid on;
ylim([-10 50]);
legend('location', 'bestoutside');
xlabel('\Delta_{gap} / \lambda_{0}');
ylabel('Z / \Omega');
title(['Infinite Slot Impedance Components @ \epsilon_{r,1} = 1, \epsilon_{r,2} = 11.7, ' ...
    'w_{s} / \lambda_{0} = 0.01']);

% FUNCTION OF SLOT WIDTH
ws_vec = linspace(1e-6, 0.01, 101) * wlen;
Zdyn = NaN(1, length(ws_vec));
Zin = NaN(1, length(ws_vec));
for ws_idx = 1 : 1 : length(ws_vec)
    Zdyn(ws_idx) = compute_Zdyn(f, ws_vec(ws_idx), 11.7, 1, 0.001 * wlen);
    Zin(ws_idx) = compute_Zin(f, ws_vec(ws_idx), 11.7, 1, 0.001 * wlen);
end
Zqs = compute_Zqs(f, ws_vec, 11.7, 1, 0.001 * wlen);
Zin_sum = Zdyn + Zqs;

figure('Name', 'Infinite Slot Impedance Components', 'Position', [150 150 800 500]);

plot(ws_vec / wlen, real(Zin), '-+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(ws_vec), 'MarkerSize', 5, 'DisplayName', '\Re\{Z_{in}\}, num');
hold on;
plot(ws_vec / wlen, imag(Zin), '--+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(ws_vec), 'MarkerSize', 5, 'DisplayName', '\Im\{Z_{in}\}, num');

hold on;
plot(ws_vec / wlen, real(Zin_sum), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{Z_{in}\}');
hold on;
plot(ws_vec / wlen, imag(Zin_sum), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{Z_{in}\}');
hold on;
plot(ws_vec / wlen, real(Zqs), '-.', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Im\{Z_{qs}\}');
hold on;
plot(ws_vec / wlen, imag(Zdyn), ':', 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560], ...
    'DisplayName', '\Im\{Z_{dyn}\}');
% plot(ws_vec / wlen, imag(Zin - Zqs), ':', 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560], ...
%     'DisplayName', '\Im\{Z_{dyn}\}');

grid on;
ylim([0 40]);
legend('location', 'bestoutside');
xlabel('w_{s} / \lambda_{0}');
ylabel('Z / \Omega');
title(['Infinite Slot Impedance Components @ \epsilon_{r,1} = 1, \epsilon_{r,2} = 11.7, ' ...
    '\Delta / \lambda_{0} = 0.001']);

%% FUNCTIONS
function Zin = compute_Zin(f0, ws, er1, er2, dw)
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
    integrand = @(kx) (sinc(kx * dw / (2 * pi)) .^ 2) ./ D(kx);
    
    % Needs larger integration domain, due to not convergence for small delta gap
    start_int = - 1000 * k0 - 1j * delta;
    end_int = 1000 * k0 + 1j * delta;

    % Try Dunja's waypoints
    waypts = [-(1 + 1j) (1 + 1j)] * delta;

    Zin = integral(integrand, start_int, end_int, 'Waypoints', waypts);
    Zin = Zin / (2 * pi);
end

function Zdyn = compute_Zdyn(f0, ws, er1, er2, dw)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    const = 1 / (2 * k0 * eta0);
    const_inf = const * 4 / (pi * ws);
    delta = 0.01;

    D1 = @(kx) (k1 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4);
    D2 = @(kx) (k2 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4);
    D = @(kx) const * ( D1(kx) + D2(kx) );
    Dinf = @(kx) const_inf * (k1 ^ 2 - k2 ^ 2) ./ ( - 1j * sqrt(- k1 ^ 2 + kx .^ 2) + 1j * sqrt(- k2 ^ 2 + kx .^ 2) );
    integrand = @(kx) (sinc(kx * dw / (2 * pi)) .^ 2) .* (1 ./ D(kx) - 1 ./ Dinf(kx));
    
    % Needs larger integration domain, due to not convergence for small delta gap
    start_int = - 150 * k0 - 1j * delta;
    end_int = 150 * k0 + 1j * delta;

    % Try Dunja's waypoints
    waypts = [-(1 + 1j) (1 + 1j)] * delta;

    Zdyn = integral(integrand, start_int, end_int, 'Waypoints', waypts);
    Zdyn = Zdyn / (2 * pi);
end

function Zqs = compute_Zqs(f0, ws, er1, er2, dw)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    k0 = 2 * pi * f0 / c0;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    const = eta0 * k0 * pi * ws / 2;
    eug = double(eulergamma);
    eug = exp(eug);

    % PAPER IMPLEMENTATION
    Zqs_p = - 1j * ( (k1 ^ 2) * log(k1 / 2) - (k2 ^ 2) * log(k2 / 2) ) / (pi * (k1 ^ 2 - k2 ^ 2));
    Zqs_p = Zqs_p + (1 - 1j * (2 / pi) * (log(eug * dw) - 3 / 2)) / 2;
%     Zqs_p = const * Zqs_p;

    Zqs_pp = besselh(0, 2, k1 * dw) - besselh(0, 2, k2 * dw) + 1j * (2 / pi) * log(k1 / k2);
    Zqs_pp = Zqs_pp ./ (2 * (dw .^ 2) * (k1 ^ 2 - k2 ^ 2));
%     Zqs_pp = const * Zqs_pp;

%     Zqs = Zqs_p + Zqs_pp;
    Zqs = const * (Zqs_p + Zqs_pp);

    % PAPER IMPLEMENTATION RE-WRITTEN
%     Zqs_p = - 2j * ( (k1 ^ 2) * log(k1 / 2) - (k2 ^ 2) * log(k2 / 2) );
%     Zqs_p = Zqs_p - pi * (k1 ^ 2 - k2 ^ 2) * ( 1 - 1j * (2 / pi) * log(eug * dw - 3 / 2) );
%     Zqs_p = Zqs_p / (2 * pi);

%     Zqs_pp = besselh(0, 2, k1 * dw) - besselh(0, 2, k2 * dw) + 1j * (2 / pi) * log(k1 / k2);
%     Zqs_pp = Zqs_pp ./ (2 * (dw .^ 2));

%     Zqs = const * (Zqs_p - Zqs_pp) / (k1 ^ 2 - k2 ^ 2);

    % CROSS IMPLEMENTATION
%     Zqs_p = - 1j * ( (k1 ^ 2) * log(k1 / 2) - (k2 ^ 2) * log(k2 / 2) ) / (pi * (k1 ^ 2 - k2 ^ 2));
%     Zqs_p = Zqs_p + (1 - 1j * (2 / pi) * log(eug * dw) - (2 / pi) * log(eug / 2)) / 2;
% 
%     Zqs_pp = (besselh(0, 2, k1 * dw) - besselh(0, 2, k2 * dw)) ./ (2 * (dw .^ 2) * (k1 ^ 2 - k2 ^ 2));
% 
%     Zqs = const * (Zqs_p + Zqs_pp);

    % DERIVATION
% %     Zqs_p = - 1j * ( (k1 ^ 2) * log(k1) - (k2 ^ 2) * log(k2) ) / (pi * (k1 ^ 2 - k2 ^ 2));
%     Zqs_p = - 1j * ( (k1 ^ 2) * log(k1 / 2) - (k2 ^ 2) * log(k2 / 2) ) / (pi * (k1 ^ 2 - k2 ^ 2));
% %     Zqs_p = Zqs_p + (1 - (2 / pi) * log(eug / 2) - 1j * (2 / pi) * dw .* (log(dw) - 3 / 2)) / 2;
%     Zqs_p = Zqs_p + (1 - 1j * (2 / pi) * (log(eug * dw) - 3 / 2)) / 2;
% 
% %     Zqs_pp = (k1 * besselh(0, 2, k1 * dw) - k2 * besselh(0, 2, k2 * dw)) ./ (2 * (dw .^ 2) * (k1 ^ 2 - k2 ^ 2));
%     Zqs_pp = (besselh(0, 2, k1 * dw) - besselh(0, 2, k2 * dw) + 1j * (2 / pi) * log(k1 / k2)) ./ (2 * (dw .^ 2) * (k1 ^ 2 - k2 ^ 2));
% 
%     Zqs = const * (Zqs_p + Zqs_pp);
end

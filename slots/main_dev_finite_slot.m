close all;
clear;
clc;

addpath('../utils');

%% PARAMETERS
% FREQUENCY
f0 = 10 * 1e9;
wlen = physconst('LightSpeed') / f0;
% SLOT WIDTH
ws = wlen / 50;
% SOURCE WIDTH
dw = wlen / 40;
% TERMINATION DISTANCE
d = wlen / 4;

%% INPUT IMPEDANCE OF SEMI-INFINITE SLOT IN FREE SPACE
f = linspace(0.5, 2, 101) * f0;
Zin = NaN(1, length(f));
for f_idx = 1 : 1 : length(f)
    Zin(f_idx) = compute_Zin_single_finite(f(f_idx), ws, 1, 1, dw, d);
end

figure('Name', 'Semi-Finite Slot Impedance', 'Position', [150 150 700 400]);

plot(f / f0, real(Zin), 'LineWidth', 1.5, 'DisplayName', '\Re\{Z_{in}\}');
hold on;
plot(f / f0, imag(Zin), 'LineWidth', 1.5, 'DisplayName', '\Im\{Z_{in}\}');

grid on;
ylim([-50 150]);
legend('Location', 'bestoutside');
xlabel('f / f_{0}');
ylabel('Z_{in} / \Omega');
title(['Semi-Infinite Slot Impedance @ free space, d = \lambda_{0} / 4, w_{s} = \lambda_{0} / 50, ' ...
    '\Delta = \lambda_{0} / 40']);

%% VOLTAGE AS A FUNCTION OF X
x = (-1 : 0.001 : 1) * d;
v = NaN(1, length(x));
for x_idx = 1 : 1 : length(x)
    v(x_idx) = compute_v(x(x_idx), f0, ws, 1, 1, dw, d);
end

figure('Name', 'Semi-Finite Slot Normalized Magnetic Current', 'Position', [150 150 700 400]);

plot(x / d, real(v), 'LineWidth', 1.5, 'DisplayName', '\Re\{v\}');
hold on;
plot(x / d, imag(v), 'LineWidth', 1.5, 'DisplayName', '\Im\{v\}');

grid on;
% ylim([-50 150]);
legend('Location', 'bestoutside');
xlabel('x / d');
ylabel('v / V');
title(['Semi-Infinite Slot Normalized Magnetic Current @ free space, d = \lambda_{0} / 4, ' ...
    'w_{s} = \lambda_{0} / 50, \Delta = \lambda_{0} / 40']);

%% INPUT IMPEDANCE OF SEMI-INFINITE SLOT IN FREE SPACE
% NOT IMPLEMENTED (MUST BE FIXED AND COMPLETED)
% f = linspace(0.5, 2, 101) * f0;
% Zin = NaN(1, length(f));
% for f_idx = 1 : 1 : length(f)
%     Zin(f_idx) = compute_Zin_single_finite(f(f_idx), ws, 1, 4, dw, d);
% end
% 
% figure('Name', 'Semi-Finite Slot Impedance', 'Position', [150 150 700 400]);
% 
% plot(f / f0, real(Zin), 'LineWidth', 1.5, 'DisplayName', '\Re\{Z_{in}\}');
% hold on;
% plot(f / f0, imag(Zin), 'LineWidth', 1.5, 'DisplayName', '\Im\{Z_{in}\}');
% 
% grid on;
% ylim([-50 150]);
% legend('Location', 'bestoutside');
% xlabel('f / f_{0}');
% ylabel('Z_{in} / \Omega');
% title(['Semi-Infinite Slot Impedance @ \epsilon_{r,2} = 4, \epsilon_{r,1} = 1, d = \lambda_{0} / 4, ' ...
%     'w_{s} = \lambda_{0} / 50, \Delta = \lambda_{0} / 40']);

%% FUNCTIONS
% NOT IMPLEMENTED (IN CASE OF DIELECTRIC SLAB)
function Zin = compute_Zin_dielectric(f0, ws, er, dw, d, h)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    wlen = c0 / f0;
    k0 = 2 * pi / wlen;
end

% NOT IMPLEMENTED (IN CASE OF DIELECTRIC SLAB)
function D = compute_D(kx, f0, ws, er, h)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');

    eta_d = eta0 / sqrt(er);
    
    wlen = c0 / f0;
    k0 = 2 * pi / wlen;
    kd = k0 * sqrt(er);

    I1 = @(Z0, Z1) 1 ./ (Z0 + Z1);
    Zeq = @(Z0, Zd) Zd .* (Z0 + 1j * Zd * tan(kd * h)) ./ (Zd + 1j * Z0 * tan(kd * h));

    eta0_tm = @(ky) eta0 * sqrt(k0 ^ 2 - kx ^ 2 - ky .^ 2) / k0;
    eta0_te = @(ky) eta0 * k0 / sqrt(k0 ^ 2 - kx ^ 2 - ky .^ 2);
    eta_d_tm = @(ky) eta_d * sqrt(kd ^ 2 - kx ^ 2 - ky .^ 2) / kd;
    eta_d_te = @(ky) eta_d * kd / sqrt(kd ^ 2 - kx ^ 2 - ky .^ 2);

    itm = @(ky) I1(eta0_tm(kx, ky), Zeq(eta0_tm(kx, ky), eta_d_tm(kx, ky)));
    ite = @(ky) I1(eta0_te(kx, ky), Zeq(eta0_te(kx, ky), eta_d_te(kx, ky)));

    % GREEN'S FUNCTION IS WRONG, MUST INCLUDE MEDIUM 2 FROM MFIE
    ghm = @(ky) - ( (kx ^ 2) * ite(kx, ky) + (ky .^ 2) * itm(kx, ky) ) ./ (kx ^ 2 + ky .^ 2);
    integrand = @(ky) ghm(ky) .* besselj(0, ky * ws / 2);

    start_int = - 1e6;
    end_int = 1e6;

    D = integral(integrand, start_int, end_int);
    D = D / (2 * pi);
end

% MUST BE CHECKED, VERIFIED AND FIXED
function v = compute_v(x, f0, ws, er1, er2, dw, d)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    wlen = c0 / f0;
    k0 = 2 * pi / wlen;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    const = 1 / (2 * k0 * eta0);
    g = (5 / 3) * sqrt(ws * wlen);
    ig = compute_ig(f0, ws, er1, er2, dw, d);
    delta = 0.01;

    D1 = @(kx) (k1 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4);
    D2 = @(kx) (k2 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4);
    D = @(kx) const * ( D1(kx) + D2(kx) );
    Fd = @(kx) sinc(kx * dw / (2 * pi));
    Fg = @(kx) exp(1j * kx * g / 2) .* (besselj(0, kx * g / 2) - 1j * StruveH0(kx * g / 2) - (2 / pi) * sinc(kx * g / (4 * pi)) .* exp(- 1j * kx * g / 4));
    integrand = @(kx) Fd(kx) + ig * Fg(kx) .* exp(- 1j * kx * x) ./ D(kx);

    start_int = - 50 * k0 - 1j * delta;
    end_int = 50 * k0 + 1j * delta;

    % Try Dunja's waypoints
    waypts = [-(1 + 1j) (1 + 1j)] * delta;

    v = integral(integrand, start_int, end_int, 'Waypoints', waypts);
    v = v / (2 * pi);
end

function ig = compute_ig(f0, ws, er1, er2, dw, d)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    wlen = c0 / f0;
    k0 = 2 * pi / wlen;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    const = 1 / (2 * k0 * eta0);
    g = (5 / 3) * sqrt(ws * wlen);
    delta = 0.01;

    D1 = @(kx) (k1 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4);
    D2 = @(kx) (k2 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4);
    D = @(kx) const * ( D1(kx) + D2(kx) );
    Fd = @(kx) sinc(kx * dw / (2 * pi));
    Fg = @(kx) exp(1j * kx * g / 2) .* (besselj(0, kx * g / 2) - 1j * StruveH0(kx * g / 2) - (2 / pi) * sinc(kx * g / (4 * pi)) .* exp(- 1j * kx * g / 4));
    integrand_gd = @(kx) Fd(kx) .* Fg(- kx) .* exp(1j * kx * d) ./ D(kx);
    integrand_gg = @(kx) Fg(kx) .* Fg(- kx) ./ D(kx);
    
    % Needs larger integration domain, due to not convergence for small delta gap
    start_int = - 100 * k0 - 1j * delta;
    end_int = 100 * k0 + 1j * delta;

    % Try Dunja's waypoints
    waypts = [-(1 + 1j) (1 + 1j)] * delta;

    Zgd = integral(integrand_gd, start_int, end_int, 'Waypoints', waypts);
    Zgd = Zgd / (2 * pi);

    Zgg = integral(integrand_gg, start_int, end_int, 'Waypoints', waypts);
    Zgg = Zgg / (2 * pi);

    ig = - Zgd / Zgg;
end

function Zin = compute_Zin_single_finite(f0, ws, er1, er2, dw, d)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    wlen = c0 / f0;
    k0 = 2 * pi / wlen;
    k1 = k0 * sqrt(er1);
    k2 = k0 * sqrt(er2);
    
    const = 1 / (2 * k0 * eta0);
    g = (5 / 3) * sqrt(ws * wlen);
    delta = 0.01;

    D1 = @(kx) (k1 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k1 ^ 2 + kx .^ 2) / 4);
    D2 = @(kx) (k2 ^ 2 - kx .^ 2) .* besselj(0, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4) .* besselh(0, 2, - 1j * ws * sqrt(- k2 ^ 2 + kx .^ 2) / 4);
    D = @(kx) const * ( D1(kx) + D2(kx) );
    Fd = @(kx) sinc(kx * dw / (2 * pi));
    Fg = @(kx) exp(1j * kx * g / 2) .* (besselj(0, kx * g / 2) - 1j * StruveH0(kx * g / 2) - (2 / pi) * sinc(kx * g / (4 * pi)) .* exp(- 1j * kx * g / 4));
    integrand_dd = @(kx) Fd(kx) .* Fd(- kx) ./ D(kx);
    integrand_dg = @(kx) Fg(kx) .* Fd(- kx) .* exp(- 1j * kx * d) ./ D(kx);
    integrand_gd = @(kx) Fd(kx) .* Fg(- kx) .* exp(1j * kx * d) ./ D(kx);
    integrand_gg = @(kx) Fg(kx) .* Fg(- kx) ./ D(kx);
    
    % Needs larger integration domain, due to not convergence for small delta gap
    start_int = - 100 * k0 - 1j * delta;
    end_int = 100 * k0 + 1j * delta;

    % Try Dunja's waypoints
    waypts = [-(1 + 1j) (1 + 1j)] * delta;

    Zdd = integral(integrand_dd, start_int, end_int, 'Waypoints', waypts);
    Zdd = Zdd / (2 * pi);

    Zdg = integral(integrand_dg, start_int, end_int, 'Waypoints', waypts);
    Zdg = Zdg / (2 * pi);

    Zgd = integral(integrand_gd, start_int, end_int, 'Waypoints', waypts);
    Zgd = Zgd / (2 * pi);

    Zgg = integral(integrand_gg, start_int, end_int, 'Waypoints', waypts);
    Zgg = Zgg / (2 * pi);

    Zin = Zdd - Zdg * Zgd / Zgg;
end

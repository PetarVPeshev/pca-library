close all;
clear;
clc;

addpath('..\slots', '..\utils');

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
dx = 0.5 * 1e-6;
x  = (100 : 100 : 300) * 1e-6;
% NUMBER OF POINTS
Nf = length(f);
Nt = length(t);
Nx = length(x);

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% WEIGHT
% frequency domain
W  = NaN(1, Nf);
Wm = W;
for f_idx = 1 : 1 : Nf
    W(f_idx)  = slot.compute_zin(f(f_idx));
    Wm(f_idx) = evaluate_Zx(0, f(f_idx), slot, dx);
end
W  = (1 ./ W) .^ 2;
Wm = (1 ./ Wm) .^ 2;

% time domain
w  = 2 * real(eval_IFT(t, f, W));
wx = 2 * real(eval_IFT(t, f, Wm));

figure('Position', [680 558 700 420]);

plot(t * 1e12, w, 'Color', "#0072BD", 'LineWidth', 1.5, 'DisplayName', 'w');
hold on;
plot(t * 1e12, wx, '--', 'Color', "#EDB120", 'LineWidth', 1.5, 'DisplayName', 'w_{x}');

box off;
grid on;
ylim([min(w) max(w)]);
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('w [S^{2}]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
    num2str(dx * 1e6) ' \mum']);

%% IMPEDANCE
Zin = NaN(1, Nf);
Zm  = NaN(Nx, Nf);
for f_idx = 1 : 1 : Nf
    Zin(f_idx)   = slot.compute_zin(f(f_idx));
    Zm(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx);
end

%% ADMITTANCE
Hin = 1 ./ Zin;
Hm  = Zm .* Wm;

% time domain
hin = 2 * real(eval_IFT(t, f, Hin));
hm  = NaN(Nx, Nt);
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

box off;
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

    D  = @(kx) slot.compute_D(kx, f, 'TopSheet');
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

function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

% close all;
clear;
clc;

addpath('..', '..\slots', '..\utils');

% er_up : medium 2, er_dn : medium 1
config_slot  = struct('er_up', 11.7, 'er_dn', 1);
config_laser = struct('wlen', 780 * 1e-9, 'T', 12.5 * 1e-9, 'tau_p', 100 * 1e-15, 'R_3db', 5 * 1e-6);
config_pcm   = struct('wz', 2 * 1e-6, 'er', 12.96, 'tau_rec', 300 * 1e-15, 'tau_s', 8.5 * 1e-15, ...
                      'me_coef', 0.067, 'alpha', 1 * 1e-6);

% Plot time indecies
t_plt = [2066 2501 3001 4001 5001];

% Plot x indecies
x_plt = [401 501 601 701];

%% PARAMETERS
f = (0.05 : 0.005 : 2) * 1e12;
t = linspace(-2, 5, 7001) * 1e-12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% BIAS
Vb = 30;
% OPTICAL POWER
P = 5 * 1e-3;
% POINTS ALONG SLOT
dx = 1 * 1e-6;
x  = (- 400 * 1e-6 : dx : 400 * 1e-6);
% NUMBER OF POINTS
Nf = length(f);
Nt = length(t);
Nx = length(x);
% TIME STEP
dt = t(2) - t(1);
% ADMITTANCE TIME VECTOR
t_ga = (0 : 1 : Nt - 1) * dt;

%% CREATE SLOT, LASER, PCM, AND TIME-STEP ALGORITHM OBJECTS
[slot, laser, pcm] = create_objects(config_slot, config_laser, config_pcm, 'd_gap', d_gap, 'ws', ws, 'P', P);
time_step          = TimeStepAlgorithmWeight(t);

%% INPUT TIME-DOMAIN ADMITTANCE
% Frequency-domain impedance
Zin = NaN(1, Nf);
for f_idx = 1 : 1 : Nf
    Zin(f_idx) = slot.compute_zin(f(f_idx));
end

% Time-domain admittance
gin = 2 * real(eval_IFT(t_ga, f, 1 ./ Zin));
w   = 2 * real(eval_IFT(t_ga, f, (1 ./ Zin) .^ 2));

%% CREATE PHOTO-CONDUCTIVE ANTENNA OBJECT
pca = PhotoConductiveAntenna(laser, pcm, Vb, gin, 't_vec', t, 'eta_opt', 1);

%% RADIATING CURRENT
% Set time-step algorithm parameters
time_step.K       = pca.K;
time_step.Vb      = Vb;
time_step.ga      = gin;
time_step.w       = w;
time_step.tau_c   = pcm.tau_rec;
time_step.tau_s   = pcm.tau_s;
time_step.sigma_t = laser.sigma_t;

% Current
i     = NaN(1, Nt);
i_int = i;

for m = 1 : 1 : Nt
    [~, ~, i(m), i_int(m)] = step(time_step);
end

i_impr = time_step.i_impr;

% Plot currents
Title = ['@ w_{y} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P_{opt} = ' ...
         num2str(P * 1e3) ' mW, V_{b} = ' num2str(Vb) ' V'];

figure('Position', [680 558 650 370]);

plot(t * 1e12, i, 'LineWidth', 1.5, 'DisplayName', 'i');
hold on;
plot(t * 1e12, i_impr, 'LineWidth', 1.5, 'DisplayName', 'i_{impr}');
hold on;
plot(t * 1e12, i_int, 'LineWidth', 1.5, 'DisplayName', 'i_{int}');

box off;
grid on;
legend('location', 'bestoutside');
xlim([-0.5 2.5]);

xlabel('t [ps]');
ylabel('[A]');
title(Title);

%% MUTUAL TIME-DOMAIN ADMITTANCE
% Frequency-domain mutual impedance
Zm = NaN(Nx, Nf);
for f_idx = 1 : 1 : Nf
    Zm(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx);
end

% Time-domain mutual admittance
gx = NaN(Nx, Nt);
for x_idx = 1 : 1 : Nx
    gx(x_idx, :) = 2 * real(eval_IFT(t_ga, f, 1 ./ Zm(x_idx, :)));
end

% Time-domain mutual impedance
zx = NaN(Nx, Nt);
for x_idx = 1 : 1 : Nx
    zx(x_idx, :) = 2 * real(eval_IFT(t_ga, f, Zm(x_idx, :)));
end

%% MUTUAL WEIGHT AND WEIGHTED ADMITTANCE
% Frequency-domain
x0_idx = find(x == 0, 1);
Wx     = (1 ./ Zm(x0_idx, :)) .^ 2;
H      = Zm .* Wx;

% Time-domain
wx = 2 * real(eval_IFT(t_ga, f, Wx));
h  = NaN(Nx, Nt);
for x_idx = 1 : 1 : Nx
    h(x_idx, :) = 2 * real(eval_IFT(t_ga, f, H(x_idx, :)));
end

%% MUTUAL IMPEDANCE, ADMITTANCE, AND WEIGHT PLOTS
Position = [680 558 650 370];

% Frequency-domain plots
Zm_plt = NaN(length(x_plt), Nf);
Legend = cell(1, length(x_plt));
for idx = 1 : 1 : length(x_plt)
    Zm_plt(idx, :) = Zm(x_plt(idx), :);
    Legend(idx)    = cellstr(['x = ' num2str(x(x_plt(idx)) * 1e6) ' \mum']);
end
Gx_plt = 1 ./ Zm_plt;

XLim = [0 2];

plot_many(f * 1e-12, real(Zm_plt), Legend, 'f [THz]', 'Z_{x} [\Omega]', Title, XLim, Position);
plot_many(f * 1e-12, real(Gx_plt), Legend, 'f [THz]', 'G_{x} [S]', Title, XLim, Position);

% Time-domain plots
gx_plt = NaN(length(x_plt), Nt);
zx_plt = gx_plt;
h_plt  = gx_plt;
for idx = 1 : 1 : length(x_plt)
    gx_plt(idx, :) = gx(x_plt(idx), :);
    zx_plt(idx, :) = zx(x_plt(idx), :);
    h_plt(idx, :)  = h(x_plt(idx), :);
end

XLim   = [min(t_ga) max(t_ga)] * 1e12;

plot_many(t_ga * 1e12, zx_plt, Legend, 't [ps]', 'z_{x} [\Omega]', Title, XLim, Position);
plot_many(t_ga * 1e12, gx_plt, Legend, 't [ps]', 'g_{x} [S]', Title, XLim, Position);
plot_many(t_ga * 1e12, h_plt, Legend, 't [ps]', 'h_{x} [S]', Title, XLim, Position);
plot_many(t_ga * 1e12, wx, "w_{x}", 't [ps]', 'w_{x} [S^{2}]', Title, XLim, Position);

%% VOLTAGE ALONG SLOT
vx    = evaluate_vx_weight(i, h, w);
vx_zm = evaluate_vx_impedance(dt, i, zx);

% Plot comparison
LineLegend = ["w/ weight", "w/o weight"];
YLabel     = 'v_{x} [V]';
Position   = [680 558 700 420];

vx_plt      = NaN(length(t_plt), Nx);
vx_zm_plt   = vx_plt;
ColorLegend = cell(1, length(t_plt));
for idx = 1 : 1 : length(t_plt)
    vx_plt(idx, :)    = vx(:, t_plt(idx))';
    vx_zm_plt(idx, :) = vx_zm(:, t_plt(idx))';
    ColorLegend(idx)  = cellstr(['t = ' num2str(t(t_plt(idx)) * 1e12) ' ps']);
end

XLim   = [-200 200];
XLabel = 'x [\mum]';
plot_comparison(x * 1e6, vx_plt, vx_zm_plt, ColorLegend, LineLegend, 'YLabel', YLabel, 'XLabel', XLabel, ...
                'Title', Title, 'XLim', XLim, 'Position', Position);
grid on;

vx_plt      = NaN(length(x_plt), Nt);
vx_zm_plt   = vx_plt;
ColorLegend = cell(1, length(x_plt));
for idx = 1 : 1 : length(x_plt)
    vx_plt(idx, :)    = vx(x_plt(idx), :);
    vx_zm_plt(idx, :) = vx_zm(x_plt(idx), :);
    ColorLegend(idx)  = cellstr(['x = ' num2str(x(x_plt(idx)) * 1e6) ' \mum']);
end

XLim   = [-0.5 5];
XLabel = 't [ps]';
plot_comparison(t * 1e12, vx_plt, vx_zm_plt, ColorLegend, LineLegend, 'YLabel', YLabel, 'XLabel', XLabel, ...
                'Title', Title, 'XLim', XLim, 'Position', Position);
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

function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

function vx = evaluate_vx(dt, i, gm)
    % row - x, col - t
    Nt  = length(i);
    den = dt * gm(:, 1);
    vx  = NaN(size(gm));
    
    % First time index
    vx(:, 1) = i(1) ./ den;
    
    % All other time indecies
    for m = 2 : 1 : Nt
        conv_part = dt * sum(vx(:, 1 : m - 1) .* fliplr(gm(:, 2 : m)), 2);
        vx(:, m)  = (i(m) - conv_part) ./ den;
    end
end

function vx = evaluate_vx_weight(i, h, w)
    % row - x, col - t
    Nt = length(i);
    vx = NaN(size(h));
    
    % First time index
    vx(:, 1) = sum(i(:, 1) .* fliplr(h(:, 1)), 2) ./ w(:, 1);
    
    % All other time indecies
    for m = 2 : 1 : Nt
        conv_part_1 = sum(i(:, 1 : m) .* fliplr(h(:, 1 : m)), 2);
        conv_part_2 = sum(vx(:, 1 : m - 1) .* fliplr(w(:, 2 : m)), 2);
        vx(:, m)  = (conv_part_1 - conv_part_2) ./ w(:, 1);
    end
end

function vx = evaluate_vx_impedance(dt, i, zm)
    Nt = length(i);
    vx = NaN(size(zm));

    for m = 1 : 1 : Nt 
        vx(:, m) = dt * sum(i(1 : m) .* fliplr(zm(:, 1 : m)), 2);
    end
end

% FUNCTIONS USED TO CREATE OBJECTS
function [slot, laser, pcm] = create_objects(slot_config, laser_config, pcm_config, params)
    arguments
        slot_config  (1,1) struct
        laser_config (1,1) struct
        pcm_config   (1,1) struct
        params.d_gap (1,1) double
        params.ws    (1,1) double
        params.P     (1,1) double
    end

    slot  = SlotInDielectrics(params.d_gap, params.ws, slot_config.er_up, slot_config.er_dn);
    laser = Laser(laser_config.wlen, laser_config.T, params.P, 'tau_p', laser_config.tau_p, ...
                  'R_3db', laser_config.R_3db);
    pcm   = PhotoConductor([params.ws params.d_gap pcm_config.wz], pcm_config.er, ...
                           'tau_rec', pcm_config.tau_rec, 'tau_s', pcm_config.tau_s, ...
                           'me_coef', pcm_config.me_coef, 'absorp_len', pcm_config.alpha);
end

% PLOT FUNCTIONS
function fig = plot_comparison(x, y1, y2, ColorLegend, LineLegend, options)
    arguments
        x
        y1
        y2
        ColorLegend
        LineLegend
        options.XLabel
        options.YLabel
        options.Title
        options.XLim (1,:) double = NaN
        options.YLim (1,:) double = NaN
        options.Position (1,:) double = NaN
    end

    Colors = ["#0072BD", "#A2142F", "#EDB120", "#7E2F8E", "#77AC30"];
    row_num = size(y1, 1);

    fig = figure();
    if ~isnan(options.Position)
        fig.Position = options.Position;
    end
    
    AX1 = axes(fig);
    L = gobjects(1, row_num);
    for row = 1 : 1 : row_num
        L(row) = plot(x, y1(row, :), 'LineWidth', 1.5, ...
                      'Color', Colors(row), 'Parent', AX1);
        hold on;
                 plot(x, y2(row, :), '--', 'LineWidth', 1.5, ...
                      'Color', Colors(row), 'Parent', AX1);
        hold on;
    end
    
    set(AX1, 'Box', 'off');
    if ~isnan(options.XLim)
        xlim(options.XLim);
    end
    if ~isnan(options.YLim)
        ylim(options.YLim);
    end
    LGND1 = legend(AX1, L, ColorLegend, 'location', 'bestoutside');
    hold on;
    
    AX2 = copyobj(AX1, gcf);
    delete( get(AX2, 'Children') );
    H1 = plot(NaN, 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
    H2 = plot(NaN, '--', 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);

    set(AX2, 'Color', 'none', 'XTick', [], 'Visible', 'off');
    LGND2 = legend([H1 H2], LineLegend, 'location', 'bestoutside');
    LGND2_bottom = LGND1.Position(2) - LGND2.Position(4) - 0.02;
    LGND2.Position = [LGND1.Position(1) LGND2_bottom LGND2.Position(3) LGND2.Position(4)];
    LGND2.Box = 'off';
    LGND2.Color = 'none';
    hold off;
    
    xlabel(options.XLabel);
    ylabel(options.YLabel);
    title(options.Title);
end

function fig = plot_many(x, y, Legend, XLabel, YLabel, Title, XLim, Position)
    Colors = ["#0072BD", "#A2142F", "#EDB120", "#7E2F8E", "#77AC30"];
    row_num = size(y, 1);

    fig = figure('Position', Position);
    for row = 1 : 1 : row_num
        plot(x, y(row, :), 'LineWidth', 1.5, 'Color', Colors(row), 'DisplayName', Legend{row});
        hold on;
    end
    hold off;

    box off;
    grid on;
    if row ~= 1
        legend('location', 'bestoutside');
    end
    xlim(XLim);

    xlabel(XLabel);
    ylabel(YLabel);
    title(Title);
end

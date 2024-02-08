close all;
clear;
clc;

addpath('..', '..\slots', '..\utils');

% er_up : medium 2, er_dn : medium 1
config_slot  = struct('er_up', 11.7, 'er_dn', 1);
config_laser = struct('wlen', 780 * 1e-9, 'T', 12.5 * 1e-9, 'tau_p', 100 * 1e-15, 'R_3db', 5 * 1e-6);
config_pcm   = struct('wz', 2 * 1e-6, 'er', 12.96, 'tau_rec', 300 * 1e-15, 'tau_s', 8.5 * 1e-15, ...
                      'me_coef', 0.067, 'alpha', 1 * 1e-6);

% Plot time indecies
t_plt = [3650 8000 13000 23001 33001];

% Plot x indecies
x_plt = [1 2 3 4];

%% PARAMETERS
f = (0.05 : 0.005 : 2) * 1e12;
t = -2.9995e-13 : 1e-16 : 1.2e-11;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% BIAS
Vb = 30;
% OPTICAL POWER
P = 5 * 1e-3;
% POINTS ALONG SLOT
dx = 9.9875e-07;
x  = [0 99.8752 199.7503 299.6255] * 1e-6;
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
time_step          = TimeStepAlgorithmParallel(t);

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
i = NaN(1, Nt);

for m = 1 : 1 : Nt
    [~, ~, i(m), ~] = step(time_step);
end

%% MUTUAL TIME-DOMAIN ADMITTANCE
% Frequency-domain mutual impedance
Zm = NaN(Nx, Nf);
for f_idx = 1 : 1 : Nf
    Zm(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx);
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
h  = NaN(Nx, Nt);
for x_idx = 1 : 1 : Nx
    h(x_idx, :) = 2 * real(eval_IFT(t_ga, f, H(x_idx, :)));
end

%% VOLTAGE ALONG SLOT
vx = evaluate_vx_weight(i, h, w);

% Plot comparison
Title    = ['@ w_{y} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P_{opt} = ' ...
            num2str(P * 1e3) ' mW, V_{b} = ' num2str(Vb) ' V'];
YLabel   = 'v_{x} [V]';
Position = [680 558 700 420];

vx_plt      = NaN(length(t_plt), Nx);
ColorLegend = cell(1, length(t_plt));
for idx = 1 : 1 : length(t_plt)
    vx_plt(idx, :)    = vx(:, t_plt(idx))';
    ColorLegend(idx)  = cellstr(['t = ' num2str(t(t_plt(idx)) * 1e12) ' ps']);
end
vx_plt = [fliplr(vx_plt(:, 2 : end)) vx_plt];
xx_plt = [- fliplr(x(:, 2 : end)) x];

XLim   = [-400 400];
XLabel = 'x [\mum]';
plot_comparison(xx_plt * 1e6, vx_plt, ColorLegend, 'YLabel', YLabel, 'XLabel', XLabel, ...
                'Title', Title, 'XLim', XLim, 'Position', Position);
grid on;

vx_plt      = NaN(length(x_plt), Nt);
ColorLegend = cell(1, length(x_plt));
for idx = 1 : 1 : length(x_plt)
    vx_plt(idx, :)    = vx(x_plt(idx), :);
    ColorLegend(idx)  = cellstr(['x = ' num2str(x(x_plt(idx)) * 1e6) ' \mum']);
end

XLim   = [-0.3 5];
XLabel = 't [ps]';
plot_comparison(t * 1e12, vx_plt, ColorLegend, 'YLabel', YLabel, 'XLabel', XLabel, ...
                'Title', Title, 'XLim', XLim, 'Position', Position);
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
        vx(:, m)    = (conv_part_1 - conv_part_2) ./ w(:, 1);
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
function fig = plot_comparison(x, y, ColorLegend, options)
    arguments
        x
        y
        ColorLegend
        options.XLabel
        options.YLabel
        options.Title
        options.XLim (1,:) double = NaN
        options.YLim (1,:) double = NaN
        options.Position (1,:) double = NaN
    end

    Colors = ["#0072BD", "#A2142F", "#EDB120", "#7E2F8E", "#77AC30"];
    row_num = size(y, 1);

    fig = figure();
    if ~isnan(options.Position)
        fig.Position = options.Position;
    end
    
    AX1 = axes(fig);
    L = gobjects(1, row_num);
    for row = 1 : 1 : row_num
        L(row) = plot(x, y(row, :), 'LineWidth', 1.5, 'Color', Colors(row), 'Parent', AX1);
        hold on;
    end
    
    set(AX1, 'Box', 'off');
    if ~isnan(options.XLim)
        xlim(options.XLim);
    end
    if ~isnan(options.YLim)
        ylim(options.YLim);
    end
    legend(AX1, L, ColorLegend, 'location', 'bestoutside');
    
    xlabel(options.XLabel);
    ylabel(options.YLabel);
    title(options.Title);
end

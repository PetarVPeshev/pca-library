close all;
clear;
clc;

addpath('..', '..\slots', '..\utils');

% er_up : medium 2, er_dn : medium 1
config_slot  = struct('er_up', 11.7, 'er_dn', 1);
config_laser = struct('wlen', 780 * 1e-9, 'T', 12.5 * 1e-9, 'tau_p', 100 * 1e-15, 'R_3db', 5 * 1e-6);
config_pcm   = struct('wz', 2 * 1e-6, 'er', 12.96, 'tau_rec', 300 * 1e-15, 'tau_s', 8.5 * 1e-15, ...
                      'me_coef', 0.067, 'alpha', 1 * 1e-6);

LineLegend = {'\Re', '\Im'};
Location   = 'bestoutside';
LineWidth  = 1.5;
Color      = ["#0072BD", "#EDB120", "#77AC30", "#A2142F"];
Position   = [680 558 700 420];

%% PARAMETERS
f = linspace(eps, 2, 4001) * 1e12;
t = linspace(-2, 8, 4001) * 1e-12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws     = 10 * 1e-6;
ws_vec = (10 : - 2.5 : 2.5) * 1e-6;
% BIAS
Vb = 30;
% OPTICAL POWER
P     = 10 * 1e-3;
P_vec = (5 : 10 : 35) * 1e-3;
% NUMBER OF POINTS
Nf  = length(f);
Nt  = length(t);
Nws = length(ws_vec);
NP  = length(P_vec);
% ADMITTANCE TIME VECTOR
t_ga = (0 : 1 : length(t) - 1) * (t(2) - t(1));

%% CREATE LASER, AND TIME-STEP ALGORITHM OBJECTS
laser     = create_laser_object(config_laser, 'P', P);
time_step = TimeStepAlgorithmWeight(t);

%% INTERNAL IMPEDANCE
% In terms of slot width
Z_int       = NaN(Nws, Nf);
ColorLegend = cell(1, Nws);
for ws_idx = 1 : 1 : Nws
    % Create slot, and pcm objects
    slot = SlotInDielectrics(d_gap, ws_vec(ws_idx), config_slot.er_up, config_slot.er_dn);
    pcm  = create_pcm_object(config_pcm, 'd_gap', d_gap, 'ws', ws_vec(ws_idx));

    % Frequency-domain impedance
    Zin = NaN(1, Nf);
    for f_idx = 1 : 1 : Nf
        Zin(f_idx) = slot.compute_zin(f(f_idx));
    end

    % Time-domain admittance
    gin = 2 * real(eval_IFT(t_ga, f, 1 ./ Zin));
    w   = 2 * real(eval_IFT(t_ga, f, (1 ./ Zin) .^ 2));

    % Create photo-conductive antenna object
    pca = PhotoConductiveAntenna(laser, pcm, Vb, gin, 't_vec', t, 'eta_opt', 1);

    % Set time-step algorithm parameters
    time_step.K       = pca.K;
    time_step.Vb      = Vb;
    time_step.ga      = gin;
    time_step.w       = w;
    time_step.tau_c   = pcm.tau_rec;
    time_step.tau_s   = pcm.tau_s;
    time_step.sigma_t = laser.sigma_t;

    % Voltage and currents
    v     = NaN(1, Nt);
    i     = v;
    i_int = v;

    for m = 1 : 1 : Nt
        [v(m), ~, i(m), i_int(m)] = step(time_step);
    end

    i_impr = time_step.i_impr;

    % Plot currents
    figure('Position', [680 558 650 370]);

    plot(t * 1e12, i, 'LineWidth', LineWidth, 'DisplayName', 'i');
    hold on;
    plot(t * 1e12, i_impr, 'LineWidth', LineWidth, 'DisplayName', 'i_{impr}');
    hold on;
    plot(t * 1e12, i_int, 'LineWidth', LineWidth, 'DisplayName', 'i_{int}');

    box off;
    grid on;
    legend('location', 'bestoutside');
    xlim([-0.5 2.5]);

    xlabel('t [ps]');
    ylabel('[A]');
    title(['@ w_{y} = ' num2str(ws_vec(ws_idx) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ...
           ' \mum, P_{opt} = ' num2str(P * 1e3) ' mW, V_{b} = ' num2str(Vb) ' V']);

    % Frequency domain of transient voltage and internal current
    V     = eval_FT(t, f, v);
    I_int = eval_FT(t, f, i_int);

    % Internal impedance
    Z_int(ws_idx, :) = V ./ I_int;

    % Release and reset algorithm
    release(time_step);
    reset(time_step);

    % Legend
    ColorLegend(ws_idx) = cellstr(['w_{y} = ' num2str(ws_vec(ws_idx) * 1e6) ' \mum']);
end

XLabel = 'f [THz]';
YLabel = 'Z_{int} [\Omega]';
Title = ['@ \Delta = ' num2str(d_gap * 1e6) ' \mum, P_{opt} = ' num2str(P * 1e3) ' mW, V_{b} = ' ...
         num2str(Vb) ' V'];

plot_two(f * 1e-12, real(Z_int), imag(Z_int), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'Position', Position);
grid on;

% Create slot, and pcm objects
slot = SlotInDielectrics(d_gap, ws, config_slot.er_up, config_slot.er_dn);
pcm  = create_pcm_object(config_pcm, 'd_gap', d_gap, 'ws', ws);

% Frequency-domain admittance
Zin = NaN(1, Nf);
for f_idx = 1 : 1 : Nf
    Zin(f_idx) = slot.compute_zin(f(f_idx));
end

% Time-domain admittance
gin = 2 * real(eval_IFT(t_ga, f, 1 ./ Zin));
w   = 2 * real(eval_IFT(t_ga, f, (1 ./ Zin) .^ 2));

% Set time-step algorithm parameters
time_step.Vb      = Vb;
time_step.ga      = gin;
time_step.w       = w;
time_step.tau_c   = pcm.tau_rec;
time_step.tau_s   = pcm.tau_s;
time_step.sigma_t = laser.sigma_t;

% Create pca object
pca = PhotoConductiveAntenna(laser, pcm, Vb, gin, 't_vec', t, 'eta_opt', 1);

% In terms of optical power
Z_int       = NaN(NP, Nf);
ColorLegend = cell(1, NP);
for P_idx = 1 : 1 : NP
    % Set laser power
    laser.P = P_vec(P_idx);

    % Set time-step algorithm K constant
    time_step.K = pca.K;

    % Voltage and currents
    v     = NaN(1, Nt);
    i     = v;
    i_int = v;

    for m = 1 : 1 : Nt
        [v(m), ~, i(m), i_int(m)] = step(time_step);
    end

    i_impr = time_step.i_impr;

    % Plot voltage and currents
    figure('Position', [680 558 650 370]);

    plot(t * 1e12, v, 'LineWidth', LineWidth, 'DisplayName', 'v');

    box off;
    grid on;
    xlim([-0.5 2.5]);

    xlabel('t [ps]');
    ylabel('v [V]');
    title(['@ w_{y} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P_{opt} = ' ...
           num2str(P_vec(P_idx) * 1e3) ' mW, V_{b} = ' num2str(Vb) ' V']);

    figure('Position', [680 558 650 370]);

    plot(t * 1e12, i, 'LineWidth', LineWidth, 'DisplayName', 'i');
    hold on;
    plot(t * 1e12, i_impr, 'LineWidth', LineWidth, 'DisplayName', 'i_{impr}');
    hold on;
    plot(t * 1e12, i_int, 'LineWidth', LineWidth, 'DisplayName', 'i_{int}');

    box off;
    grid on;
    legend('location', 'bestoutside');
    xlim([-0.5 2.5]);

    xlabel('t [ps]');
    ylabel('[A]');
    title(['@ w_{y} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P_{opt} = ' ...
           num2str(P_vec(P_idx) * 1e3) ' mW, V_{b} = ' num2str(Vb) ' V']);

    % Frequency domain of transient voltage and internal current
    V     = eval_FT(t, f, v);
    I_int = eval_FT(t, f, i_int);

    % Internal impedance
    Z_int(P_idx, :) = V ./ I_int;

    % Release and reset algorithm
    release(time_step);
    reset(time_step);

    % Legend
    ColorLegend(P_idx) = cellstr(['P_{opt} = ' num2str(P_vec(P_idx) * 1e3) ' mW']);
end

XLabel = 'f [THz]';
YLabel = 'Z_{int} [\Omega]';
Title = ['@ w_{y} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, V_{b} = ' ...
         num2str(Vb) ' V'];

plot_two(f * 1e-12, real(Z_int), imag(Z_int), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'Position', Position);
grid on;

%% FUNCTIONS
function Y = eval_FT(t, f, y)
    N = length(f);
    dt = t(2) - t(1);

    Y = NaN(1, N);
    for f_idx = 1 : 1 : N
        Y(f_idx) = sum(y .* exp(- 1j * 2 * pi * f(f_idx) * t)) * dt;
    end
end

function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

% FUNCTIONS USED TO CREATE OBJECTS
function laser = create_laser_object(laser_config, params)
    arguments
        laser_config (1,1) struct
        params.P     (1,1) double
    end

    laser = Laser(laser_config.wlen, laser_config.T, params.P, 'tau_p', laser_config.tau_p, ...
                  'R_3db', laser_config.R_3db);
end

function pcm = create_pcm_object(pcm_config, params)
    arguments
        pcm_config   (1,1) struct
        params.d_gap (1,1) double
        params.ws    (1,1) double
    end

    pcm = PhotoConductor([params.ws params.d_gap pcm_config.wz], pcm_config.er, ...
                         'tau_rec', pcm_config.tau_rec, 'tau_s', pcm_config.tau_s, ...
                         'me_coef', pcm_config.me_coef, 'absorp_len', pcm_config.alpha);
end

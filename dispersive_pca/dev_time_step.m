close all;
clear;
clc;

addpath('..', '..\slots', '..\utils');

% er_up : medium 2, er_dn : medium 1
config_slot  = struct('er_up', 11.7, 'er_dn', 1);
config_laser = struct('wlen', 780 * 1e-9, 'T', 12.5 * 1e-9, 'tau_p', 100 * 1e-15, 'R_3db', 5 * 1e-6);
config_pcm   = struct('wz', 2 * 1e-6, 'er', 12.96, 'tau_rec', 300 * 1e-15, 'tau_s', 8.5 * 1e-15, ...
                      'me_coef', 0.067, 'alpha', 1 * 1e-6);

%% PARAMETERS
f = linspace(eps, 2, 4001) * 1e12;
t = linspace(-2, 8, 4001) * 1e-12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% BIAS
Vb = 30;
% OPTICAL POWER
P = 10 * 1e-3;
% NUMBER OF POINTS
Nf = length(f);
Nt = length(t);
% ADMITTANCE TIME VECTOR
t_ga = (0 : 1 : length(t) - 1) * (t(2) - t(1));

%% CREATE SLOT, LASER, AND PCM OBJECTS
[slot, laser, pcm] = create_objects(config_slot, config_laser, config_pcm, 'd_gap', d_gap, 'ws', ws, 'P', P);

%% INPUT ADMITTANCE AND WEIGHT
Zin = NaN(1, Nf);
for f_idx = 1 : 1 : Nf
    Zin(f_idx) = slot.compute_zin(f(f_idx));
end

% Time-domain admittance
gin = 2 * real(eval_IFT(t_ga, f, 1 ./ Zin));
w   = 2 * real(eval_IFT(t_ga, f, (1 ./ Zin) .^ 2));

%% CREATE PHOTO-CONDUCTIVE ANTENNA AND TIME-STEP ALGORITHM OBJECTS
pca                           = PhotoConductiveAntenna(laser, pcm, Vb, gin, 't_vec', t, 'eta_opt', 1);
[time_step, time_step_weight] = create_time_step_objects(t, gin, w, pca);

%% CURRENTS AND VOLTAGES
v       = NaN(1, Nt);
v_w     = v;
vg      = v;
vg_w    = v;
i       = v;
i_w     = v;
i_int   = v;
i_int_w = v;

for m = 1 : 1 : Nt
    [v(m), vg(m), i(m), i_int(m)]         = step(time_step);
    [v_w(m), vg_w(m), i_w(m), i_int_w(m)] = step(time_step_weight);
end

i_impr   = time_step.i_impr;
i_impr_w = time_step_weight.i_impr;

figure('Position', [250 250 1000 570]);

subplot(4, 1, 1);
plot(t * 1e12, v, 'LineWidth', 1.5, 'DisplayName', 'w/o weight');
hold on;
plot(t * 1e12, v_w, '--', 'Color', "#A2142F", 'LineWidth', 1.5, 'DisplayName', 'w/ weight');

box off;
grid on;
legend('location', 'bestoutside');
xlim([-1 5]);
ylim([0 30]);
ylabel('v [V]');

subplot(4, 1, 2);
plot(t * 1e12, vg, 'LineWidth', 1.5);
hold on;
plot(t * 1e12, vg_w, '--', 'Color', "#A2142F", 'LineWidth', 1.5);
hold on;
VBLine = yline(Vb, 'k--', 'LineWidth', 1.5);

box off;
grid on;
legend(VBLine, 'V_{b}', 'location', 'bestoutside');
xlim([-1 5]);
ylim([0 30]);
ylabel('v_{g} [V]');

subplot(4, 1, 3);
plot(t * 1e12, i, 'LineWidth', 1.5);
hold on
plot(t * 1e12, i_w, '--', 'Color', "#A2142F", 'LineWidth', 1.5);

box off
grid on;
xlim([-1 5]);
ylim([0 0.25]);
ylabel('i [A]');

subplot(4, 1, 4);
plot(t * 1e12, i_int, 'LineWidth', 1.5);
hold on
plot(t * 1e12, i_int_w, '--', 'Color', "#A2142F", 'LineWidth', 1.5);

box off;
grid on;
xlim([-1 5]);
ylim([min(i_int) 0.25]);
ylabel('i_{int} [A]');

xlabel('t [ps]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P = ' ...
         num2str(P * 1e3) ' mW'], 'FontSize', 11, 'FontWeight', 'bold');

%% VERIFICATION
I       = NaN(2, Nf);
I(1, :) = eval_FT(t, f, i);
I(2, :) = eval_FT(t, f, i_w);

VZin       = NaN(2, Nf);
VZin(1, :) = eval_FT(t, f, v) ./ Zin;
VZin(2, :) = eval_FT(t, f, v_w) ./ Zin;

Position = [680 558 700 420];
ColorLegend = {'w/o weight', 'w/ weight'};
LineLegend = {'I', 'V / Z_{in}'};
Location = 'bestoutside';
LineWidth = 1.5;
Color = ["#0072BD", "#A2142F"];
XLabel = 'f [THz]';

YLabel = '\Re\{I\} [A / Hz]';
YLim = [-5 25] * 1e-14;
Title = ['\Re @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P = ' ...
         num2str(P * 1e3) ' mW'];
plot_two(f * 1e-12, real(I), real(VZin), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'YLim', YLim, 'Title', Title, 'Position', Position);
grid on;

YLabel = '\Im\{I\} [A / Hz]';
YLim = [-1 0] * 1e-13;
Title = ['\Im @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P = ' ...
         num2str(P * 1e3) ' mW'];
plot_two(f * 1e-12, imag(I), imag(VZin), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'YLim', YLim, 'Title', Title, 'Position', Position);
grid on;

%% FUNCTIONS
function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

function Y = eval_FT(t, f, y)
    N = length(f);
    dt = t(2) - t(1);

    Y = NaN(1, N);
    for f_idx = 1 : 1 : N
        Y(f_idx) = sum(y .* exp(- 1j * 2 * pi * f(f_idx) * t)) * dt;
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

function [time_step, time_step_weight] = create_time_step_objects(t, gin, w, pca)
    time_step        = TimeStepAlgorithm(t, 'ga', gin, 'K', pca.K, 'Vb', pca.Vb, ...
                                         'sigma_t', pca.laser.sigma_t, ...
                                         'tau_c', pca.pc_material.tau_rec, ...
                                         'tau_s', pca.pc_material.tau_s);
    time_step_weight = TimeStepAlgorithmWeight(t, 'ga', gin, 'w', w, 'K', pca.K, 'Vb', pca.Vb, ...
                                               'sigma_t', pca.laser.sigma_t, ...
                                               'tau_c', pca.pc_material.tau_rec, ...
                                               'tau_s', pca.pc_material.tau_s);
end

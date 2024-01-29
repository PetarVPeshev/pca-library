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
Location = 'bestoutside';
LineWidth = 1.5;
Color = ["#0072BD", "#EDB120"];
Position = [680 558 700 420];

%% PARAMETERS
% f = 0.05e12 : 5e9 : 2e12;             % Lawrence vectors
% t = -0.3e-12 : 0.0001e-12 : 3e-12;    % Lawrence vectors
f = (0.05 : 0.005 : 2) * 1e12;
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

%% EVALUATE INPUT ADMITTANCE AND WEIGHT
Zin = NaN(1, Nf);
for f_idx = 1 : 1 : Nf
    Zin(f_idx) = slot.compute_zin(f(f_idx));
end

% Time-domain admittance
gin = 2 * real(eval_IFT(t_ga, f, 1 ./ Zin));
w   = 2 * real(eval_IFT(t_ga, f, (1 ./ Zin) .^ 2));

% Plot
% figure('Position', [50 50 700 420]);
% 
% plot(t_ga * 1e12, gin, 'LineWidth', 1.5);
% grid on;
% xlim([min(t_ga) max(t_ga)] * 1e12);
% 
% xlabel('t [ps]');
% ylabel('g_{in} [S/s]');
% title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum']);
% 
% figure('Position', [50 50 700 420]);
% 
% plot(t_ga * 1e12, w, 'LineWidth', 1.5);
% grid on;
% xlim([min(t_ga) max(t_ga)] * 1e12);
% 
% xlabel('t [ps]');
% ylabel('w [S^{2}/s]');
% title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum']);

%% CREATE PHOTO-CONDUCTIVE ANTENNA AND TIME-STEP ALGORITHM OBJECTS
pca                           = PhotoConductiveAntenna(laser, pcm, Vb, gin, 't_vec', t, 'eta_opt', 1);
[time_step, time_step_weight] = create_time_step_objects(t, gin, w, pca);

%% EVALUATE CURRENTS AND VOLTAGES
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

figure('Position', [50 50 700 420]);

plot(t * 1e12, v, 'LineWidth', 1.5, 'DisplayName', 'v, no weight');
hold on;
plot(t * 1e12, v_w, '--', 'LineWidth', 1.5, 'DisplayName', 'v, weight');

grid on;
legend('location', 'bestoutside');
xlim([min(t) max(t)] * 1e12);

xlabel('t [ps]');
ylabel('v [V]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum']);

%% FUNCTIONS
function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

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

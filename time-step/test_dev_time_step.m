close all;
clear;
clc;

run('..\dispersive_pca\coupling\env.m');
run('..\dispersive_pca\coupling\sim_params.m');

%% ANTENNA
d_feeds   = 100 * 1e-6;      % FEED DISTANCE   : 100 um (distance b/n feeds)
num_feeds = 3;               % NUMBER OF FEEDS : 3
Vb        = [30 0 0];        % BIAS VOLTAGE    : 30 V @ first feed, other feeds are unbiased
feedDelay = [0 1 2] * 1e-12; % FEED DELAY TIME : 0 ps, 1 ps, 2 ps

%% FREQUENCY POINTS
Nf     = 4001;          % FREQUENCY POINTS : 4001
f_lims = [1e9 3e12];    % FREQUENCY LIMITS : 1 MHz to 3 THz

%% SIMULATION TIME
t_lims = [-3 5] * 1e-12;  % TIME LIMITS : -1 ps to 10 ps
dt     = 0.5 * 1e-15;     % TIME STEP   : 0.5 fs

%% OBJECTS
slot  = create_slot(params_slot);
pcm   = create_pcm(params_pcm, 'ws', params_slot.ws, 'd_gap', params_slot.d_gap);
laser = create_laser(params_laser);

%% TIME STRUCT & FREQUENCY ARRAY
% Frequency
mustBeIncreasing(f_lims);
f    = linspace(f_lims(1), f_lims(2), Nf);
% Time
time = create_time_arrays(dt, t_lims, 'struct');

%% ANTENNA IMPULSE RESPONSE
% impedMatrix = create_imped_matrix(f, slot, 'DistFeeds', d_feeds, 'NumFeeds', num_feeds);
% [impulResp, weight] = create_impul_resp_matrix(time.t_res, f, impedMatrix);

%% TEST TIME STEP STATIC FUNCTIONS
K = calculate_K(laser, pcm);
% Fs = compute_Fs(t, K, sigma_t, tau_c, tau_d)
Fs = TimeStep1N.compute_Fs(time.t_sim, K, laser.sigma_t, pcm.tau_rec, feedDelay);
Fs_check = TimeStep1.compute_Fs(time.t_sim, K, laser.sigma_t, pcm.tau_rec, feedDelay);

figure('Name', 'Photo-Conductor Impulse Response', 'Position', [250 250 950 450]);

plot(time.t_sim * 1e12, Fs(1, :), 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', 'FEED 1');
hold on;
plot(time.t_sim * 1e12, Fs_check(1, :), '--', 'LineWidth', 1.5, 'Color', '#FF00FF', ...
     'DisplayName', 'FEED 1, test');
hold on;
plot(time.t_sim * 1e12, Fs(2, :), 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', 'FEED 2');
hold on;
plot(time.t_sim * 1e12, Fs_check(2, :), '--', 'LineWidth', 1.5, 'Color', '#00FF00', ...
     'DisplayName', 'FEED 1, test');
hold on;
plot(time.t_sim * 1e12, Fs(3, :), 'LineWidth', 1.5, 'Color', '#77AC30', ...
     'DisplayName', 'FEED 3');

grid on;
legend('Location', 'bestoutside');

xlabel('t [ps]');
ylabel('F_{s} [-]');
title(['@ K = ' num2str(round(K * 1e-27, 2)) '* 10^{27}, \tau_{p} = ' num2str(laser.tau_p * 1e15) ' fs, ' ...
       '\tau_{c} = ' num2str(pcm.tau_rec * 1e15) ' fs']);

%% TIME STEP OBJECT
% timeStep = TimeStep1N(time.t_sim);
% 
% timeStep.numFeeds = num_feeds;
% timeStep.K        = calculate_K(laser, pcm);
% timeStep.Vb       = Vb;
% timeStep.w        = weight.time_domain;
% timeStep.h        = impulResp.time_domain;
% timeStep.tau_c    = pcm.tau_rec;
% timeStep.tau_s    = pcm.tau_s;
% timeStep.tau_d    = feedDelay;
% timeStep.sigma_t  = laser.sigma_t;
% 
% [v, i] = timeStep.simulate();

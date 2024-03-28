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
impedMatrix = create_imped_matrix(f, slot, 'DistFeeds', d_feeds, 'NumFeeds', num_feeds);
[impulResp, weight] = create_impul_resp_matrix(time.t_res, f, impedMatrix);

%% TIME STEP OBJECT
timeStep = TimeStepByImpressedCurrent(time.t_sim);

timeStep.numFeeds     = num_feeds;
timeStep.constK       = calculate_K(laser, pcm);
timeStep.bias         = Vb;
timeStep.weight       = weight.time_domain;
timeStep.impulResp    = impulResp.time_domain;
timeStep.recTime      = pcm.tau_rec;
timeStep.scatTime     = pcm.tau_s;
timeStep.feedDelay    = feedDelay;
timeStep.laserTimeStd = laser.sigma_t;

[voltages, currents] = timeStep.simulate();

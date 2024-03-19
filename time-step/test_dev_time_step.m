close all;
clear;
clc;

run('..\dispersive_pca\coupling\env.m');
run('..\dispersive_pca\coupling\sim_params.m');

%% ANTENNA
d_feeds = 100 * 1e-6;  % FEED DISTANCE : 100 um (distance b/n feeds)

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
Z = eval_Z(f, slot, d_feeds);

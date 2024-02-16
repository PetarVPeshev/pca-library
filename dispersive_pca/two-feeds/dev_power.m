close all;
clear;
clc;

env;
sim_params;

%% PARAMETERS
f       = linspace(50 * 1e9, 2 * 1e12, 4001);
dt      = 1 * 1e-15;
[t, tt] = create_t_vec(dt, -0.5 * 1e-12, 5 * 1e-12);
% FEED
d_gap = 10 * 1e-6;
d     = 200 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% LASER POWER
P = 5 * 1e-3;
% BIAS VOLTAGE
Vb = 30;
% NUMBER OF POINTS
Nf  = length(f);
Nt  = length(t);

%% FUNCTIONS
% FUNCTIONS USED TO CREATE OBJECTS AND VECTORS
function slot = create_slot(d_gap, ws, params_slot)
    slot = SlotInDielectrics(d_gap, ws, params_slot.er_up, params_slot.er_dn);
end

function pcm = create_pcm(d_gap, ws, pcm_config)
    pcm = PhotoConductor([ws d_gap pcm_config.wz], pcm_config.er, ...
                         'tau_rec', pcm_config.tau_rec, 'tau_s', pcm_config.tau_s, ...
                         'me_coef', pcm_config.me_coef, 'absorp_len', pcm_config.alpha);
end

function laser = create_laser(P, laser_config)
    laser = Laser(laser_config.wlen, laser_config.T, P, 'tau_p', laser_config.tau_p, ...
                  'R_3db', laser_config.R_3db);
end

function [t, tt] = create_t_vec(dt, tstart, tend)
    t  = tstart : dt : tend;
    tt = (0 : 1 : length(t) - 1) * dt;
end

% SIMULATION PARAMETERS
% This script contains typical simulation parameters for photo-conductive antennas.
% Run this script to load default simulation parameters.

% SLOT
% Infinite slot b/n two dielectrics
% Upward dielectric   (medium 2) : 11.7 (silicon)
% Downward dielectric (medium 1) : 1 (free-space)
% Feed length                    : 10 um
% Slot width                     : 5 um
params_slot = struct('er_up', 11.7, 'er_dn', 1, 'd_gap', 10 * 1e-6, 'ws', 5 * 1e-6);

% LASER
% Gaussian laser distributed over space and time
% Wavelength                    : 780 nm
% Repetition rate               : 12.5 ns
% Average power                 : 180 mW
% Pulse half width power (FWHM) : 100 fs
% Radius half width             : 5 um
params_laser = struct('wlen', 780 * 1e-9, 'T', 12.5 * 1e-9, 'P', 180 * 1e-3, ...
                      'tau_p', 100 * 1e-15, 'R_3db', 5 * 1e-6);

% PHOTO-CONDUCTIVE MATERIAL
% Gallium-Arsenide (GaAs)
% Carrier recombination time : 300 fs
% Carrier scattering time    : 8.5 fs
% Dielectric constant        : 12.96
% Electron effective mass    : 0.067
% Absorption coefficient     : 1 um
% Height                     : 2 um
params_pcm = struct('wz', 2 * 1e-6, 'er', 12.96, 'tau_rec', 300 * 1e-15, ...
                    'tau_s', 8.5 * 1e-15, 'me_coef', 0.067, 'alpha', 1 * 1e-6);

% FREQUENCY DOMAIN
% Parameters specifying the minimum frequency pointa and frequency
% resolution
% Minimum frequency    : 50 GHz
% Frequency resolution : 
params_fd = struct('f_min', 50 * 1e9, 'df', 0);

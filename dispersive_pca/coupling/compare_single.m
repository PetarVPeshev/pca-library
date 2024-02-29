close all;
clear;
clc;

env;
sim_params;

%% PARAMETERS
% SLOT
ws     = 10 * 1e-6;          % SLOT WIDTH      : 10 um
d_gap  = 10 * 1e-6;          % FEED WIDTH      : 10 um
x_feed = -50 * 1e-6;         % FEED LOCATION   : -50 um
x_lims = [-100 100] * 1e-6;  % DISTANCE LIMITS : -200 um to 200 um
dx     = 10 * 1e-6;          % DISTANCE STEP   : 10 um
% PHOTO-CONDUCTORS
Vb    = 30;                  % BIAS VOLTAGE     : 30 V & 0 V
% LASER
P = 10 * 1e-3;               % AVERAGE POWER : 10 mW
% SIMULATION TIME
t_lims = [-1 10] * 1e-12;         % TIME LIMITS : -1 ps to 10 ps
dt     = 0.5 * 1e-15;             % TIME STEP   : 0.5 fs
t_plot = [5e-2 1 1.5 2] * 1e-12;  % PLOT TIME   : 5 fs, 1 ps, and 1.5 ps
% FREQUENCY DOMAIN
Nf     = 4001;          % FREQUENCY POINTS : 4001
f_lims = [1e9 3e12];    % FREQUENCY LIMITS : 1 MHz to 3 THz

%% FREQUENCY, TIME, & SLOT POSITION VECTORS
mustBeIncreasing(f_lims);
mustBeIncreasing(x_lims);
f            = linspace(f_lims(1), f_lims(2), Nf);
[tsim, tres] = create_time_arrays(dt, t_lims);
Nt           = length(tsim);
x            = x_lims(1) : dx : x_lims(2);
Nx           = length(x);

%% TIME INDEX PLOT
Nt_plot = length(t_plot);
t_idx   = NaN(1, Nt_plot);
for idx = 1 : Nt_plot
    t_idx(idx) = find(round(tsim * 1e12, 4) == t_plot(idx) * 1e12, 1);
end

%% FEED INDEX
x_idx = NaN(1, 2);
for idx = 1 : 2
    x_idx(idx) = find(round(x * 1e6) == round(((-1) ^ idx) * x_feed * 1e6), 1);
end

%% OBJECTS
slot  = create_slot(params_slot, 'd_gap', d_gap, 'ws', ws);
pcm   = create_pcm(params_pcm, 'd_gap', d_gap, 'ws', ws);
laser = create_laser(params_laser, 'P', P);

%% INPUT TIME-DOMAIN ADMITTANCE
% Frequency-domain impedance
Zin = slot.compute_zin(f);

% Time-domain admittance
gin = 2 * real(eval_IFT(tres, f, 1 ./ Zin));
w   = 2 * real(eval_IFT(tres, f, (1 ./ Zin) .^ 2));

%% RADIATING CURRENT
time_step          = TimeStepAlgorithmWeight(tsim);

% Set time-step algorithm parameters
time_step.K       = calculate_K(laser, pcm);
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

%% WAVE ALONG SLOT
[vx, ~, ~] = eval_vx(i, 'x', x, 't_res', tres, 'f', f, 'slot', slot, 'x_feed', x_feed);

figure('Position', [250 250 850 550]);

for idx = 1 : Nt_plot
    subplot(Nt_plot, 1, idx);
    plot(x * 1e6, vx(:, t_idx(idx))', '--', 'LineWidth', 1.5, 'Color', '#0072BD');
    hold on;
    plot(x_feed * 1e6, 0, 'm*', 'LineWidth', 3.0);
    
    grid on;
    ylim([-3 12.5]);
    yticks(-3 : 3 : 12);
    
    ylabel('v [V]');
    title(['@ t = ' num2str(t_plot(idx) * 1e12) ' ps']);
end

xlabel('x [\mum]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ...
         ' \mum, x_{F} = ' num2str(x_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], ...
         'FontSize', 11, 'FontWeight', 'bold');

figure('Position', [250 250 950 550]);

for idx = 1 : 2
    subplot(2, 1, idx);
    plot(tsim * 1e12, vx(x_idx(idx), :), 'LineWidth', 1.5, 'Color', '#0072BD', ...
         'DisplayName', 'ONE FEED');
    
    grid on;
    xlim([-0.3 2.5]);
    ylim([-3 12.5]);
    yticks(-3 : 3 : 12);
    legend('location', 'bestoutside');
    
    ylabel('v [V]');
    title(['x = ' num2str(x(x_idx(idx)) * 1e6) ' \mum']);
end

xlabel('t [ps]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, x_{F} = ' ...
        num2str(x_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], 'FontSize', 11, ...
        'FontWeight', 'bold');

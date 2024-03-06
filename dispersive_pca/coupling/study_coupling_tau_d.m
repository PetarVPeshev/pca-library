close all;
clear;
clc;

env;
sim_params;

c0 = get_phys_const('LightSpeed');

%% PARAMETERS
% SLOT
ws     = 10 * 1e-6;          % SLOT WIDTH      : 10 um
d_gap  = 10 * 1e-6;          % FEED WIDTH      : 10 um
d_feed = 100 * 1e-6;         % FEED DISTANCE   : 100 um (distance b/n feeds)
% PHOTO-CONDUCTORS
Vb    = [30; 0];            % BIAS VOLTAGE                      : 30 V & 0 V
Ntau  = 101;                % NUMBER OF EXCITATION DELAY POINTS : 101
% LASER
P = 10 * 1e-3;  % AVERAGE POWER : 10 mW
% SIMULATION TIME
t_lims = [-1 10] * 1e-12;         % TIME LIMITS : -1 ps to 10 ps
dt     = 0.5 * 1e-15;             % TIME STEP   : 0.5 fs
% FREQUENCY DOMAIN
Nf     = 4001;          % FREQUENCY POINTS : 4001
f_lims = [1e9 3e12];    % FREQUENCY LIMITS : 1 MHz to 3 THz

%% ALLOCATE SPACE FOR ENERGY VECTORS
E_in   = NaN(1, Ntau);
E1_ind = E_in;
E2_ind = E_in;
E_m    = E_in;

%% FREQUENCY, TIME, & SLOT POSITION VECTORS
mustBeIncreasing(f_lims);
f            = linspace(f_lims(1), f_lims(2), Nf);
[tsim, tres] = create_time_arrays(dt, t_lims);
Nt           = length(tsim);

%% OBJECTS
slot  = create_slot(params_slot, 'd_gap', d_gap, 'ws', ws);
pcm   = create_pcm(params_pcm, 'd_gap', d_gap, 'ws', ws);
laser = create_laser(params_laser, 'P', P);

%% EFFECTIVE PERMITTIVITY & WAVE PROPAGATION
n_eff = real(slot.find_kxp(f)) ./ (2 * pi * f / c0);
tau_w = d_feed * n_eff / c0;

%% EXCITATION DELAY @ FEED 2
tau_d = linspace(-3, 3, Ntau) * laser.tau_p + round(mean(tau_w), 14);

%% IMPEDANCE
Z = eval_Z(f, slot, d_feed);

%% ADMITTANCE
[H, W] = compute_H(Z);

%% IMPULSE RESPONSE
h = compute_h_td(tres, f, H);
w = 2 * real(eval_IFT(tres, f, W));

W_struct = struct('w', w, 'W', W);

parfor idx = 1 : Ntau
    %% TIME-STEP ALGORITHM
    time_step = TimeStep1(tsim);

    % Set time-step algorithm parameters
    time_step.K       = calculate_K(laser, pcm);
    time_step.Vb      = Vb;
    time_step.w       = w;
    time_step.h       = h;
    time_step.tau_c   = pcm.tau_rec;
    time_step.tau_s   = pcm.tau_s;
    time_step.tau_d   = [0 tau_d(idx)];
    time_step.sigma_t = laser.sigma_t;
    
    % Voltage and currents
    [v, i] = time_step.simulate();
    
    %% VOLTAGE CONTRIBUTIONS
    v1s = eval_vx(i(1, :), 'h', struct('h', permute(h(1, 1, :), [1 3 2])), 'w', W_struct);
    v2s = eval_vx(i(1, :), 'h', struct('h', permute(h(2, 1, :), [1 3 2])), 'w', W_struct);
    v1m = eval_vx(i(2, :), 'h', struct('h', permute(h(1, 2, :), [1 3 2])), 'w', W_struct);
    v2m = eval_vx(i(2, :), 'h', struct('h', permute(h(2, 2, :), [1 3 2])), 'w', W_struct);
    
    %% INSTANTANEOUS POWER COMPONENTS
    p_in    = v1s .* i(1, :);
    p1m_ind = v1m .* i(1, :);
    p2s_ind = v2s .* i(2, :);
    pm_rad  = v2m .* i(2, :);
    
    E_in(idx)   = sum(p_in, 2) * dt;
    E1_ind(idx) = sum(p1m_ind, 2) * dt;
    E2_ind(idx) = sum(p2s_ind, 2) * dt;
    E_m(idx)    = sum(pm_rad, 2) * dt;
end

%% PLOTS
% EFFECTIVE REFRACTIVE INDEX
figure('Position', [250 250 950 450]);

subplot(2, 2, 1);
plot(f * 1e-12, n_eff .^ 2, 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', '\epsilon_{r,eff}');

grid on;
ylim([6 10]);

xlabel('f [THz]');
ylabel('\epsilon_{r,eff} [-]');

subplot(2, 2, 2);
plot(f * 1e-12, n_eff, 'LineWidth', 1.5, 'Color', '#77AC30', ...
     'DisplayName', 'n_{eff}')

grid on;
ylim([2.5 3.1]);

xlabel('f [THz]');
ylabel('n_{eff} [-]');

subplot(2, 2, [3, 4]);
plot(f * 1e-12, tau_w * 1e12, 'LineWidth', 1.5, 'Color', '#7E2F8E', ...
     'DisplayName', ['\tau_{w}, \tau_{w,av} = ' num2str(round(mean(tau_w) * 1e12, 2)) ' ps']);

grid on;
ylim([0.85 1.05]);

xlabel('f [THz]');
ylabel('\tau_{w} [ps]');
title(['\tau_{w,av} = ' num2str(round(mean(tau_w) * 1e12, 2)) ' ps @ d_{F} = ' ...
       num2str(d_feed * 1e6) ' \mum']);
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum'], ...
        'FontSize', 11, 'FontWeight', 'bold');

% ENERGY COMPONENTS
figure('Position', [250 250 1250 550]);

subplot(2, 2, 1);
plot(tau_d * 1e12, round(E_in * 1e12, 3), 'LineWidth', 1.5, 'Color', '#0000FF', ...
     'DisplayName', 'E_{in}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);
ylim([0 1]);

ylabel('E [pJ]');
title('E_{in}');

subplot(2, 2, 2);
plot(tau_d * 1e12, E1_ind * 1e12, 'LineWidth', 1.5, 'Color', '#FF0000', ...
     'DisplayName', 'E_{1}^{ind}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);
ylim([-0.1 0]);

title('E_{1}^{ind}');

subplot(2, 2, 3);
plot(tau_d * 1e12, E2_ind * 1e12, 'LineWidth', 1.5, 'Color', '#00FF00', ...
     'DisplayName', 'E_{2}^{ind}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);

ylabel('E [pJ]');
xlabel('\tau_{d} [ps]');
title('E_{2}^{ind}');

subplot(2, 2, 4);
plot(tau_d * 1e12, E_m * 1e12, 'LineWidth', 1.5, 'Color', '#FF00FF', ...
     'DisplayName', 'E_{m}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);

xlabel('\tau_{d} [ps]');
title('E_{m}');
    
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{F} = ' ...
        num2str(d_feed * 1e6) ' \mum, V_{b} = ' num2str(Vb(1)) ' V, \delta_{t} = ' num2str(dt * 1e15) ...
        ' fs'], 'FontSize', 11, 'FontWeight', 'bold');

% RADIATED & DISSIPATED ENERGY
figure('Position', [250 250 1250 550]);

subplot(2, 2, 1);
plot(tau_d * 1e12, (E_in + E2_ind) * 1e12, 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', 'E_{s}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);

ylabel('E [pJ]');
title('E_{s}');

subplot(2, 2, 2);
plot(tau_d * 1e12, E_m * 1e12, 'LineWidth', 1.5, 'Color', '#77AC30', ...
     'DisplayName', 'E_{m}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);

title('E_{m}');

subplot(2, 2, 3);
plot(tau_d * 1e12, (E_in + E2_ind + E_m) * 1e12, 'LineWidth', 1.5, 'Color', '#7E2F8E', ...
     'DisplayName', 'E_{rad}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);
ylim([0.54 0.64]);

ylabel('E [pJ]');
xlabel('\tau_{d} [ps]');
title('E_{rad}');

subplot(2, 2, 4);
plot(tau_d * 1e12, (E2_ind + E_m) * 1e12, 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', 'E_{dis}');

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);

xlabel('\tau_{d} [ps]');
title('E_{dis}');
    
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{F} = ' ...
        num2str(d_feed * 1e6) ' \mum, V_{b} = ' num2str(Vb(1)) ' V, \delta_{t} = ' num2str(dt * 1e15) ...
        ' fs'], 'FontSize', 11, 'FontWeight', 'bold');

% ENERGY COMPARISSON
figure('Position', [250 250 1250 550]);

subplot(2, 1, 1);
plot(tau_d * 1e12, (E_in + E2_ind) * 1e12, 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', 'E_{s}');
hold on;
plot(tau_d * 1e12, (E_in + E2_ind + E_m) * 1e12, 'LineWidth', 1.5, 'Color', '#7E2F8E', ...
     'DisplayName', 'E_{rad}');
hold on;
xline(round(mean(tau_w) * 1e12, 2), 'k--', 'LineWidth', 1.5, ...
     'DisplayName', ['\tau_{w,av} = ' num2str(round(mean(tau_w) * 1e12, 2)) ' ps']);
hold on;
xline(tau_d(find((E_in + E2_ind) == min(E_in + E2_ind), 1)) * 1e12, 'b--', 'LineWidth', 1.5, ...
     'DisplayName', ['min E_{s}, t = ' ...
     num2str(round(tau_d(find((E_in + E2_ind) == min(E_in + E2_ind), 1)) * 1e12, 2)) ' ps']);
hold on;
xline(tau_d(find((E_in + E2_ind + E_m) == min(E_in + E2_ind + E_m), 1)) * 1e12, 'm--', 'LineWidth', 1.5, ...
     'DisplayName', ['min E_{rad}, t = ' ...
     num2str(round(tau_d(find((E_in + E2_ind + E_m) == min(E_in + E2_ind + E_m), 1)) * 1e12, 2)) ' ps']);

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);
legend('location', 'bestoutside');

ylabel('E [pJ]');

subplot(2, 1, 2);
plot(tau_d * 1e12, E_m * 1e12, 'LineWidth', 1.5, 'Color', '#77AC30', ...
     'DisplayName', 'E_{m}');
hold on;
plot(tau_d * 1e12, abs(E2_ind + E_m) * 1e12, 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', 'E_{dis}');
hold on;
xline(round(mean(tau_w) * 1e12, 2), 'k--', 'LineWidth', 1.5, ...
     'DisplayName', ['\tau_{w,av} = ' num2str(round(mean(tau_w) * 1e12, 2)) ' ps']);
hold on;
xline(tau_d(find((E_m) == max(E_m), 1)) * 1e12, 'g--', 'LineWidth', 1.5, ...
     'DisplayName', ['max E_{m}, t = ' ...
     num2str(round(tau_d(find((E_m) == max(E_m), 1)) * 1e12, 2)) ' ps']);
hold on;
xline(tau_d(find((E2_ind + E_m) == min(E2_ind + E_m), 1)) * 1e12, 'r--', 'LineWidth', 1.5, ...
     'DisplayName', ['max E_{dis}, t = ' ...
     num2str(round(tau_d(find((E2_ind + E_m) == min(E2_ind + E_m), 1)) * 1e12, 2)) ' ps']);

grid on;
xlim([tau_d(1) tau_d(end)] * 1e12);
legend('location', 'bestoutside');

ylabel('E [pJ]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{F} = ' ...
        num2str(d_feed * 1e6) ' \mum, V_{b} = ' num2str(Vb(1)) ' V, \delta_{t} = ' num2str(dt * 1e15) ...
        ' fs'], 'FontSize', 11, 'FontWeight', 'bold');

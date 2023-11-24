close all;
clear;
clc;

addpath([pwd() '\utils']);

%% DEFAULT CONFIGURATION
[laser, GaAs, pca] = get_default_config(eta_opt = 0.35);

%% ANTENNA RESPONSE
[v, vg, i_impr, i_int, i] = pca.compute_response();

%% PLOT RESPONSE
figure('Position', [250 250 1400 700]);
annotation('textbox', [.91 .60 .3 .3], 'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'String', {'----QO LINK----',['\eta_{opt} = ' num2str(pca.qo_link.eta_opt)], ...
    ['P_{L} = ' num2str(laser.P * 1e3) ' mW'], ['T_{L} = ' num2str(laser.T * 1e9) ' ns'], ...
    ['\tau_{p} = ' num2str(laser.tau_p * 1e15) ' fs'], '------GaAs------', ...
    ['\tau_{s} = ' num2str(GaAs.tau_s * 1e15) ' fs'], ['\tau_{c} = ' num2str(GaAs.tau_rec * 1e15) ' fs'], ...
    ['W_{x} = ' num2str(GaAs.dimensions.Wx * 1e6) ' \mum'], ...
    ['W_{y} = ' num2str(GaAs.dimensions.Wy * 1e6) ' \mum'], ...
    ['W_{z} = ' num2str(GaAs.dimensions.Wz * 1e6) ' \mum']})
sgtitle('Photo-Conductive Antenna Simulations', 'FontWeight', 'bold');

% LASER
subplot(3, 1, 1);
plot(pca.time_step.t_vec * 1e12, exp(- 0.5 * (pca.time_step.t_vec / laser.sigma_t) .^ 2), ...
    'LineWidth', 2.0, 'DisplayName', ['P_{L}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
hold on;
xline(- laser.tau_p * 1e12, '--', 'LineWidth', 2.0, 'Color', [1 0 1], 'DisplayName', '-\tau_{p}');
hold on;
xline(laser.tau_p * 1e12, '--', 'LineWidth', 2.0, 'Color', [1 0 1], 'DisplayName', '\tau_{p}');
grid on;
xlim([-1 2]);
legend('location', 'bestoutside');
ylabel('P_{L} / kW');

% CURRENTS
subplot(3, 1, 2);
plot(pca.time_step.t_vec * 1e12, i_impr, 'LineWidth', 2.0, ...
    'DisplayName', ['i_{impr}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
hold on;
plot(pca.time_step.t_vec * 1e12, i_int, 'LineWidth', 2.0, ...
    'DisplayName', ['i_{int}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
hold on;
plot(pca.time_step.t_vec * 1e12, i, 'LineWidth', 2.0, ...
    'DisplayName', ['i, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
grid on;
xlim([-1 2]);
ylim([0 3]);
legend('location', 'bestoutside');
xlabel('t / ps');
ylabel('i / A');

% VOLTAGES
subplot(3, 1, 3);
yline(pca.Vb, '--', 'LineWidth', 2.0, 'DisplayName', 'V_{b}');
hold on;
plot(pca.time_step.t_vec * 1e12, v, 'LineWidth', 2.0, ...
    'DisplayName', ['v, \delta_t = ' num2str(pca.time_step.dt * 1e15) ' fs']);
hold on;
plot(pca.time_step.t_vec * 1e12, vg, 'LineWidth', 2.0, ...
    'DisplayName', ['v_{g}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
grid on;
xlim([-1 2]);
ylim([0 35]);
legend('location', 'bestoutside');
xlabel('t / ps');
ylabel('v / V');

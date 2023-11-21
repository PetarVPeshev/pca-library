close all;
clear;
clc;

[laser, photoconductor] = get_default_config();
march_on = TimeStepAlgorithm(laser.get_t_vec(laser.tau_p / 100));
laser.plot_time_envelope();
% eff = Efficiency(laser, photoconductor, eff_opt = 0.35);

antenna = PhotoConductiveAntenna(laser, photoconductor, march_on, eff_opt = 0.35);


%% PARAMETERS
dt = laser.pulse_half_pwr / 1000;
t_vec = laser.get_time_vec(dt, max_t = 15 * laser.pulse_half_pwr);
k = antenna.k_const;
sigma_t = laser.temporal_std;
tau_rec = photoconductor.rec_time;
tau_s = photoconductor.scat_time;
Vb = 30;

%% h_{m}
% v = NaN(1, length(t_vec));
% v(1) = 0;

i_impr = compute_i_impr(t_vec, k, Vb, tau_s, tau_rec, sigma_t);

figure('Position', [250 250 700 400]);
plot(t_vec * 1e12, i_impr, 'LineWidth', 2.0, ...
    'DisplayName', ['time step, \delta_{t} = ' num2str(dt * 1e15) ' fs']);
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('t / ps');
ylabel('i(t) / A');
title(['i_{impr} @ V_{b} = ' num2str(Vb) ' V, \tau_{rec} = ' num2str(tau_rec * 1e15) ...
    ' fs, \tau_{s} = ' num2str(tau_s * 1e15) ' fs, \sigma_{t} = ' num2str(sigma_t * 1e15) ' fs']);

%% FUNCTIONS
function i_impr = compute_i_impr(t_vec, k, vb, tau_s, tau_rec, sigma_t)
    m_max = length(t_vec);
    dt = t_vec(2) - t_vec(1);

    [N, M] = meshgrid(1 : 1 : m_max, 1 : 1 : m_max);

    i_impr = tril( exp(- 0.5 * (t_vec(N) / sigma_t) .^ 2) .* exp(- (t_vec(M) - t_vec(N)) / tau_rec) ...
        .* (1 - exp(- (t_vec(M) - t_vec(N)) / tau_s)) );
    i_impr = k * vb * tau_s * dt * sum(i_impr, 2)';
end

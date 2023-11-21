close all;
clear;
clc;

addpath([pwd() '\utils']);

[laser, GaAs] = get_default_config();
[~] = laser.plot_time_envelope();

qo_link = QuasiOpticalLink(laser, GaAs, eta_opt = 0.35);

march_on = TimeStepAlgorithm(laser.get_t_vec(laser.tau_p / 100));
% antenna = PhotoConductiveAntenna(laser, GaAs, march_on, eff_opt = 0.35);

%% PARAMETERS
% LASER
sigma_t = laser.sigma_t;
% PHOTOCONDUCTOR
tau_rec = GaAs.tau_rec;
tau_s = GaAs.tau_s;
% ANTENNA
Vb = 30;
% k = antenna.k_const;
% SIMULATION
dt = laser.tau_p / 1000;
t_vec = laser.get_t_vec(dt);

%% h_{m}
% v = NaN(1, length(t_vec));
% v(1) = 0;
% 
% i_impr = compute_i_impr(t_vec, k, Vb, tau_s, tau_rec, sigma_t);
% 
% figure('Position', [250 250 700 400]);
% plot(t_vec * 1e12, i_impr, 'LineWidth', 2.0, ...
%     'DisplayName', ['time step, \delta_{t} = ' num2str(dt * 1e15) ' fs']);
% grid on;
% legend show;
% legend('location', 'bestoutside');
% xlabel('t / ps');
% ylabel('i(t) / A');
% title(['i_{impr} @ V_{b} = ' num2str(Vb) ' V, \tau_{rec} = ' num2str(tau_rec * 1e15) ...
%     ' fs, \tau_{s} = ' num2str(tau_s * 1e15) ' fs, \sigma_{t} = ' num2str(sigma_t * 1e15) ' fs']);

%% FUNCTIONS
function i_impr = compute_i_impr(t_vec, k, vb, tau_s, tau_rec, sigma_t)
    m_max = length(t_vec);
    dt = t_vec(2) - t_vec(1);

    [N, M] = meshgrid(1 : 1 : m_max, 1 : 1 : m_max);

    i_impr = tril( exp(- 0.5 * (t_vec(N) / sigma_t) .^ 2) .* exp(- (t_vec(M) - t_vec(N)) / tau_rec) ...
        .* (1 - exp(- (t_vec(M) - t_vec(N)) / tau_s)) );
    i_impr = k * vb * tau_s * dt * sum(i_impr, 2)';
end

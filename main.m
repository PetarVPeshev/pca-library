close all;
clear;
clc;

[laser, photoconductor] = get_default_config();
march_on = TimeStepAlgorithm(laser.get_time_vec(laser.pulse_half_pwr / 100));
% laser.plot_time_envelope(laser.pulse_half_pwr / 100);
% eff = Efficiency(laser, photoconductor, eff_opt = 0.35);

antenna = PhotoConductiveAntenna(laser, photoconductor, march_on, eff_opt = 0.35);


%% PARAMETERS
t_vec = laser.get_time_vec(laser.pulse_half_pwr / 100, max_t = 15 * laser.pulse_half_pwr);
k = antenna.k_const;
sigma_t = laser.temporal_std;
tau_rec = photoconductor.rec_time;
tau_s = photoconductor.scat_time;
Vb = 30;

%% h_{m}
% hm = compute_hm(t_vec, k, sigma_t, tau_rec);
% a = compute_Fm_impr(t_vec, hm, tau_s);

v = NaN(1, length(t_vec));
i_impr = NaN(1, length(t_vec));
v(1) = 0;
i_impr(1) = 0;
for m = 2 : 1 : length(v)
    [v, i_impr(m)] = do_time_step(t_vec, v, Vb, k, sigma_t, tau_rec, tau_s);
end

figure('Position', [250 250 700 400]);
plot(t_vec * 1e12, i_impr, 'LineWidth', 2.0, 'DisplayName', 'i_{impr}[n]');
grid on;
% xlim([-2 max(t_vec * 1e12)]);
legend show;
legend('location', 'bestoutside');
xlabel('t / ps');
ylabel('i(t) / A');
title(['i_{impr} @ V_{b} = ' num2str(Vb) ' V, \tau_{rec} = ' num2str(tau_rec * 1e15) ...
    ' fs, \tau_{s} = ' num2str(tau_s * 1e15) ' fs, \sigma_{t} = ' num2str(sigma_t * 1e15) ' fs']);

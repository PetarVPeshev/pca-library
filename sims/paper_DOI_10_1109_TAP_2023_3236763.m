% A. F. Bernardis, P. M. Sberna, J. Bueno, H. Zhang, N. Llombart and A. Neto, "Time-Domain Modelling of 
% Pulsed Photoconducting Sources—Part II: Characterization of an LT GaAs Bow-Tie Antenna," in IEEE 
% Transactions on Antennas and Propagation, vol. 71, no. 3, pp. 2536-2545, March 2023, 
% doi: 10.1109/TAP.2023.3236763.
%
% This simulation replicates figures 3, 4, and 10 from the cited paper.
% 

close all;
clear;
clc;

addpath([pwd() '/..']);

%% PARAMETERS
% LASER
wlen       = 780 * 1e-9;
T          = 12.5 * 1e-9;
P          = [5 30 200] * 1e-3;
R_3db      = 5 * 1e-6;
tau_p      = 100 * 1e-15;
% GaAs
dimensions = [10 10 2] * 1e-6;
er         = 12.96;
tau_rec    = 300 * 1e-15;
tau_s      = 8.5 * 1e-15;
me_coef    = 0.067;
absorp_len = 1 * 1e-6;
eta_opt    = 0.35;
% ANTENNA
Vb         = [30 20 10];
Ra         = 70;
% SIMULATION
t          = (-1 : 6.25 * 1e-4 : 3) * 1e-12;
pwr_pts    = linspace(0.1, 60, 101) * 1e-3 / eta_opt;
eta_qo     = 0.27;

%% ALLOCATE SPACE FOR CURRENT AND VOLTAGES
i          = NaN(length(P), length(t));
v          = NaN(length(P), length(t));

%% PHOTO-CONDUCTIVE ANTENNA
laser = Laser(wlen, T, P(1), R_3db = R_3db, tau_p = tau_p);  % Created at initial laser power
pc    = PhotoConductor(dimensions, er, tau_rec = tau_rec, tau_s = tau_s, me_coef = me_coef, ...
                       absorp_len = absorp_len);
pca   = PhotoConductiveAntenna(laser, pc, Vb(1), 1 / Ra, t_vec = t, eta_opt = eta_opt);

%% ANTENNA RESPONSE SIMULATIONS
[v(1, :), ~, ~, ~, i(1, :)] = pca.compute_response();  % Computed for initial laser power
for pwr_idx = 2 : 1 : length(P)
    laser.P = P(pwr_idx);
    [v(pwr_idx, :), ~, ~, ~, i(pwr_idx, :)] = pca.compute_response();
end

%% PLOT FIGURE 3
figure('Position', [50 50 1200 350]);
sgtitle('DOI: 10.1109/TAP.2023.3236763, FIG 3', 'FontWeight', 'bold');
line_styles = ["-", "--", ":"];

subplot(2, 1, 1);
for pwr_idx = 1 : 1 : length(P)
    plot(t * 1e12, i(pwr_idx, :), 'LineStyle', line_styles(pwr_idx), 'LineWidth', 1.5, ...
        'DisplayName', ['P_{L} = ' num2str(P(pwr_idx) * 1e3) ' mW']);
    hold on;
end
hold off;
grid on;
yticks(0 : 0.1 : 0.4);
legend('Location', 'bestoutside');
ylabel('i(t) / A');

subplot(2, 1, 2);
for pwr_idx = 1 : 1 : length(P)
    plot(t * 1e12, Vb(1) - v(pwr_idx, :), 'LineStyle', line_styles(pwr_idx), 'LineWidth', 1.5, ...
        'DisplayName', ['P_{L} = ' num2str(P(pwr_idx) * 1e3) ' mW']);
    hold on;
end
hold off;
grid on;
yticks(0 : 10 : 30);
legend('Location', 'bestoutside');
ylabel('v_{g}(t) / V');
xlabel('t / ps');

%% POWER SPECTRA
% Number of FFT points and sampling frequency
NFFT = size(i, 2);
fs = (1 / pca.time_step.dt) * NFFT;
% Allocate space for power spectra
Erad = NaN(length(P), floor(NFFT / 2) + 1);
% Perform power spectra computation using Welch method the whole signal without overlap
for pwr_idx = 1 : 1 : length(P)
    [pxx, f] = pwelch(i(pwr_idx, :), NFFT, 0, NFFT, fs);
    Erad(pwr_idx, :) = Ra * pxx';
end
f = transpose(f);

%% PLOT FIGURE 4
figure('Position', [50 50 700 350]);
for pwr_idx = 1 : 1 : length(P)
    plot(f * 1e-12, 10 * log10(Erad(pwr_idx, :)), 'LineStyle', line_styles(pwr_idx), 'LineWidth', 1.5, ...
        'DisplayName', ['P_{L} = ' num2str(P(pwr_idx) * 1e3) ' mW']);
    hold on;
end
grid on;
xlim([0 1e5]);
ylim([-260 -140]);
legend('location', 'bestoutside');
ylabel('E_{rad}(f) / dBJ/Hz');
xlabel('f / THz');
title('DOI: 10.1109/TAP.2023.3236763, FIG 4');

%% ALLOCATE SPACE FOR POWER AND EFFICIENCY
Pdet        = NaN(length(Vb), length(pwr_pts));

%% DETECTED POWER SIMULATIONS
for Vb_idx = 1 : 1 : length(Vb)
    pca.Vb = Vb(Vb_idx);

    for pwr_idx = 1 : 1 : length(pwr_pts)
        laser.P = pwr_pts(pwr_idx);
        [v_pwr, ~, ~, ~, i_pwr] = pca.compute_response();
    
        Pdet(Vb_idx, pwr_idx)    = pca.compute_rad_power(pca.time_step.dt, laser.T, i_pwr, v_pwr) * eta_qo;
    end
end

%% PLOT FIGURE 10
figure('Position', [50 50 700 350]);
for Vb_idx = 1 : 1 : length(Vb)
    plot(pwr_pts * eta_opt * 1e3, 10 * log10(Pdet(Vb_idx, :) * 1e3), 'LineWidth', 1.5, ...
        'DisplayName', ['V_{b} = ' num2str(Vb(Vb_idx)) ' V']);
    hold on;
end
grid on;
% xlim([0 1e5]);
ylim([-30 -10]);
legend('location', 'bestoutside');
ylabel('P_{det} / dBm');
xlabel('P_{opt} / mW');
title('DOI: 10.1109/TAP.2023.3236763, FIG 10');

% A. Neto, N. L. Juan and A. Freni, "Time-Domain Modelling of Pulsed  Photoconducting Sources—Part I: The 
% Norton Equivalent Circuit," in IEEE Transactions on Antennas and Propagation, vol. 71, no. 3, pp. 2527-2535,
% March 2023, doi: 10.1109/TAP.2022.3184517.
% 
% This simulation replicates figures 9, 10, and 11 from the cited paper.
% 

close all;
clear;
clc;

addpath([pwd() '/..']);

%% PARAMETERS
% LASER
wlen       = 780 * 1e-9;
T          = 12.5 * 1e-9;
P          = [5 10 20 40 80] * 1e-3;
R_3db      = 5 * 1e-6;
tau_p      = 0.1 * 1e-12;
% GaAs
dimensions = [10 10 2] * 1e-6;
er         = 12.96;
tau_rec    = 0.3 * 1e-12;
tau_s      = 8 * 1e-15;
me_coef    = 0.067;
absorp_len = 1 * 1e-6;
% ANTENNA
Vb         = 30;
Ra         = [35 70 140 280];
% SIMULATION TIME
t          = linspace(-0.5, 2, 4001) * 1e-12;
pwr_pts    = linspace(0.1, 100, 101) * 1e-3;

%% ALLOCATE SPACE FOR CURRENT AND VOLTAGES
i          = NaN(length(P), length(t));
i_impr     = NaN(length(P), length(t));
i_int      = NaN(length(P), length(t));
v          = NaN(length(P), length(t));

%% PHOTO-CONDUCTIVE ANTENNA
laser = Laser(wlen, T, P(1), R_3db = R_3db, tau_p = tau_p);  % Created at initial laser power
pc    = PhotoConductor(dimensions, er, tau_rec = tau_rec, tau_s = tau_s, me_coef = me_coef, ...
                    absorp_len = absorp_len);
pca   = PhotoConductiveAntenna(laser, pc, Vb, 1 / Ra(2), t_vec = t, eta_opt = 1);

%% ANTENNA RESPONSE SIMULATIONS
[v(1, :), ~, i_impr(1, :), i_int(1, :), i(1, :)] = pca.compute_response();  % Computed for initial laser power
for pwr_idx = 2 : 1 : length(P)
    laser.P = P(pwr_idx);
    [v(pwr_idx, :), ~, i_impr(pwr_idx, :), i_int(pwr_idx, :), i(pwr_idx, :)] = pca.compute_response();
end

%% PLOT FIGURE 9
figure('Position', [50 50 1400 700]);
sgtitle('DOI: 10.1109/TAP.2022.3184517, FIG 9', 'FontWeight', 'bold');

subplot(4, 1, 1);
for pwr_idx = 1 : 1 : length(P)
    plot(t * 1e12, i_impr(pwr_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', ['P_{opt} = ' num2str(P(pwr_idx) * 1e3) ' mW']);
    hold on;
end
hold off;
grid on;
yticks(0 : 0.5 : 3);
legend('Location', 'bestoutside');
ylabel('i_{impr} / A');

subplot(4, 1, 2);
for pwr_idx = 1 : 1 : length(P)
    plot(t * 1e12, i(pwr_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', ['P_{opt} = ' num2str(P(pwr_idx) * 1e3) ' mW']);
    hold on;
end
hold off;
grid on;
yticks(0 : 0.05 : 0.4);
legend('Location', 'bestoutside');
ylabel('i / A');

subplot(4, 1, 3);
for pwr_idx = 1 : 1 : length(P)
    plot(t * 1e12, v(pwr_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', ['P_{opt} = ' num2str(P(pwr_idx) * 1e3) ' mW']);
    hold on;
end
hold off;
grid on;
yticks(0 : 5 : 30);
legend('Location', 'bestoutside');
ylabel('v / V');

subplot(4, 1, 4);
for pwr_idx = 1 : 1 : length(P)
    plot(t * 1e12, Vb - v(pwr_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', ['P_{opt} = ' num2str(P(pwr_idx) * 1e3) ' mW']);
    hold on;
end
hold off;
grid on;
yticks(0 : 5 : 30);
legend('Location', 'bestoutside');
ylabel('V_{b} - v / V');
xlabel('t / ps');

%% ALLOCATE SPACE FOR POWER AND EFFICIENCY
Prad        = NaN(length(Ra), length(pwr_pts));
Ptot        = NaN(length(Ra), length(pwr_pts));
Ppseudo     = NaN(length(Ra), length(pwr_pts));
eta_bias    = NaN(length(Ra), length(pwr_pts));
eta_pca     = NaN(length(Ra), length(pwr_pts));

%% POWER & EFFICIENCY SIMULATIONS
for ra_idx = 1 : 1 : length(Ra)
    pca.ga = 1 / Ra(ra_idx);

    for pwr_idx = 1 : 1 : length(pwr_pts)
        laser.P = pwr_pts(pwr_idx);
        [v_pwr, ~, ~, i_int_pwr, i_pwr] = pca.compute_response();
    
        Prad(ra_idx, pwr_idx)    = pca.compute_rad_power(pca.time_step.dt, laser.T, i_pwr, v_pwr);
        Ptot(ra_idx, pwr_idx)    = pca.compute_total_power(pca.time_step.dt, laser.T, i_pwr, pca.Vb);
        Ppseudo(ra_idx, pwr_idx) = pca.compute_pseudo_power(pca.time_step.dt, laser.T, i_int_pwr, v_pwr);

        [eta_pca(ra_idx, pwr_idx), eta_bias(ra_idx, pwr_idx)] ...
            = pca.compute_efficiencies(Prad(ra_idx, pwr_idx), Ptot(ra_idx, pwr_idx), laser.P);
    end
end

%% PLOT FIGURE 10
figure('Position', [50 50 750 350]);
plot(pwr_pts * 1e3, Prad(2, :) * 1e3, 'LineWidth', 1.5, 'DisplayName', 'P_{l}');
hold on;
plot(pwr_pts * 1e3, Ptot(2, :) * 1e3, '--', 'LineWidth', 1.5, 'DisplayName', 'P_{tot}');
hold on;
plot(pwr_pts * 1e3, Ppseudo(2, :) * 1e3, ':', 'LineWidth', 1.5, 'DisplayName', 'P_{pseudo}');
grid on;
yticks(0 : 0.2 : 1.6);
ylim([0 1.6]);
legend('location', 'bestoutside');
xlabel('P_{opt} / mW');
ylabel('P / mW');
title('DOI: 10.1109/TAP.2022.3184517, FIG 10');

%% PLOT FIGURE 11
figure('Position', [50 50 1350 350]);
sgtitle('DOI: 10.1109/TAP.2022.3184517, FIG 11', 'FontWeight', 'bold');

subplot(1, 2, 1);
for ra_idx = 1 : 1 : length(Ra)
    plot(pwr_pts * 1e3, eta_bias(ra_idx, :) * 100, 'LineWidth', 1.5, ...
        'DisplayName', ['R_{a} = ' num2str(Ra(ra_idx)) ' \Omega']);
    hold on;
end
grid on;
legend('location', 'bestoutside');
xlabel('P_{opt} / mW');
ylabel('\eta_{bias} / %');

subplot(1, 2, 2);
for ra_idx = 1 : 1 : length(Ra)
    plot(pwr_pts * 1e3, eta_pca(ra_idx, :) * 100, 'LineWidth', 1.5, ...
        'DisplayName', ['R_{a} = ' num2str(Ra(ra_idx)) ' \Omega']);
    hold on;
end
grid on;
legend('location', 'bestoutside');
xlabel('P_{opt} / mW');
ylabel('\eta_{pca} / %');

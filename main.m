close all;
clear;
clc;

addpath([pwd() '\utils']);
addpath([pwd() '\slots']);

%% PCA DEFAULT CONFIGURATION
[laser, GaAs, pca] = get_default_config(eta_opt = 0.35);

%% ANTENNA RESPONSE
[v, vg, i_impr, i_int, i] = pca.compute_response();

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

plot(pca.time_step.t * 1e12, exp(- 0.5 * (pca.time_step.t / laser.sigma_t) .^ 2), ...
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

plot(pca.time_step.t * 1e12, i_impr, 'LineWidth', 2.0, ...
    'DisplayName', ['i_{impr}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
hold on;
plot(pca.time_step.t * 1e12, i_int, 'LineWidth', 2.0, ...
    'DisplayName', ['i_{int}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
hold on;
plot(pca.time_step.t * 1e12, i, 'LineWidth', 2.0, ...
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
plot(pca.time_step.t * 1e12, v, 'LineWidth', 2.0, ...
    'DisplayName', ['v, \delta_t = ' num2str(pca.time_step.dt * 1e15) ' fs']);
hold on;
plot(pca.time_step.t * 1e12, vg, 'LineWidth', 2.0, ...
    'DisplayName', ['v_{g}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);

grid on;
xlim([-1 2]);
ylim([0 35]);
legend('location', 'bestoutside');

xlabel('t / ps');
ylabel('v / V');

%% PCA POWER AND EFFICIENCY
[Prad, ~]           = pca.compute_rad_power(pca.time_step.dt, laser.T, i, v);
[Pt, ~]             = pca.compute_total_power(pca.time_step.dt, laser.T, i, pca.Vb);
[Ppseudo, ~]        = pca.compute_pseudo_power(pca.time_step.dt, laser.T, i_int, v);
[eta_pca, eta_bias] = pca.compute_efficiencies(Prad, Pt, pca.P_opt);

fprintf('Prad = %.2f mW\n', Prad * 1e3);
fprintf('Pt = %.2f mW\n', Pt * 1e3);
fprintf('Ppseudo = %.2f mW\n', Ppseudo * 1e3);
fprintf('Popt = %.2f mW\n', pca.P_opt * 1e3);
fprintf('eta_pca = %.2f\n', eta_pca);
fprintf('eta_bias = %.2f\n', eta_bias);

%% INFINITE SLOT PARAMETERS
% DIELECTRICS
er_up = 11.7;
er_dn = 1;
% FREQUENCY
f = 10 * 1e9;
wlen = physconst('LightSpeed') / f;
k0 = 2 * pi / wlen;
% SLOT WIDTH
ws = 0.03 * wlen;
% FEED DELTA-GAP
d_gap = 0.06 * wlen;
% SLOT POINTS
x = linspace(-1, 1, 1001) * wlen;

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% NORMALIZED MAGNETIC CURRENT / VOLTAGE
v = slot.compute_v(x, f);

figure('Name', 'Normalized Magnetic Current', 'Position', [150 150 800 350]);

plot(x / wlen, real(v), 'LineWidth', 1.5, 'DisplayName', '\Re\{v\}');
hold on;
plot(x / wlen, imag(v), 'LineWidth', 1.5, 'DisplayName', '\Im\{v\}');

grid on;
legend('location', 'bestoutside');

xlabel('x / \lambda_{0}');
ylabel('v / V');
title(['Normalized Magnetic Current @ \Delta = 0.06\lambda_{0}, w_{s} = 0.03\lambda_{0}, ' ...
    '\epsilon_{r,2} = 11.7, \epsilon_{r,1} = 1']);

%% PROPAGATION CONSTANT
% DIELECTRIC ARRAY
er_up = [11.7 4 2.54];
% SLOT WIDTH ARRAY
ws = linspace(1e-5, 0.1, 101) * wlen;
% POLES
kxp = NaN(3, 101);
kxp_analyt = NaN(3, 101);
for er_idx = 1 : 1 : length(er_up)
    for ws_idx = 1 : 1 : length(ws)
        slot = SlotInDielectrics(d_gap, ws(ws_idx), er_up(er_idx), er_dn);
    
        kxp(er_idx, ws_idx) = slot.find_kxp(f);
        kxp_analyt(er_idx, ws_idx) = slot.evaluate_kxp_analyt(f);
    end
end

figure('Name', 'Leaky-Wave Propagation Constant', 'Position', [150 150 800 500]);

subplot(2, 1, 1);

plot(ws / wlen, real(kxp(1, :) / k0), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\epsilon_{r,2} = 11.7, numerical');
hold on;
plot(ws / wlen, real(kxp_analyt(1, :) / k0), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\epsilon_{r,2} = 11.7, analytical');
hold on;
plot(ws / wlen, real(kxp(2, :) / k0), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\epsilon_{r,2} = 4, numerical');
hold on;
plot(ws / wlen, real(kxp_analyt(2, :) / k0), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\epsilon_{r,2} = 4, analytical');
hold on;
plot(ws / wlen, real(kxp(3, :) / k0), 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\epsilon_{r,2} = 2.54, numerical');
hold on;
plot(ws / wlen, real(kxp_analyt(3, :) / k0), '--', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\epsilon_{r,2} = 2.54, analytical');

grid on;
legend('location', 'bestoutside');

ylabel('\Re\{k_{x}^{LW} / k_{0}\}');
ylim([1 3.5]);
xticks(0 : 0.02 : 0.1);

subplot(2, 1, 2);

plot(ws / wlen, imag(kxp(1, :) / k0), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\epsilon_{r,2} = 11.7, numerical');
hold on;
plot(ws / wlen, imag(kxp_analyt(1, :) / k0), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\epsilon_{r,2} = 11.7, analytical');
hold on;
plot(ws / wlen, imag(kxp(2, :) / k0), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\epsilon_{r,2} = 4, numerical');
hold on;
plot(ws / wlen, imag(kxp_analyt(2, :) / k0), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\epsilon_{r,2} = 4, analytical');
hold on;
plot(ws / wlen, imag(kxp(3, :) / k0), 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\epsilon_{r,2} = 2.54, numerical');
hold on;
plot(ws / wlen, imag(kxp_analyt(3, :) / k0), '--', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\epsilon_{r,2} = 2.54, analytical');

grid on;
legend('location', 'bestoutside');

ylabel('\Im\{k_{x}^{LW} / k_{0}\}');
xlabel('w_{s} / \lambda_{0}');
ylim([-0.5 0]);
xticks(0 : 0.02 : 0.1);

sgtitle('Leaky-Wave Propagation Constant', 'FontSize', 11, 'FontWeight', 'bold');

%% LEAKY-WAVE MAGNETIC CURRENT / VOLTAGE
% DIELECTRIC
er_up = 11.7;
% SLOT WIDTH
ws = 0.03 * wlen;
% SLOT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);
% LEAKY-WAVE MAGNETIC CURRENT
vlw = slot.compute_vlw(x, f);

figure('Name', 'Normalized Leaky-Wave Magnetic Current', 'Position', [150 150 800 350]);

plot(x / wlen, real(vlw), 'LineWidth', 1.5, 'DisplayName', '\Re\{v_{LW}\}');
hold on;
plot(x / wlen, imag(vlw), 'LineWidth', 1.5, 'DisplayName', '\Im\{v_{LW}\}');

grid on;
legend('location', 'bestoutside');

xlabel('x / \lambda_{0}');
ylabel('v_{LW} / V');
title(['Normalized Leaky-Wave Magnetic Current @ \Delta = 0.06\lambda_{0}, w_{s} = 0.03\lambda_{0}, ' ...
    '\epsilon_{r,2} = 11.7, \epsilon_{r,1} = 1']);

%% INPUT IMPEDANCE
% SLOT WIDTH
ws = [1/30 1/100 1/1000] * wlen;
% FEED DELTA-GAP
d_gap = linspace(0, 3, 1001) * wlen / 100;
% INPUT IMPEDANCE
Zin = NaN(length(ws), length(d_gap));
for ws_idx = 1 : 1 : length(ws)
    for d_gap_idx = 1 : 1 : length(d_gap)
        slot = SlotInDielectrics(d_gap(d_gap_idx), ws(ws_idx), er_up, er_dn);

        Zin(ws_idx, d_gap_idx) = slot.compute_zin(f);
    end
end

figure('Name', 'Infinite Slot Impedance', 'Position', [150 150 800 500]);

plot(100 * d_gap / wlen, real(Zin(1, :)), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{Z_{in}\}, w_{s}=\lambda_{0}/30');
hold on;
plot(100 * d_gap / wlen, imag(Zin(1, :)), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{Z_{in}\}, w_{s}=\lambda_{0}/30');
hold on;

plot(100 * d_gap / wlen, real(Zin(2, :)), '--', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Re\{Z_{in}\}, w_{s}=\lambda_{0}/100');
hold on;
plot(100 * d_gap / wlen, imag(Zin(2, :)), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', '\Im\{Z_{in}\}, w_{s}=\lambda_{0}/100');
hold on;

plot(100 * d_gap / wlen, real(Zin(3, :)), '--', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Re\{Z_{in}\}, w_{s}=\lambda_{0}/1000');
hold on;
plot(100 * d_gap / wlen, imag(Zin(3, :)), 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Im\{Z_{in}\}, w_{s}=\lambda_{0}/1000');

grid on;
ylim([-20 100]);
legend('location', 'bestoutside');

xlabel('100 \Delta_{gap} / \lambda_{0}');
ylabel('Z / \Omega');
title('Infinite Slot Impedance @ \epsilon_{r,1} = 1, \epsilon_{r,2} = 11.7');

%% QUASI-STATIC AND DYNAMIC IMPEDANCE
% SLOT WIDTH
ws = linspace(1e-6, 0.01, 101) * wlen;
% FEED DELTA-GAP
d_gap = 0.001 * wlen;
% IMPEDANCE
Zdyn = NaN(1, length(ws));
Zqs = NaN(1, length(ws));
Zin = NaN(1, length(ws));
for ws_idx = 1 : 1 : length(ws)
    slot = SlotInDielectrics(d_gap, ws(ws_idx), er_up, er_dn);

    Zdyn(ws_idx) = slot.compute_Zdyn(f);
    Zqs(ws_idx) = slot.compute_Zqs(f);
    Zin(ws_idx) = slot.compute_zin(f);
end
Zin_sum = Zdyn + Zqs;

figure('Name', 'Infinite Slot Impedance Components', 'Position', [150 150 800 500]);

plot(ws / wlen, real(Zin), '-+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(ws), 'MarkerSize', 5, 'DisplayName', '\Re\{Z_{in}\}, num');
hold on;
plot(ws / wlen, imag(Zin), '--+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(ws), 'MarkerSize', 5, 'DisplayName', '\Im\{Z_{in}\}, num');

hold on;
plot(ws / wlen, real(Zin_sum), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{Z_{in}\}');
hold on;
plot(ws / wlen, imag(Zin_sum), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{Z_{in}\}');
hold on;
plot(ws / wlen, real(Zqs), '-.', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Im\{Z_{qs}\}');
hold on;
plot(ws / wlen, imag(Zdyn), ':', 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560], ...
    'DisplayName', '\Im\{Z_{dyn}\}');

grid on;
ylim([0 40]);
legend('location', 'bestoutside');

xlabel('w_{s} / \lambda_{0}');
ylabel('Z / \Omega');
title(['Infinite Slot Impedance Components @ \epsilon_{r,1} = 1, \epsilon_{r,2} = 11.7, ' ...
    '\Delta = 0.001\lambda_{0}']);

% SLOT WIDTH
ws = 0.01 * wlen;
% FEED DELTA_GAP
d_gap = linspace(1e-6, 0.1, 101) * wlen;
% IMPEDANCE
Zdyn = NaN(1, length(d_gap));
Zqs = NaN(1, length(d_gap));
Zin = NaN(1, length(d_gap));
for d_gap_idx = 1 : 1 : length(d_gap)
    slot = SlotInDielectrics(d_gap(d_gap_idx), ws, er_up, er_dn);

    Zdyn(d_gap_idx) = slot.compute_Zdyn(f);
    Zqs(d_gap_idx) = slot.compute_Zqs(f);
    Zin(d_gap_idx) = slot.compute_zin(f);
end
Zin_sum = Zdyn + Zqs;

figure('Name', 'Infinite Slot Impedance Components', 'Position', [150 150 800 500]);

plot(d_gap / wlen, real(Zin), '-+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(d_gap), 'MarkerSize', 5, 'DisplayName', '\Re\{Z_{in}\}, num');
hold on;
plot(d_gap / wlen, imag(Zin), '--+', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880], ...
    'MarkerIndices', 1 : 5 : length(d_gap), 'MarkerSize', 5, 'DisplayName', '\Im\{Z_{in}\}, num');
hold on;

plot(d_gap / wlen, real(Zin_sum), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Re\{Z_{in}\}');
hold on;
plot(d_gap / wlen, imag(Zin_sum), '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], ...
    'DisplayName', '\Im\{Z_{in}\}');
hold on;
plot(d_gap / wlen, imag(Zqs), '-.', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250], ...
    'DisplayName', '\Im\{Z_{qs}\}');
hold on;
plot(d_gap / wlen, imag(Zdyn), ':', 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560], ...
    'DisplayName', '\Im\{Z_{dyn}\}');

grid on;
ylim([-10 50]);
legend('location', 'bestoutside');

xlabel('\Delta_{gap} / \lambda_{0}');
ylabel('Z / \Omega');
title(['Infinite Slot Impedance Components @ \epsilon_{r,1} = 1, \epsilon_{r,2} = 11.7, ' ...
    'w_{s} = 0.01\lambda_{0}']);

%% SEMI-FINITE SLOT PARAMETERS
% DIELECTRICS
er_up = 1;
er_dn = 1;
% FREQUENCY
f0 = 10 * 1e9;
f = linspace(0.5, 2, 101) * f0;
wlen = physconst('LightSpeed') / f0;
% SLOT WIDTH
ws = wlen / 50;
% FEED DELTA-GAP
d_gap = wlen / 40;
% TERMINATION DISTANCE
ds = wlen / 4;
% SLOT POINTS
x = linspace(-1, 1, 101) * wlen;

%% SLOT OBJECT
slot = SlotSemiFinite(d_gap, ws, ds, er_up, er_dn);

%% INPUT IMPEDANCE
Zin = NaN(1, length(f));
for f_idx = 1 : 1 : length(f)
    Zin(f_idx) = slot.compute_zin(f(f_idx));
end

figure('Name', 'Semi-Finite Slot Impedance', 'Position', [150 150 700 400]);

plot(f / f0, real(Zin), 'LineWidth', 1.5, 'DisplayName', '\Re\{Z_{in}\}');
hold on;
plot(f / f0, imag(Zin), 'LineWidth', 1.5, 'DisplayName', '\Im\{Z_{in}\}');

grid on;
ylim([-50 150]);
legend('Location', 'bestoutside');

xlabel('f / f_{0}');
ylabel('Z_{in} / \Omega');
title(['Semi-Infinite Slot Impedance @ free space, d = \lambda_{0} / 4, w_{s} = \lambda_{0} / 50, ' ...
    '\Delta = \lambda_{0} / 40']);

%% NORMALIZED MAGNETIC CURRENT
v = slot.compute_v(x, f0);

figure('Name', 'Semi-Finite Slot Normalized Magnetic Current', 'Position', [150 150 700 400]);

plot(x / wlen, real(v), 'LineWidth', 1.5, 'DisplayName', '\Re\{v\}');
hold on;
plot(x / wlen, imag(v), 'LineWidth', 1.5, 'DisplayName', '\Im\{v\}');

grid on;
legend('Location', 'bestoutside');

xlabel('x / d');
ylabel('v / V');
title(['Semi-Infinite Slot Normalized Magnetic Current @ free space, d = \lambda_{0} / 4, ' ...
    'w_{s} = \lambda_{0} / 50, \Delta = \lambda_{0} / 40']);

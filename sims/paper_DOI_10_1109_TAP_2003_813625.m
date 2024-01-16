% A. Neto and S. Maci, "Green's function for an infinite slot printed between two homogeneous dielectrics. 
% I. Magnetic currents," in IEEE Transactions on Antennas and Propagation, vol. 51, no. 7, pp. 1572-1581, 
% July 2003, doi: 10.1109/TAP.2003.813625.
% 
% This simulation replicates figures 2, 3, and 6 from the cited paper.
% 

close all;
clear;
clc;

addpath([pwd() '/..']);
addpath([pwd() '/../slots']);
addpath([pwd() '/../utils']);

%% ARAMETERS
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

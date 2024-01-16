% R. M. van Schelven, D. Cavallo and A. Neto, "Equivalent Circuit Models of Finite Slot Antennas," in 
% IEEE Transactions on Antennas and Propagation, vol. 67, no. 7, pp. 4367-4376, July 2019, 
% doi: 10.1109/TAP.2019.2908112.
% 
% This simulation replicates figures 5a from the cited paper.
% 

close all;
clear;
clc;

addpath([pwd() '/..']);
addpath([pwd() '/../slots']);
addpath([pwd() '/../utils']);

%% PARAMETERS
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

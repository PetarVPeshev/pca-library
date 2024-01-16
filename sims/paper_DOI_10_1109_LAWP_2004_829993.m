% Neto and Maci, "Input impedance of slots printed between two dielectric media and fed by a 
% small /spl Delta/-gap," in IEEE Antennas and Wireless Propagation Letters, vol. 3, pp. 113-116, 2004, 
% doi: 10.1109/LAWP.2004.829993.
% 
% This simulation replicates figures 2 and 3 from the cited paper.
% 

close all;
clear;
clc;

addpath([pwd() '/..']);
addpath([pwd() '/../slots']);
addpath([pwd() '/../utils']);

%% PARAMETERS
% DIELECTRICS
er_up = 11.7;
er_dn = 1;
% FREQUENCY
f = 10 * 1e9;
wlen = physconst('LightSpeed') / f;
k0 = 2 * pi / wlen;
% SLOT WIDTH
ws = [1/30 1/100 1/1000] * wlen;
% FEED DELTA-GAP
d_gap = linspace(0, 3, 1001) * wlen / 100;
% SLOT POINTS
x = linspace(-1, 1, 1001) * wlen;

%% INPUT IMPEDANCE
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

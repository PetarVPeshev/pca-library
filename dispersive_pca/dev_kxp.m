close all;
clear;
clc;

addpath('..\slots');
addpath('..\utils');

c0 = get_phys_const('LightSpeed');

LineLegend = {'\Re', '\Im'};
Location = 'bestoutside';
LineWidth = 1.5;
Color = ["#0072BD", "#EDB120", "#77AC30", "#A2142F"];
Position = [680 558 700 420];

%% PARAMETERS
f = linspace(eps, 2e12, 1001);
% FEED GAP
d_gap = 4.5 * 1e-6;
d_gap_vec = (4.5 : - 1 : 1.5) * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
ws_vec = (10 : - 2.5 : 2.5) * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1

%% POLE / LEAKY-WAVE PROPAGATION CONSTANT
kxp = NaN(1, length(f));
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

for f_idx = 1 : 1 : length(f)
    kxp(f_idx) = slot.find_kxp(f(f_idx));
end

fig = figure('Position', Position);

SPLT1 = subplot(2, 1, 1);
plot(f * 1e-12, real(kxp * c0 ./ (2 * pi * f)), 'LineWidth', 1.5);

set(SPLT1, 'Box', 'off');
ylim([2.55 2.9]);
yticks(2.55 : 0.05 : 2.9);
grid on;
ylabel('\Re(k_{xp}/k_{0})');

SPLT2 = subplot(2, 1, 2);
plot(f * 1e-12, imag(kxp * c0 ./ (2 * pi * f)), 'LineWidth', 1.5);

set(SPLT2, 'Box', 'off');
ylim([-0.5 -0.15]);
yticks(-0.5 : 0.05 : -0.15);
grid on;
ylabel('\Im(k_{xp}/k_{0})');

xlabel('f [THz]');
sgtitle(['@ w_{y} = ' num2str(ws * 1e6) ' \mum'], 'FontSize', 11, 'FontWeight', 'bold');

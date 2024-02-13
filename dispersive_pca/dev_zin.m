close all;
clear;
clc;

addpath('..\slots', '..\utils');

LineLegend = {'\Re', '\Im'};
Location   = 'bestoutside';
LineWidth  = 1.5;
Color      = ["#0072BD", "#EDB120", "#77AC30", "#A2142F"];
Position   = [680 558 700 420];

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

%% OBJECTS
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% INPUT IMPEDANCE OF INFINITE SLOT
% FOR DIFFERENT SLOT WIDTHS
Zin         = NaN(length(ws_vec), length(f));
ColorLegend = cell(1, 4);

for ws_idx = 1 : length(ws_vec)
    slot.ws        = ws_vec(ws_idx);
    Zin(ws_idx, :) = slot.compute_zin(f);

    ColorLegend(ws_idx) = cellstr( ['w_{y} = ' num2str(ws_vec(ws_idx) * 1e6) ' \mum'] );
end

XLabel = 'f [THz]';
YLabel = 'z_{in} [\Omega]';
Title = ['@ \Delta = ' num2str(d_gap * 1e6) ' \mum'];
YLim = [0 80];

plot_two(f * 1e-12, real(Zin), imag(Zin), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'YLim', YLim, 'Position', Position);
grid on;
yticks(0 : 20 : 80);

% FOR DIFFERENT FEED DELTA-GAPS
slot.ws = ws;

Zin         = NaN(length(d_gap_vec), length(f));
ColorLegend = cell(1, 4);

for d_gap_idx = 1 : length(d_gap_vec)
    slot.d_gap        = d_gap_vec(d_gap_idx);
    Zin(d_gap_idx, :) = slot.compute_zin(f);

    ColorLegend(d_gap_idx) = cellstr( ['\Delta = ' num2str(d_gap_vec(d_gap_idx) * 1e6) ' \mum'] );
end

XLabel = 'f [THz]';
YLabel = 'z_{in} [\Omega]';
Title = ['@ w_{s} = ' num2str(ws * 1e6) ' \mum'];
YLim = [0 120];

plot_two(f * 1e-12, real(Zin), imag(Zin), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
    'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
    Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'YLim', YLim, 'Position', Position);
grid on;
yticks(0 : 20 : 120);

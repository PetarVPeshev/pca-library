close all;
clear;
clc;

env;
sim_params;

ws = 10 * 1e-6;
d_gap = 10 * 1e-6;

slot  = create_slot(params_slot, 'd_gap', 10 * 1e-6, 'ws', ws);
pcm   = create_pcm(params_pcm, 'd_gap', 10 * 1e-6, 'ws', ws);
laser = create_laser(params_laser, 'P', 10 * 1e-3);

[t, tt] = create_time_arrays(1 * 1e-15, [-0.5 5] * 1e-12);

x_metal = (- 200 : 1 : 200) * 1e-6;
d = 100 * 1e-6;

x_metal = [100 100 2 2] * 1e-6;
y_metal = [2 * ws, ws / 2, ws / 2, 2 * ws];
x_feed = - d / 2 - [-d_gap / 2, -d_gap / 2, d_gap / 2, d_gap / 2];
y_feed = [ws / 2, - ws / 2, - ws / 2, ws / 2];
metal_1 = polyshape(x_metal, y_metal);
metal_2 = polyshape(-x_metal, y_metal);
metal_3 = polyshape(x_metal, -y_metal);
metal_4 = polyshape(-x_metal, -y_metal);
feed_1 = polyshape(- x_feed, y_feed);
feed_2 = polyshape(x_feed, y_feed);
COLOR = '#404040';

figure('Position', [250 250 950 450]);
% y = [25 25];  % ; -ws / 2, -25; -25 - ws / 2
% x = [-2 -0.2];
plot(metal_1, 'FaceColor', COLOR, 'LineStyle', 'none');
hold on;
plot(metal_2, 'FaceColor', COLOR, 'LineStyle', 'none');
hold on;
plot(metal_3, 'FaceColor', COLOR, 'LineStyle', 'none');
hold on;
plot(metal_4, 'FaceColor', COLOR, 'LineStyle', 'none');
hold on;
plot(feed_1, 'FaceColor', '#B73905', 'LineStyle', 'none');
hold on;
plot(feed_2, 'FaceColor', '#B73905', 'LineStyle', 'none');

xlim([-100 100] * 1e-6);
ylim([-50 50] * 1e-6);

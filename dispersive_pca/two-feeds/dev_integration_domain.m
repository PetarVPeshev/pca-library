close all;
clear;
clc;

env;
sim_params;

c = physconst('LightSpeed');

%% PARAMETERS
f     = 2 * 1e12;
d_gap = 4.5 * 1e-6;
kx    = linspace(-800, 800, 1001);
% DEPENDENT
wlen = c / f;
k0   = 2 * pi / wlen;

%% SINC
s  = sinc(kx * k0 * d_gap / (2 * pi));

a  = wlen / d_gap;
af = get_integration_domain(d_gap, f);

% PLOT
figure('Position', [250 250 750 450]);

SP  = plot(kx, s, 'LineWidth', 1.5, 'DisplayName', 'sin(x) / x');
hold on;
AP  = xline(a, ':', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', '1st zero');
hold on;
AFP = xline(af, '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', '15th zero');
hold on;
xline(- a, ':', 'LineWidth', 1.5, 'Color', "#7E2F8E", 'DisplayName', '1st zero');
hold on;
xline(- af, '--', 'LineWidth', 1.5, 'Color', "#A2142F", 'DisplayName', '15th zero');

grid on;
legend([SP AP AFP], 'location', 'bestoutside');
annotation('doublearrow', [0.23 0.632], [0.8 0.8], 'LineWidth', 1.5, 'Color', '#D95319');
annotation('textbox', [0.235 0.76 0.22 0.1], 'String', 'integration domain', 'Color', '#D95319', ...
           'EdgeColor', 'none', 'FontWeight', 'bold');

xlabel('\Re \{k_{x} / k_{0}\}');
ylabel('sinc(k_{x}\Delta / 2)');
title(['@ \Delta = ' num2str(d_gap * 1e6) ' \mum, \lambda_{0} = ' num2str(round(wlen * 1e6, 2)) ' \mum']);

%% FUNCTIONS
function a = get_integration_domain(d_gap, f, N)
    arguments
        d_gap (1,1) double {mustBePositive}
        f     (1,1) double {mustBePositive}
        N     (1,1) double {mustBePositive} = 15
    end

    c = physconst('LightSpeed');
    wlen = c / f;

    a = N * wlen / d_gap;
end

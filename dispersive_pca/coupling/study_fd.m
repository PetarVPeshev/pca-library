close all;
clear;
clc;

env;
sim_params;

%% PARAMETERS
% SLOT
ws     = 10 * 1e-6;   % SLOT WIDTH    : 10 um
d_gap  = 10 * 1e-6;   % FEED WIDTH    : 10 um
d_feed = 100 * 1e-6;  % FEED DISTANCE : 100 um
% PHOTO-CONDUCTORS
Vb    = [30 0];    % BIAS VOLTAGE  : 30 V & 0 V
tau_d = [0 0];  % EXCITATION DELAY : 0 ps & 0 ps
% LASER
P = 10 * 1e-3;  % AVERAGE POWER : 10 mW
% SIMULATION TIME
t_lims = [-1 10] * 1e-12;   % TIME LIMITS : -1 ps to 10 ps
dt     = 0.5 * 1e-15;       % TIME STEP   : 0.5 fs
% FREQUENCY DOMAIN
Nf     = 4001;          % FREQUENCY POINTS : 4001
f_lims = [1e9 3e12];    % FREQUENCY LIMITS : 1 MHz to 3 THz

%% FREQUENCY & TIME VECTORS
mustBeIncreasing(f_lims);
f            = linspace(f_lims(1), f_lims(2), Nf);
[tsim, tres] = create_time_arrays(dt, t_lims);
Nt           = length(tsim);

%% OBJECTS
slot  = create_slot(params_slot, 'd_gap', d_gap, 'ws', ws);
pcm   = create_pcm(params_pcm, 'd_gap', d_gap, 'ws', ws);
laser = create_laser(params_laser, 'P', P);

%% IMPEDANCE
Z = eval_Z(f, slot, d_feed);

fig = figure('Position', [250 250 950 450]);

AX1  = axes(fig);
PLTS = gobjects(1, 2);
PLTS(1) = plot(f * 1e-12, real(permute(Z(1, 1, :), [1 3 2])), 'LineWidth', 1.5, 'Color', '#0072BD');
hold on;
          plot(f * 1e-12, imag(permute(Z(1, 1, :), [1 3 2])), '--', 'LineWidth', 1.5, 'Color', '#0072BD');
hold on;
PLTS(2) = plot(f * 1e-12, real(permute(Z(1, 2, :), [1 3 2])), 'LineWidth', 1.5, 'Color', '#A2142F');
hold on;
          plot(f * 1e-12, imag(permute(Z(1, 2, :), [1 3 2])), '--', 'LineWidth', 1.5, 'Color', '#A2142F');

grid on;
LGND1 = legend(PLTS, {'Z_{11} = Z_{22}', 'Z_{12} = Z_{21}'}, 'location', 'bestoutside');
hold on;

AX2 = copyobj(AX1, gcf);
delete( get(AX2, 'Children') );
H1 = plot(NaN, 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
H2 = plot(NaN, '--', 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);

set(AX2, 'Color', 'none', 'XTick', [], 'Visible', 'off');
LGND2 = legend([H1 H2], {'\Re', '\Im'}, 'location', 'bestoutside');
LGND2.Position(1) = LGND1.Position(1);
LGND2.Position(2) = LGND1.Position(2) - LGND2.Position(4) - 0.02;
LGND2.Box = 'off';
LGND2.Color = 'none';
hold off;
    
xlabel('f [THz]');
ylabel('Z [\Omega]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
       num2str(d_feed * 1e6) ' \mum']);

%% ADMITTANCE
H = compute_H(Z);

fig = figure('Position', [250 250 950 450]);

AX1  = axes(fig);
PLTS = gobjects(1, 2);
PLTS(1) = plot(f * 1e-12, real(permute(H(1, 1, :), [1 3 2])), 'LineWidth', 1.5, 'Color', '#0072BD');
hold on;
          plot(f * 1e-12, imag(permute(H(1, 1, :), [1 3 2])), '--', 'LineWidth', 1.5, 'Color', '#0072BD');
hold on;
PLTS(2) = plot(f * 1e-12, real(permute(H(1, 2, :), [1 3 2])), 'LineWidth', 1.5, 'Color', '#A2142F');
hold on;
          plot(f * 1e-12, imag(permute(H(1, 2, :), [1 3 2])), '--', 'LineWidth', 1.5, 'Color', '#A2142F');

grid on;
LGND1 = legend(PLTS, {'Y_{11} = Y_{22}', 'Z_{12} / Z_{11} = Z_{21} / Z_{22}'}, 'location', 'bestoutside');
hold on;

AX2 = copyobj(AX1, gcf);
delete( get(AX2, 'Children') );
H1 = plot(NaN, 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
H2 = plot(NaN, '--', 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);

set(AX2, 'Color', 'none', 'XTick', [], 'Visible', 'off');
LGND2 = legend([H1 H2], {'\Re', '\Im'}, 'location', 'bestoutside');
LGND2.Position(1) = LGND1.Position(1);
LGND2.Position(2) = LGND1.Position(2) - LGND2.Position(4) - 0.02;
LGND2.Box = 'off';
LGND2.Color = 'none';
hold off;
    
xlabel('f [THz]');
ylabel('Z [\Omega]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
       num2str(d_feed * 1e6) ' \mum']);

%% IMPULSE RESPONSE
h = compute_h(tres, f, H);

figure('Position', [250 250 950 450]);

subplot(2, 1, 1);
plot(tres * 1e12, permute(h(1, 1, :), [1 3 2]), 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', 'y_{11} = y_{22}');

grid on;
xlim([0 4]);
legend('location', 'bestoutside');

ylabel('y [S/s]');

subplot(2, 1, 2);
plot(tres * 1e12, permute(h(1, 2, :), [1 3 2]), 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', 'z_{12} / z_{11} = z_{21} / z_{22}');

grid on;
xlim([0 4]);
legend('location', 'bestoutside');

ylabel('h [-]');
xlabel('t [ps]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
        num2str(d_feed * 1e6) ' \mum'], 'FontSize', 11, 'FontWeight', 'bold');

%% TIME-STEP ALGORITHM
time_step = TimeStepCoupling(tsim);

% Set time-step algorithm parameters
time_step.K       = calculate_K(laser, pcm);
time_step.Vb      = Vb;
time_step.h       = h;
time_step.tau_c   = pcm.tau_rec;
time_step.tau_s   = pcm.tau_s;
time_step.tau_d   = tau_d;
time_step.sigma_t = laser.sigma_t;

% Voltage and currents
[v, i] = time_step.simulate();

figure('Position', [250 250 950 450]);

subplot(2, 1, 1);
plot(tsim * 1e12, v(1, :), 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'v_{1}');
hold on;
plot(tsim * 1e12, v(2, :), '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'v_{2}');

grid on;
xlim([-1 4]);
ylim([-0.1 2]);
legend('location', 'bestoutside');

ylabel('v [V]');

subplot(2, 1, 2);
plot(tsim * 1e12, i(1, :), 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'i_{1}');
hold on;
plot(tsim * 1e12, i(2, :), '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'i_{2}');

grid on;
xlim([-1 4]);
ylim([0 0.025]);
legend('location', 'bestoutside');

ylabel('i [A]');
xlabel('t [ps]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
        num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], 'FontSize', 11, ...
        'FontWeight', 'bold');

%% FUNCTIONS
function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

function h = compute_h(t, f, H)
    Nt = length(t);
    h = NaN(2, 2, Nt);

    h(1, 1, :) = 2 * real(eval_IFT(t, f, permute(H(1, 1, :), [1 3 2])));
    h(2, 2, :) = h(1, 1, :);
    h(1, 2, :) = 2 * real(eval_IFT(t, f, permute(H(1, 2, :), [1 3 2])));
    h(2, 1, :) = h(1, 1, :);
end

function h_filt = filter_h(h, dt)
    N  = size(h, 3);
%     fs = 1 / dt;
%     Fstop = 500 * 1e6;
%     Fpass = 750 * 1e6;
%     Apass = 0.01;
%     Astop = 80;
    window_size = 500;

    b = (1 / window_size) * ones(1, window_size);

%     FIReq = fdesign.highpass(Fstop, Fpass, Astop, Apass, fs);
%     filterCoeff = design(FIReq, 'equiripple', 'SystemObject', true);

    h_filt = NaN(2, 2, N);
    h_filt(1, 1, :) = filter(b, 1, permute(h(1, 1, :), [1 3 2]));
    h_filt(1, 2, :) = filter(b, 1, permute(h(1, 2, :), [1 3 2]));
%     h_filt(1, 1, :) = filterCoeff(permute(h(1, 1, :), [1 3 2]));
%     h_filt(1, 2, :) = filterCoeff(permute(h(1, 2, :), [1 3 2]));
    h_filt(2, 1, :) = h_filt(1, 2, :);
    h_filt(2, 2, :) = h_filt(1, 1, :);
end

function h_ifft = compute_h_ifft(H)
    N = size(H, 3);

%     H11 = [permute(H(1, 1, :), [1 3 2]) conj(fliplr(permute(H(1, 1, 2 : end), [1 3 2])))];
%     H12 = [permute(H(1, 2, :), [1 3 2]) conj(fliplr(permute(H(1, 2, 2 : end), [1 3 2])))];
    H11 = permute(H(1, 1, :), [1 3 2]);
    H12 = permute(H(1, 2, :), [1 3 2]);

    h_ifft = NaN(2, 2, N);
    h_ifft(1, 1, :) = 2 * real(ifft(H11));
    h_ifft(1, 2, :) = 2 * real(ifft(H12));
    h_ifft(2, 1, :) = h_ifft(1, 2, :);
    h_ifft(2, 2, :) = h_ifft(1, 1, :);
end

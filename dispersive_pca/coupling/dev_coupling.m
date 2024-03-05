close all;
clear;
clc;

env;
sim_params;

%% PARAMETERS
% SLOT
ws     = 10 * 1e-6;          % SLOT WIDTH      : 10 um
d_gap  = 10 * 1e-6;          % FEED WIDTH      : 10 um
d_feed = 100 * 1e-6;         % FEED DISTANCE   : 100 um (distance b/n feeds)
x_lims = [-100 100] * 1e-6;  % DISTANCE LIMITS : -200 um to 200 um
dx     = 10 * 1e-6;          % DISTANCE STEP   : 10 um
% PHOTO-CONDUCTORS
Vb    = [30; 0];            % BIAS VOLTAGE     : 30 V & 0 V
tau_d = [0 1] * 1e-12;      % EXCITATION DELAY : 0 ps & 0 ps
% LASER
P = 10 * 1e-3;  % AVERAGE POWER : 10 mW
% SIMULATION TIME
t_lims = [-1 10] * 1e-12;         % TIME LIMITS : -1 ps to 10 ps
dt     = 0.5 * 1e-15;             % TIME STEP   : 0.5 fs
t_plot = [5e-2 1 1.5 2] * 1e-12;  % PLOT TIME   : 5 fs, 1 ps, and 1.5 ps
% FREQUENCY DOMAIN
Nf     = 4001;          % FREQUENCY POINTS : 4001
f_lims = [1e9 3e12];    % FREQUENCY LIMITS : 1 MHz to 3 THz

%% FREQUENCY, TIME, & SLOT POSITION VECTORS
mustBeIncreasing(f_lims);
mustBeIncreasing(x_lims);
f            = linspace(f_lims(1), f_lims(2), Nf);
[tsim, tres] = create_time_arrays(dt, t_lims);
Nt           = length(tsim);
x            = x_lims(1) : dx : x_lims(2);
Nx           = length(x);

%% TIME INDEX PLOT
Nt_plot = length(t_plot);
t_idx   = NaN(1, Nt_plot);
for idx = 1 : Nt_plot
    t_idx(idx) = find(round(tsim * 1e12, 4) == t_plot(idx) * 1e12, 1);
end

%% FEED INDEX
x_idx = NaN(1, 2);
for idx = 1 : 2
    x_idx(idx) = find(round(x * 1e6) == round(((-1) ^ idx) * d_feed * 1e6 / 2), 1);
end

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
[H, W] = compute_H(Z);

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
LGND1 = legend(PLTS, {'H_{11} = H_{22}', 'H_{12} = H_{21}'}, 'location', 'bestoutside');
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

fig = figure('Position', [250 250 950 450]);

AX1 = axes(fig);
PLT = plot(f * 1e-12, real(W), 'LineWidth', 1.5, 'Color', '#77AC30');
hold on;
      plot(f * 1e-12, imag(W), '--', 'LineWidth', 1.5, 'Color', '#77AC30');

grid on;
LGND1 = legend(PLT, 'W', 'location', 'bestoutside');
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
ylabel('W [S^{2}]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
       num2str(d_feed * 1e6) ' \mum']);

%% IMPULSE RESPONSE
h = compute_h_td(tres, f, H);
w = 2 * real(eval_IFT(tres, f, W));

figure('Position', [250 250 950 450]);

subplot(2, 1, 1);
plot(tres * 1e12, permute(h(1, 1, :), [1 3 2]), 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', 'h_{11} = h_{22}');
hold on;
plot(tres * 1e12, permute(h(1, 2, :), [1 3 2]), 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', 'h_{12} = h_{21}');

grid on;
xlim([0 4]);
legend('location', 'bestoutside');

ylabel('h [S/s]');

subplot(2, 1, 2);
plot(tres * 1e12, w, 'LineWidth', 1.5, 'Color', '#77AC30', ...
     'DisplayName', 'w');

grid on;
xlim([0 4]);
legend('location', 'bestoutside');

ylabel('w [S^{2}/s]');
xlabel('t [ps]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
        num2str(d_feed * 1e6) ' \mum'], 'FontSize', 11, 'FontWeight', 'bold');

%% TIME-STEP ALGORITHM
time_step = TimeStep1(tsim);

% Set time-step algorithm parameters
time_step.K       = calculate_K(laser, pcm);
time_step.Vb      = Vb;
time_step.w       = w;
time_step.h       = h;
time_step.tau_c   = pcm.tau_rec;
time_step.tau_s   = pcm.tau_s;
time_step.tau_d   = tau_d;
time_step.sigma_t = laser.sigma_t;

% Voltage and currents
[v, i] = time_step.simulate();

laser_envelope = exp(- 0.5 * ( (tsim - tau_d') / laser.sigma_t) .^ 2) * laser.P0;

figure('Position', [250 250 1050 550]);

subplot(3, 1, 1);
plot(tsim * 1e12, laser_envelope(1, :) * 1e-3, 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', 'P_{L1}');
hold on;
plot(tsim * 1e12, laser_envelope(2, :) * 1e-3, '--', 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', 'P_{L2}');

grid on;
xlim([-0.3 2.5]);
ylim([0 1.5]);
legend('location', 'bestoutside');

ylabel('P_{L} [kW]');

subplot(3, 1, 2);
plot(tsim * 1e12, v(1, :), 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'v_{1}');
hold on;
plot(tsim * 1e12, v(2, :), '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'v_{2}');

grid on;
xlim([-0.3 2.5]);
% ylim([-1 8]);
legend('location', 'bestoutside');

ylabel('v [V]');

subplot(3, 1, 3);
plot(tsim * 1e12, i(1, :), 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'i_{1}');
hold on;
plot(tsim * 1e12, i(2, :), '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'i_{2}');

grid on;
xlim([-0.3 2.5]);
ylim([-0.1 0.4]);
legend('location', 'bestoutside');

ylabel('i [A]');
xlabel('t [ps]');
sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
        num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], 'FontSize', 11, ...
        'FontWeight', 'bold');

%% TOTAL ENERGY
p = v .* i;
E = sum(p, 2) * dt;

fig = figure('Position', [250 250 950 450]);

plot(tsim * 1e12, p(1, :), 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', ['p_{1}, E_{1} = ' num2str(round(E(1) * 1e12, 2)) ' pJ']);
hold on;
plot(tsim * 1e12, p(2, :), '--', 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', ['p_{2}, E_{2} = ' num2str(round(E(2) * 1e12, 2)) ' pJ']);

grid on;
xlim([-0.3 2.5]);
ylim([-0.5 3]);
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('p [W]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
       num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs']);

%% VOLTAGE CONTRIBUTIONS
W_struct = struct('w', w, 'W', W);

v_comp          = NaN(2, Nt, 2);
v_comp(1, :, 1) = eval_vx(i(1, :), 'h', struct('h', permute(h(1, 1, :), [1 3 2])), 'w', W_struct);
v_comp(2, :, 1) = eval_vx(i(1, :), 'h', struct('h', permute(h(2, 1, :), [1 3 2])), 'w', W_struct);
v_comp(1, :, 2) = eval_vx(i(2, :), 'h', struct('h', permute(h(1, 2, :), [1 3 2])), 'w', W_struct);
v_comp(2, :, 2) = eval_vx(i(2, :), 'h', struct('h', permute(h(2, 2, :), [1 3 2])), 'w', W_struct);

figure('Position', [250 250 1050 550]);

subplot(2, 2, 1);
plot(tsim * 1e12, v_comp(1, :, 1), 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'v_{1}^{s}');

grid on;
xlim([-0.3 2.5]);
ylim([-3 15]);

ylabel('v [V]');
title("$$ \int_{-\infty}^{t}i_{1}(t^{'})z_{11}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

subplot(2, 2, 2);
plot(tsim * 1e12, v_comp(1, :, 2), '--', 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'v_{1}^{m}');

grid on;
xlim([-0.3 2.5]);
ylim([-3 15]);

title("$$ \int_{-\infty}^{t}i_{2}(t^{'})z_{12}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

subplot(2, 2, 3);
plot(tsim * 1e12, v_comp(2, :, 1), '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'v_{2}^{s}');

grid on;
xlim([-0.3 2.5]);
ylim([-3 15]);

ylabel('v [V]');
xlabel('t [ps]');
title("$$ \int_{-\infty}^{t}i_{1}(t^{'})z_{21}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

subplot(2, 2, 4);
plot(tsim * 1e12, v_comp(2, :, 2), 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'v_{2}^{m}');

grid on;
xlim([-0.3 2.5]);
ylim([-3 15]);

xlabel('t [ps]');
title("$$ \int_{-\infty}^{t}i_{2}(t^{'})z_{22}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
        num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], 'FontSize', 11, ...
        'FontWeight', 'bold');

%% INSTANTANEOUS POWER COMPONENTS
p_comp = NaN(2, Nt, 2);
p_comp(1, :, 1) = v_comp(1, :, 1) .* i(1, :);
p_comp(1, :, 2) = v_comp(1, :, 2) .* i(1, :);
p_comp(2, :, 1) = v_comp(2, :, 1) .* i(2, :);
p_comp(2, :, 2) = v_comp(2, :, 2) .* i(2, :);

E_comp = permute(sum(p_comp, 2) * dt, [1 3 2]);

figure('Position', [250 250 1250 550]);

subplot(2, 2, 1);
plot(tsim * 1e12, p_comp(1, :, 1), 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', ['p_{in}, E_{in} = ' num2str(round(E_comp(1, 1) * 1e12, 2)) ' pJ']);

grid on;
xlim([-0.3 2.5]);
ylim([-0.5 3]);
legend('location', 'bestoutside');

ylabel('p [W]');
title("$$ i_{1}(t)\int_{-\infty}^{t}i_{1}(t^{'})z_{11}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

subplot(2, 2, 2);
plot(tsim * 1e12, p_comp(1, :, 2), '--', 'LineWidth', 1.5, 'Color', '#0072BD', ...
     'DisplayName', ['p_{1}^{ind}, E_{1}^{ind} = ' num2str(round(E_comp(1, 2) * 1e12, 2)) ' pJ']);

grid on;
xlim([-0.3 2.5]);
ylim([-0.5 3]);
legend('location', 'bestoutside');

title("$$ i_{1}(t)\int_{-\infty}^{t}i_{2}(t^{'})z_{12}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

subplot(2, 2, 3);
plot(tsim * 1e12, p_comp(2, :, 1), '--', 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', ['p_{2}^{ind}, E_{2}^{ind} = ' num2str(round(E_comp(2, 1) * 1e12, 2)) ' pJ']);

grid on;
xlim([-0.3 2.5]);
ylim([-0.5 3]);
legend('location', 'bestoutside');

ylabel('p [W]');
xlabel('t [ps]');
title("$$ i_{2}(t)\int_{-\infty}^{t}i_{1}(t^{'})z_{21}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

subplot(2, 2, 4);
plot(tsim * 1e12, p_comp(2, :, 2), 'LineWidth', 1.5, 'Color', '#A2142F', ...
     'DisplayName', ['p_{m}, E_{m} = ' num2str(round(E_comp(2, 2) * 1e12, 2)) ' pJ']);

grid on;
xlim([-0.3 2.5]);
ylim([-0.5 3]);
legend('location', 'bestoutside');

xlabel('t [ps]');
title("$$ i_{2}(t)\int_{-\infty}^{t}i_{2}(t^{'})z_{22}(t-t^{'})dt^{'} $$", 'Interpreter', 'latex', ...
      'FontWeight', 'bold');

sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
        num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], 'FontSize', 11, ...
        'FontWeight', 'bold');

%% WAVE ALONG SLOT
% vx                   = NaN(Nx, Nt, 2);
% [vx(:, :, 1), ~, Wx] = eval_vx(i(1, :), 'x', x, 't_res', tres, 'f', f, 'slot', slot, 'x_feed', - d_feed / 2);
% [vx(:, :, 2), ~, ~]  = eval_vx(i(2, :), 'x', x, 't_res', tres, 'f', f, 'slot', slot, 'x_feed', d_feed / 2, ...
%                                'w', Wx);
% 
% figure('Position', [250 250 850 550]);
% 
% for idx = 1 : Nt_plot
%     subplot(Nt_plot, 1, idx);
%     plot(x * 1e6, vx(:, t_idx(idx), 1)', 'LineWidth', 1.5, 'Color', '#0072BD');
%     hold on;
%     plot(- d_feed * 1e6 / 2, 0, 'm*', 'LineWidth', 3.0);
%     hold on;
%     plot(d_feed * 1e6 / 2,  0, 'm*', 'LineWidth', 3.0);
%     
%     grid on;
%     ylim([-3 12.5]);
%     yticks(-3 : 3 : 12);
%     
%     ylabel('v_{s} [V]');
%     title(['@ t = ' num2str(t_plot(idx) * 1e12) ' ps']);
% end
% 
% xlabel('x [\mum]');
% sgtitle(['SELF VOLTAGE @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ...
%          ' \mum, d_{x} = ' num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], ...
%          'FontSize', 11, 'FontWeight', 'bold');
% 
% figure('Position', [250 250 850 550]);
% 
% for idx = 1 : Nt_plot
%     subplot(Nt_plot, 1, idx);
%     plot(x * 1e6, vx(:, t_idx(idx), 2)', 'LineWidth', 1.5, 'Color', '#A2142F');
%     hold on;
%     plot(- d_feed * 1e6 / 2, 0, 'm*', 'LineWidth', 3.0);
%     hold on;
%     plot(d_feed * 1e6 / 2,  0, 'm*', 'LineWidth', 3.0);
%     
%     grid on;
%     ylim([-3 12.5]);
%     yticks(-3 : 3 : 12);
%     
%     ylabel('v_{m} [V]');
%     title(['@ t = ' num2str(t_plot(idx) * 1e12) ' ps']);
% end
% 
% xlabel('x [\mum]');
% sgtitle(['MUTUAL VOLTAGE @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ...
%          ' \mum, d_{x} = ' num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], ...
%          'FontSize', 11, 'FontWeight', 'bold');
% 
% figure('Position', [250 250 950 550]);
% 
% for idx = 1 : 2
%     subplot(2, 1, idx);
%     plot(tsim * 1e12, vx(x_idx(idx), :, 1), 'LineWidth', 1.5, 'Color', '#0072BD', ...
%          'DisplayName', 'v_{s}');
%     hold on;
%     plot(tsim * 1e12, vx(x_idx(idx), :, 2), '--', 'LineWidth', 1.5, 'Color', '#A2142F', ...
%          'DisplayName', 'v_{m}');
%     
%     grid on;
%     xlim([-0.3 2.5]);
%     ylim([-3 12.5]);
%     yticks(-3 : 3 : 12);
%     legend('location', 'bestoutside');
%     
%     ylabel('v [V]');
%     title(['FEED ' num2str(idx)]);
% end
% 
% xlabel('t [ps]');
% sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d_{x} = ' ...
%         num2str(d_feed * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], 'FontSize', 11, ...
%         'FontWeight', 'bold');

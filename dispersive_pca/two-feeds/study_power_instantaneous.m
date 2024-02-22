close all;
clear;
clc;

env;
sim_params;

%% PARAMETERS
f       = linspace(1 * 1e6, 2 * 1e12, 4001);  % START FREQUENCY INFLUENCES RESULTS
df      = f(2) - f(1);
dt      = 1 * 1e-15;
[t, tt] = create_t_vec(dt, -0.5 * 1e-12, 5 * 1e-12);
% FEED
d_gap = 10 * 1e-6;
d     = 200 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% LASER POWER
P = 5 * 1e-3;
% BIAS VOLTAGE
Vb = 30;
% EXCITATION DELAY
tau_d = 1.6 * 1e-12;
% NUMBER OF POINTS
Nf  = length(f);
Nt  = length(t);

%% SLOT, PHOTO-CONDUCTOR, LASER, AND TIME-STEP
slot      = create_slot(d_gap, ws, params_slot);
pcm       = create_pcm(d_gap, ws, params_pcm);
laser     = create_laser(P, params_laser);
time_step = TimeStepAlgorithmTwoFeeds2(t);

%% INPUT IMPEDANCE
Zin = slot.compute_zin(f);

%% MUTUAL IMPEDANCE
Znm = eval_Znm(f, slot, d);

% TIME-DOMAIN
znm = compute_znm(tt, f, Znm);

%% WEIGHT
W = (1 ./ Zin) .^ 2;

% TIME-DOMAIN
w = 2 * real(eval_IFT(tt, f, W));

%% IMPEDANCE WITH WEIGHT
H = eval_Hnm(Znm, W);

% TIME-DOMAIN
h = compute_h(tt, f, H);

%% PHOTO-CONDUCTIVE ANTENNA
pca = PhotoConductiveAntenna(laser, pcm, Vb, h{1,1}, 't_vec', t, 'eta_opt', 1);

%% TIME-STEP
% Set time-step algorithm parameters
time_step.K       = pca.K;
time_step.Vb      = Vb;
time_step.h       = h;
time_step.w       = w;
time_step.tau_c   = pcm.tau_rec;
time_step.tau_s   = pcm.tau_s;
time_step.sigma_t = laser.sigma_t;
time_step.tau_d   = [tau_d 0];

% Voltage and currents
v     = NaN(2, Nt);
i     = v;
i_int = v;
vs    = v;
vm    = v;

for m = 1 : 1 : Nt
    [v(:, m), ~, i(:, m), i_int(:, m), vs(:, m), vm(:, m)] = step(time_step);
end

% PLOT
figure('Position', [250 250 1200 500]);

subplot(2, 2, 1);
plot(t * 1e12, vs(1, :), 'LineWidth', 1.5, 'DisplayName', 'v_{s}');
hold on;
plot(t * 1e12, vm(1, :), 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'v_{m}');
hold on;
plot(t * 1e12, vs(1, :) + vm(1, :), '--', 'LineWidth', 1.5, 'DisplayName', 'v');

grid on;
xlim([-0.5 5]);
ylim([-1 10]);
legend('Location', 'bestoutside');

ylabel('v [V]');
title('FEED 1');

subplot(2, 2, 3);
plot(t * 1e12, i(1, :), 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'i_{1}');
hold on;
xline(tau_d * 1e12, 'm--', 'LineWidth', 1.5, ...
      'DisplayName', ['\tau_{d} = ' num2str(tau_d * 1e12) ' ps']);

grid on;
xlim([-0.5 5]);
legend('Location', 'bestoutside');

ylabel('i [A]');
xlabel('t [ps]');

subplot(2, 2, 2);
plot(t * 1e12, vs(2, :), 'LineWidth', 1.5, 'DisplayName', 'v_{s}');
hold on;
plot(t * 1e12, vm(2, :), 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'v_{m}');
hold on;
plot(t * 1e12, vs(2, :) + vm(2, :), '--', 'LineWidth', 1.5, 'DisplayName', 'v');

grid on;
xlim([-0.5 5]);
ylim([-1 10]);
legend('Location', 'bestoutside');

ylabel('v [V]');
title('FEED 2');

subplot(2, 2, 4);
plot(t * 1e12, i(2, :), 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'i_{2}');
hold on;
xline(time_step.tau_d(2) * 1e12, 'm--', 'LineWidth', 1.5, ...
      'DisplayName', ['\tau_{d} = ' num2str(time_step.tau_d(2) * 1e12) ' ps']);

grid on;
xlim([-0.5 5]);
legend('Location', 'bestoutside');

ylabel('i [A]');
xlabel('t [ps]');

sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
         num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW'], 'FontSize', 11, 'FontWeight', 'bold');

%% INSTANTANEOUS POWER
[ps, pm, pl, pdiss, pt, ~] = eval_p_td(v, vm, vs, i, i_int, Vb);

% PLOT
figure('Position', [250 250 1200 500]);

subplot(2, 2, 1);
plot(t * 1e12, ps(1, :), 'LineWidth', 1.5, 'DisplayName', 'self');
hold on;
plot(t * 1e12, pm(1, :), 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'mutual');

grid on;
xlim([-0.5 5]);
ylim([0 1.2]);
yticks(0 : 0.2 : 1.2);
legend('Location', 'bestoutside');

ylabel('p [W]');
title('FEED 1');

subplot(2, 2, 3);
plot(t * 1e12, pl(1, :), 'LineWidth', 1.5, 'Color', '#EDB120', 'DisplayName', 'radiated');
hold on;
plot(t * 1e12, pdiss(1, :), 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', 'dissipated');
hold on;
plot(t * 1e12, pt(1, :), 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'total');

grid on;
xlim([-0.5 5]);
ylim([0 4.5]);
yticks(0 : 0.5 : 4.5);
legend('Location', 'bestoutside');

ylabel('p [W]');
xlabel('t [ps]');

subplot(2, 2, 2);
plot(t * 1e12, ps(2, :), 'LineWidth', 1.5, 'DisplayName', 'self');
hold on;
plot(t * 1e12, pm(2, :), 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'mutual');

grid on;
xlim([-0.5 5]);
ylim([0 1.2]);
yticks(0 : 0.2 : 1.2);
legend('Location', 'bestoutside');

ylabel('p [W]');
title('FEED 2');

subplot(2, 2, 4);
plot(t * 1e12, pl(2, :), 'LineWidth', 1.5, 'Color', '#EDB120', 'DisplayName', 'radiated');
hold on;
plot(t * 1e12, pdiss(2, :), 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', 'dissipated');
hold on;
plot(t * 1e12, pt(2, :), 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'total');

grid on;
xlim([-0.5 5]);
ylim([0 4.5]);
yticks(0 : 0.5 : 4.5);
legend('Location', 'bestoutside');

ylabel('p [W]');
xlabel('t [ps]');

sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
         num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW'], 'FontSize', 11, 'FontWeight', 'bold');

% COMPARISON PLOT
figure('Position', [250 250 750 450]);

subplot(2, 1, 1);
plot(t * 1e12, ps(1, :), 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'self, FEED 1');
hold on;
plot(t * 1e12, ps(2, :), '--', 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'self, FEED 2');
hold on;
plot(t * 1e12, pm(1, :), 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'mutual, FEED 1');

grid on;
xlim([-0.5 5]);
ylim([0 1.2]);
yticks(0 : 0.2 : 1.2);
legend('Location', 'bestoutside');

ylabel('p [W]');
title('SELF & MUTUAL POWER');

subplot(2, 1, 2);
plot(t * 1e12, pl(1, :), 'LineWidth', 1.5, 'Color', '#EDB120', 'DisplayName', 'FEED 1');
hold on;
plot(t * 1e12, pl(2, :), '--', 'LineWidth', 1.5, 'Color', '#EDB120', 'DisplayName', 'FEED 2');
hold on;
yline(max(pl(1, :)), 'm:', 'LineWidth', 1.5, 'DisplayName', 'MAX @ FEED 1');

grid on;
xlim([-0.5 5]);
ylim([0 1.2]);
yticks(0 : 0.2 : 1.2);
legend('Location', 'bestoutside');

ylabel('p [W]');
xlabel('t [ps]');
title('RADIATED / LOAD POWER');

sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
         num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW'], 'FontSize', 11, 'FontWeight', 'bold');

figure('Position', [250 250 750 450]);

subplot(2, 1, 1);
plot(t * 1e12, pdiss(1, :), 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', 'FEED 1');
hold on;
plot(t * 1e12, pdiss(2, :), '--', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', 'FEED 2');

grid on;
xlim([-0.5 5]);
ylim([0 3.5]);
yticks(0 : 0.5 : 3.5);
legend('Location', 'bestoutside');

ylabel('p [W]');
title('DISSIPATED POWER');

subplot(2, 1, 2);
plot(t * 1e12, pt(1, :), 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'FEED 1');
hold on;
plot(t * 1e12, pt(2, :), '--', 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'FEED 2');

grid on;
xlim([-0.5 5]);
ylim([0 4.5]);
yticks(0 : 0.5 : 4.5);
legend('Location', 'bestoutside');

ylabel('p [W]');
xlabel('t [ps]');
title('TOTAL POWER');

sgtitle(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
         num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW'], 'FontSize', 11, 'FontWeight', 'bold');

%% FUNCTIONS
% FUNCTIONS USED TO CREATE OBJECTS AND VECTORS
function slot = create_slot(d_gap, ws, params_slot)
    slot = SlotInDielectrics(d_gap, ws, params_slot.er_up, params_slot.er_dn);
end

function pcm = create_pcm(d_gap, ws, pcm_config)
    pcm = PhotoConductor([ws d_gap pcm_config.wz], pcm_config.er, ...
                         'tau_rec', pcm_config.tau_rec, 'tau_s', pcm_config.tau_s, ...
                         'me_coef', pcm_config.me_coef, 'absorp_len', pcm_config.alpha);
end

function laser = create_laser(P, laser_config)
    laser = Laser(laser_config.wlen, laser_config.T, P, 'tau_p', laser_config.tau_p, ...
                  'R_3db', laser_config.R_3db);
end

function [t, tt] = create_t_vec(dt, tstart, tend)
    t  = tstart : dt : tend;
    tt = (0 : 1 : length(t) - 1) * dt;
end

% FUNCTIONS USED FOR EVALUATION
function Zm = eval_Zm(f, slot, d)
    c0 = get_phys_const('LightSpeed');

    Nf     = length(f);
    k0     = 2 * pi * f / c0;
    kx_max = slot.get_integration_domain(f);

    Dfs  = @(kx, kfs) slot.compute_D(kx, kfs, 'Top', 'Top');
    feed = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    
    Zm = NaN(1, Nf);
    parfor idx = 1 : Nf
        start_pt = - (kx_max(idx) * k0(idx) + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        D         = @(kx) Dfs(kx, k0(idx));
        integrand = @(kx) (feed(kx) .^ 2) .* exp(- 1j * kx * d) ./ D(kx);

        Zm(idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts);
    end

    Zm = Zm / (2 * pi);
end

function Znm = eval_Znm(f, slot, d)
    Z11 = slot.compute_zin(f);
    Z12 = eval_Zm(f, slot, d);

    Znm{1,1} = Z11;
    Znm{1,2} = Z12;
    Znm{2,1} = Z12;
    Znm{2,2} = Z11;
end

function Hnm = eval_Hnm(Znm, W)
    Hnm{1,1} = Znm{1,1} .* W;
    Hnm{1,2} = Znm{1,2} .* W;
    Hnm{2,1} = Znm{2,1} .* W;
    Hnm{2,2} = Znm{2,2} .* W;
end

function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

function h = compute_h(t, f, Hnm)
    h{1,1} = 2 * real(eval_IFT(t, f, Hnm{1,1}));
    h{1,2} = 2 * real(eval_IFT(t, f, Hnm{1,2}));
    h{2,1} = 2 * real(eval_IFT(t, f, Hnm{2,1}));
    h{2,2} = 2 * real(eval_IFT(t, f, Hnm{2,2}));
end

function znm = compute_znm(t, f, Znm)
    znm11 = 2 * real(eval_IFT(t, f, Znm{1, 1}));
    znm12 = 2 * real(eval_IFT(t, f, Znm{1, 2}));
    znm21 = 2 * real(eval_IFT(t, f, Znm{2, 1}));
    znm22 = 2 * real(eval_IFT(t, f, Znm{2, 2}));

    znm = {znm11, znm12; znm21, znm22};
end

function [ps, pm, pl, pdiss, pt, p_pseudo] = eval_p_td(v, vm, vs, i, i_int, Vb)
    ps       = vs .* i;
    pm       = vm .* i;
    pl       = v .* i;
    pdiss    = (Vb - (vm + vs)) .* i;
    pt       = Vb * i;
    p_pseudo = i_int .* v;
end

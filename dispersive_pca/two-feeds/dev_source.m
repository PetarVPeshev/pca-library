close all;
clear;
clc;

env;
sim_params;

% COLORS = ["#0072BD", "#A2142F", "#EDB120", "#7E2F8E", "#77AC30"];
COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30"];

%% PARAMETERS
f = linspace(50 * 1e9, 2 * 1e12, 4001);
% t = linspace(-0.5, 2.5, 10001) * 1e-12;
dt      = 5 * 1e-15;
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
% NUMBER OF POINTS
Nf  = length(f);
Nt  = length(t);

%% SLOT, PHOTO-CONDUCTOR, LASER, AND TIME-STEP
slot      = create_slot(d_gap, ws, params_slot);
pcm       = create_pcm(d_gap, ws, params_pcm);
laser     = create_laser(P, params_laser);
time_step = TimeStepAlgorithmTwoFeeds2(t);

%% INPUT IMPEDANCE
tic;
Zin = slot.compute_zin(f);
toc

% PLOT - PLOT IS CORRECT
figure('Position', [250 250 750 450]);

plot(f * 1e-12, real(Zin), 'LineWidth', 1.5, 'Color', COLORS(1), 'DisplayName', '\Re');
hold on;
plot(f * 1e-12, imag(Zin), '--', 'LineWidth', 1.5, 'Color', COLORS(2), 'DisplayName', '\Im');

grid on;
legend('location', 'bestoutside');
ylim([0 80]);

xlabel('f [THz]');
ylabel('Z_{in} [\Omega]');
title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum']);

%% MUTUAL IMPEDANCE
tic;
Znm = eval_Znm(f, slot, d);
toc

% PLOT - PLOT LOOKS CORRECT
fig = figure('Position', [250 250 750 450]);

AX1 = axes(fig);
L = gobjects(1, 2);
for idx = 1 : 2
    L(idx) = plot(f * 1e-12, real(Znm{1, idx}), 'LineWidth', 1.5, 'Color', COLORS(idx), 'Parent', AX1);
    hold on;
             plot(f * 1e-12, imag(Znm{1, idx}), '--', 'LineWidth', 1.5, 'Color', COLORS(idx), 'Parent', AX1);
    hold on;
end

LGND1 = legend(AX1, L, {'Z_{11} = Z_{22}', 'Z_{12} = Z_{21}'}, 'location', 'bestoutside');
hold on;
grid on;
    
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
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum']);

%% WEIGHT
W = (1 ./ Zin) .^ 2;

% TIME-DOMAIN
w = 2 * real(eval_IFT(tt, f, W));

% PLOT - PLOT LOOKS CORRECT
figure('Position', [250 250 750 450]);

plot(tt * 1e12, w, 'LineWidth', 1.5, 'DisplayName', 'w');

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('w [S^{2} / s]');
title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum']);

%% IMPEDANCE WITH WEIGHT
H = eval_Hnm(Znm, W);

% TIME-DOMAIN
h = compute_h(tt, f, H);

% PLOT - PLOT LOOKS CORRECT
figure('Position', [250 250 750 450]);

plot(tt * 1e12, h{1,1}, 'LineWidth', 1.5, 'Color', COLORS(1), 'DisplayName', 'h_{1,1} = h_{2,2}');
hold on;
plot(tt * 1e12, h{1,2}, 'LineWidth', 1.5, 'Color', COLORS(2), 'DisplayName', 'h_{1,2} = h_{2,1}');

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('h [S / s]');
title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum']);

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

% Voltage and currents
v     = NaN(2, Nt);
i     = v;
i_int = v;

for m = 1 : 1 : Nt
    [v(:, m), ~, i(:, m), i_int(:, m)] = step(time_step);
end

i_impr = time_step.i_impr;

% CURRENT PLOT - PLOT IS CORRECT
figure('Position', [250 250 750 450]);

plot(t * 1e12, i_impr, 'LineWidth', 1.5, 'Color', COLORS(1), 'DisplayName', 'i_{impr}');
hold on;
plot(t * 1e12, i(1, :), 'LineWidth', 1.5, 'Color', COLORS(2), 'DisplayName', 'i_{1}');
hold on;
plot(t * 1e12, i_int(1, :), 'LineWidth', 1.5, 'Color', COLORS(3), 'DisplayName', 'i_{int,1}');
hold on;
plot(t * 1e12, i(2, :), '--', 'LineWidth', 1.5, 'Color', COLORS(4), 'DisplayName', 'i_{2}');
hold on;
plot(t * 1e12, i_int(2, :), '--', 'LineWidth', 1.5, 'Color', COLORS(5), 'DisplayName', 'i_{int,2}');

grid on;
xlim([-0.5 2.5]);
yticks(0 : 0.05 : 0.2);
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('i [A]');
title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs']);

% VOLTAGE PLOT - PLOT IS CORRECT
figure('Position', [250 250 750 450]);

plot(t * 1e12, v(1, :), 'LineWidth', 1.5, 'Color', COLORS(1), 'DisplayName', 'v_{1}');
hold on;
plot(t * 1e12, v(2, :), '--', 'LineWidth', 1.5, 'Color', COLORS(2), 'DisplayName', 'v_{2}');

grid on;
xlim([-0.5 4]);
ylim([-1 8]);
yticks([-1 0 2 4 6 8]);
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('v [V]');
title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs']);

%% FREQUENCY-DOMAIN
V = NaN(2, Nf);
I = V;
for idx = 1 : 2
    V(idx, :) = eval_FT(t, f, v(idx, :));
    I(idx, :) = eval_FT(t, f, i(idx, :));
end

IZ = [I(1, :) .* Znm{1, 1} + I(2, :) .* Znm{1, 2}; I(1, :) .* Znm{2, 1} + I(2, :) .* Znm{2, 2}];

% PLOT
fig = figure('Position', [250 250 800 450]);

AX1 = axes(fig);
L = gobjects(1, 2);
L(1) = plot(f * 1e-12, real(V(2, :)), 'LineWidth', 1.5, 'Color', COLORS(1), 'Parent', AX1);
hold on;
       plot(f * 1e-12, imag(V(2, :)), '--', 'LineWidth', 1.5, 'Color', COLORS(1), 'Parent', AX1);
hold on;
L(2) = plot(f * 1e-12, real(IZ(2, :)), 'LineWidth', 1.5, 'Color', COLORS(2), 'Parent', AX1);
hold on;
       plot(f * 1e-12, imag(IZ(2, :)), '--', 'LineWidth', 1.5, 'Color', COLORS(2), 'Parent', AX1);
hold on;

LGND1 = legend(AX1, L, {'FT\{ v_{2} \}', 'FT\{ i_{1} \}Z_{21} + FT\{ i_{2} \}Z_{21}'}, 'location', 'bestoutside');
hold on;
grid on;

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
ylabel('V [V / Hz]');
title(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs']);

%% ERROR
err = abs(V) - abs(IZ);
ene = sum(abs(V) .^ 2, 2);
NMSE = sum(err .^ 2, 2) ./ ene;
disp(NMSE);

% PLOT
figure('Position', [250 250 800 450]);

plot(f * 1e-12, err(1, :), 'LineWidth', 1.5, 'DisplayName', 'n = 1');
hold on;
plot(f * 1e-12, err(2, :), '--', 'LineWidth', 1.5, 'DisplayName', 'n = 2');

grid on;
legend('location', 'bestoutside');

xlabel('f [THz]');
ylabel('err [V / Hz]');
title(['|V_{n}| - |I_{1}Z_{n1} + I_{2}Z_{n2}| @ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' ...
       num2str(d_gap * 1e6) ' \mum, d = ' num2str(d * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs']);

%% FUNCTIONS
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

function Y = eval_FT(t, f, y)
    N = length(f);
    dt = t(2) - t(1);

    Y = NaN(1, N);
    for f_idx = 1 : 1 : N
        Y(f_idx) = sum(y .* exp(- 1j * 2 * pi * f(f_idx) * t)) * dt;
    end
end

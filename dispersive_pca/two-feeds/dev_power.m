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

%% MUTUAL IMPEDANCE
tic;
Znm = eval_Znm(f, slot, d);
toc

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
time_step.tau_d   = [0 0];

% Voltage and currents
v  = NaN(2, Nt);
i  = v;
vs = v;
vm = v;

for m = 1 : 1 : Nt
    [v(:, m), ~, i(:, m), ~, vs(:, m), vm(:, m)] = step(time_step);
end

% PLOT
figure('Position', [250 250 750 450]);

plot(t * 1e12, vs(1, :), 'LineWidth', 1.5, 'DisplayName', 'v_{s}');
hold on;
plot(t * 1e12, vm(1, :), 'LineWidth', 1.5, 'DisplayName', 'v_{m}');
hold on;
plot(t * 1e12, v(1, :), '--', 'LineWidth', 1.5, 'DisplayName', 'v');

grid on;
legend('Location', 'bestoutside');

xlabel('t [ps]');
ylabel('v [V]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW']);

% PLOT
figure('Position', [250 250 750 450]);

plot(t * 1e12, v(1, :), 'LineWidth', 1.5, 'DisplayName', 'v');
hold on;
plot(t * 1e12, vs(1, :) + vm(1, :), '--', 'LineWidth', 1.5, 'DisplayName', 'v_{s} + v_{m}');

grid on;
legend('Location', 'bestoutside');

xlabel('t [ps]');
ylabel('v [V]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW']);

%% FREQUENCY-DOMAIN
V = NaN(2, Nf);
I = V;
for idx = 1 : 2
    V(idx, :) = eval_FT(t, f, v(idx, :));
    I(idx, :) = eval_FT(t, f, i(idx, :));
end

%% SELF & MUTUAL VOLTAGES
% IN TIME-DOMAIN
[vs, vm] = eval_vms(i, znm, dt);

% PLOT - INCORRECT
figure('Position', [250 250 750 450]);

plot(t * 1e12, vs(1, :), 'LineWidth', 1.5, 'DisplayName', 'v_{s}');
hold on;
plot(t * 1e12, vm(1, :), 'LineWidth', 1.5, 'DisplayName', 'v_{m}');
hold on;
plot(t * 1e12, v(1, :), '--', 'LineWidth', 1.5, 'DisplayName', 'v');

grid on;
legend('Location', 'bestoutside');

xlabel('t [ps]');
ylabel('v [V]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW']);

% IN FREQUENCY-DOMAIN
[Vs, Vm] = eval_Vms(I, Znm);
% CONVERT TO TIME-DOMAIN
vs = NaN(2, Nt);
vm = vs;
for idx = 1 : 2
    vs(idx, :) = 2 * real(eval_IFT(t, f, Vs(idx, :)));
    vm(idx, :) = 2 * real(eval_IFT(t, f, Vm(idx, :)));
end

% PLOT
figure('Position', [250 250 750 450]);

plot(t * 1e12, vs(1, :), 'LineWidth', 1.5, 'DisplayName', 'v_{s}');
hold on;
plot(t * 1e12, vm(1, :), 'LineWidth', 1.5, 'DisplayName', 'v_{m}');
hold on;
plot(t * 1e12, v(1, :), '--', 'LineWidth', 1.5, 'DisplayName', 'v');

grid on;
legend('Location', 'bestoutside');

xlabel('t [ps]');
ylabel('v [V]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
       num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW']);

%% ENERGY & POWER
[Es, Em, Et] = eval_E(V, I, Znm);
% Power
[Ps, Pm, ~] = eval_P(Es, Em, Et, laser.T, df);

% PLOT DB SCALE
figure('Position', [250 250 750 450]);

plot(f * 1e-12, 10 * log10(Es(1, :)), 'LineWidth', 1.5, 'DisplayName', 'E_{s}');
hold on;
plot(f * 1e-12, 10 * log10(Em(1, :)), 'LineWidth', 1.5, 'DisplayName', 'E_{m}');
hold on;
plot(f * 1e-12, 10 * log10(Es(1, :) + Em(1, :)), '--', 'LineWidth', 1.5, 'DisplayName', 'E_{s} + E_{m}');

grid on;
legend('Location', 'bestoutside');

xlabel('f [THz]');
ylabel('E [dBJ / Hz]');
title(['@ d = ' num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW']);

% PLOT LINEAR SCALE
figure('Position', [250 250 750 450]);

plot(f * 1e-12, Es(1, :), 'LineWidth', 1.5, 'DisplayName', 'E_{s}');
hold on;
plot(f * 1e-12, Em(1, :), 'LineWidth', 1.5, 'DisplayName', 'E_{m}');
hold on;
plot(f * 1e-12, Es(1, :) + Em(1, :), '--', 'LineWidth', 1.5, 'DisplayName', 'E_{s} + E_{m}');

grid on;
legend('Location', 'bestoutside');

xlabel('f [THz]');
ylabel('E [dBJ / Hz]');
title(['@ d = ' num2str(d * 1e6) ' \mum, P = ' num2str(P * 1e3) ' mW']);

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

function Y = eval_FT(t, f, y)
    N = length(f);
    dt = t(2) - t(1);

    Y = NaN(1, N);
    for f_idx = 1 : 1 : N
        Y(f_idx) = sum(y .* exp(- 1j * 2 * pi * f(f_idx) * t)) * dt;
    end
end

function [vs, vm] = eval_vms(i, znm, dt)
    Nt = size(i, 2);

    z11 = znm{1, 1};
    z12 = znm{1, 2};
    z21 = znm{2, 1};
    z22 = znm{2, 2};

    vs1 = NaN(1, Nt);
    vs2 = vs1;
    vm1 = vs1;
    vm2 = vs1;
    parfor m = 1 : Nt
        vs1(m) = sum(i(1, 1 : m) .* fliplr(z11(1 : m)));
        vs2(m) = sum(i(2, 1 : m) .* fliplr(z22(1 : m)));
        vm1(m) = sum(i(2, 1 : m) .* fliplr(z12(1 : m)));
        vm2(m) = sum(i(1, 1 : m) .* fliplr(z21(1 : m)));
    end

    vs = [vs1; vs2] * dt;
    vm = [vm1; vm2] * dt;
end

function [Vs, Vm] = eval_Vms(I, Znm)
    Vs1 = I(1, :) .* Znm{1, 1};
    Vs2 = I(2, :) .* Znm{2, 2};

    Vm1 = I(2, :) .* Znm{1, 2};
    Vm2 = I(1, :) .* Znm{2, 1};

    Vs = [Vs1; Vs2];
    Vm = [Vm1; Vm2];
end

function [Es, Em, Et] = eval_E(V, I, Znm)  % FREQUENCY-DOMAIN
    N  = size(V, 1);  % number of elements
    Nf = size(V, 2);

    Es = NaN(N, Nf);
    Em = Es;
    Et = Es;

    Ss1 = I(1, :) .* (I(1, :) .* Znm{1, 1});
    Ss2 = I(2, :) .* (I(2, :) .* Znm{2, 2});

    Sm1 = I(1, :) .* (I(2, :) .* Znm{1, 2});
    Sm2 = I(2, :) .* (I(1, :) .* Znm{2, 1});

    St = I .* V;

    % self
    Es(1, :) = abs(I(1, :)) .* abs(I(1, :) .* Znm{1, 1});
    Es(2, :) = abs(I(2, :)) .* abs(I(2, :) .* Znm{2, 2});
    % Es(1, :) = I(1, :) .* conj(I(1, :) .* Znm{1, 1});
    % Es(2, :) = I(2, :) .* conj(I(2, :) .* Znm{2, 2});

    % mutual
    Em(1, :) = abs(I(1, :)) .* abs(I(2, :) .* Znm{1, 2});
    Em(2, :) = abs(I(2, :)) .* abs(I(1, :) .* Znm{2, 1});
    % Em(1, :) = I(1, :) .* conj(I(2, :) .* Znm{1, 2});
    % Em(2, :) = I(2, :) .* conj(I(1, :) .* Znm{2, 1});

    % total
    Et(1, :) = abs(I(1, :)) .* abs(V(1, :));
    Et(2, :) = abs(I(2, :)) .* abs(V(2, :));
    % Et(1, :) = I(1, :) .* conj(V(1, :));
    % Et(2, :) = I(2, :) .* conj(V(2, :));

    % for row_num = 1 : N
    %     Et(row_num, :) = abs(I(row_num, :)) .* abs(V(row_num, :));
    % 
    %     for col_num = 1 : N
    %         if row_num == col_num  % self
    %             Es(row_num, :) = (abs(I(row_num, :) .^ 2)) .* real(Znm{row_num, col_num});
    %         else  % mutual
    %             Em(row_num, :) = abs(I(row_num, :)) .* abs(I(col_num, :)) .* real(Znm{row_num, col_num});
    %         end
    %     end
    % end
end

function [Ps, Pm, Pt] = eval_P(Es, Em, Et, T, df)  % FREQUENCY-DOMAIN
    Ps = (2 / T) * sum(Es, 2) * df;
    Pm = (2 / T) * sum(Em, 2) * df;
    Pt = (2 / T) * sum(Et, 2) * df;
end

function [ps, pm, pt] = eval_p_td(vm, vs, i)
    ps1 = vs(1, :) .* i(1, :);
    ps2 = vs(2, :) .* i(2, :);

    pm1 = vm(2, :) .* i(1, :);
    pm2 = vm(1, :) .* i(2, :);

    pt1 = (vm(1, :) + vs(1, :)) .* i(1, :);
    pt2 = (vm(2, :) + vs(2, :)) .* i(2, :);

    ps = [ps1; ps2];
    pm = [pm1; pm2];
    pt = [pt1; pt2];
end

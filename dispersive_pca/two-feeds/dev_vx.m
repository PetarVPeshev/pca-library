close all;
clear;
clc;

env;
sim_params;

COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30"];

t_plts = [121 171 301 401 501 601];

%% PARAMETERS
f       = linspace(50 * 1e9, 2 * 1e12, 1001);
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
% X
dx = d_gap;
x  = (0 : dx : 250 * 1e-6);
% NUMBER OF POINTS
Nf  = length(f);
Nt  = length(t);

%% SLOT, PHOTO-CONDUCTOR, LASER, AND TIME-STEP
slot      = create_slot(d_gap, ws, params_slot);
pcm       = create_pcm(d_gap, ws, params_pcm);
laser     = create_laser(P, params_laser);
time_step = TimeStepAlgorithmTwoFeeds2(t);

%% MUTUAL IMPEDANCE
Znm = eval_Znm(f, slot, d);

%% WEIGHT
W = (1 ./ Znm{1, 1}) .^ 2;

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

% Radiating current
i = NaN(2, Nt);

for m = 1 : 1 : Nt
    [~, ~, i(:, m), ~, ~, ~] = step(time_step);
end

%% IMPEDANCE ALONG SLOT
Zx = evaluate_Zx(x, f, slot, dx, d);

%% WEIGHT
% Wx = evaluate_Wx(f, slot, dx);
Wx = W;

%% IMPEDANCE WITH WEIGHT
Hx = eval_Hx(Zx, Wx);

% TIME-DOMAIN
[hx, wx] = compute_hx_and_wx(tt, f, Hx, Wx);

%% VOLTAGE ALONG SLOT
vx = eval_vx(i, hx, wx);

x_plt  = [- fliplr(x(2 : end)) x];
Nx     = length(x);

% PLOT
figure('Position', [250 250 800 500]);

plt_idx = 1;
for idx = t_plts
    subplot(2, 3, plt_idx)

    vx_plt = [fliplr(vx(2 : Nx, idx)') vx(:, idx)'];
    
    plot(x_plt * 1e6, vx_plt, 'LineWidth', 1.5);
    hold on;
    plot(- d * 1e6 / 2, 0, 'm*', 'LineWidth', 3.0);
    hold on;
    plot(d * 1e6 / 2,  0, 'm*', 'LineWidth', 3.0);
    
    grid on;
    ylim([-2 8]);
    xticks(-250 : 50 : 250);

    xlabel('x [\mum]');
    ylabel('v [V]');
    title(['@ t = ' num2str(t(idx) * 1e12) ' ps']);

    plt_idx = plt_idx + 1;
end

sgtitle(['@ w_{s} = ' num2str(ws(1) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, d = ' ...
        num2str(d * 1e6) ' \mum, \delta_{t} = ' num2str(dt * 1e15) ' fs'], 'FontSize', 11, ...
        'FontWeight', 'bold');

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

function Zx = evaluate_Zx(x, f, slot, dx, d)
    c0 = get_phys_const('LightSpeed');

    Nx = length(x);
    Nf = length(f);
    k0 = 2 * pi * f / c0;
    kx_max = slot.get_integration_domain(f);

    start_pt = - (kx_max .* k0 + 1j * 0.01);
    waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

    Dfs = @(kx, kfs) slot.compute_D(kx, kfs, 'Top', 'Top');
    F   = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx  = @(kx) sinc(- kx * dx / (2 * pi));

    Zx1 = NaN(Nx, Nf);
    Zx2 = Zx1;
    parfor idx = 1 : Nf
        for x_idx = 1 : Nx
            D = @(kx) Dfs(kx, k0(idx));

            integrand_1 = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx) + d / 2)) ./ D(kx);
            Zx1(x_idx, idx) = integral(integrand_1, start_pt(idx), - start_pt(idx), 'Waypoints', waypts) / (2 * pi);
    
            integrand_2 = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx) - d / 2)) ./ D(kx);
            Zx2(x_idx, idx) = integral(integrand_2, start_pt(idx), - start_pt(idx), 'Waypoints', waypts) / (2 * pi);
        end
    end

    Zx = {Zx1, Zx2};
end

function W = evaluate_Wx(f, slot, dx)
    c0 = get_phys_const('LightSpeed');

    Nf = length(f);
    k0 = 2 * pi * f / c0;
    kx_max = slot.get_integration_domain(f);

    Dfs = @(kx, kfs) slot.compute_D(kx, kfs, 'Top', 'Top');
    F   = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx  = @(kx) sinc(- kx * dx / (2 * pi));

    W = NaN(1, Nf);
    parfor idx = 1 : Nf
        start_pt = - (kx_max(idx) * k0(idx) + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        D = @(kx) Dfs(kx, k0(idx));

        integrand = @(kx) F(kx) .* Fx(kx) ./ D(kx);
        W(1, idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
    end
end

function Hx = eval_Hx(Zx, Wx)
    Zx11 = Zx{1, 1};
    Zx12 = Zx{1, 2};
    Hx = {Zx11 .* Wx, Zx12 .* Wx};
end

function [hx, wx] = compute_hx_and_wx(t, f, Hx, Wx)
    Nx = size(Hx{1,1}, 1);
    Nt = length(t);

    wx = 2 * real(eval_IFT(t, f, Wx));

    Hx_11 = Hx{1,1};
    Hx_12 = Hx{1,2};

    hx_11 = NaN(Nx, Nt);
    hx_12 = hx_11;
    parfor idx = 1 : Nx
        hx_11(idx, :) = 2 * real(eval_IFT(t, f, Hx_11(idx, :)));
        hx_12(idx, :) = 2 * real(eval_IFT(t, f, Hx_12(idx, :)));
    end

    hx = {hx_11, hx_12};
end

function vx = eval_vx(i, hx, wx)
    hx11 = hx{1, 1};
    hx12 = hx{1, 2};

    Nt = size(i, 2);
    Nx = size(hx11, 1);

    wx0 = wx(1);

    vx = NaN(Nx, Nt);
    parfor idx = 1 : Nx
        v = NaN(1, Nt);
        for m = 1 : Nt
            Cx1 = sum(i(1, 1 : m) .* fliplr(hx11(idx, 1 : m)));
            Cx2 = sum(i(2, 1 : m) .* fliplr(hx12(idx, 1 : m)));
    
            if m == 1
                Cvx = 0;
            else
                Cvx = sum(v(1 : m - 1) .* fliplr(wx(2 : m)));
            end
    
            v(m) = (Cx1 + Cx2 - Cvx) / wx0;
        end

        vx(idx, :) = v;
    end
end

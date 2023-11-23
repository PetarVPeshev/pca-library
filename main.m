close all;
clear;
clc;

addpath([pwd() '\utils']);

[laser, GaAs] = get_default_config();
[~] = laser.plot_time_envelope();

qo_link = QuasiOpticalLink(laser, GaAs, eta_opt = 0.35);

% march_on = TimeStepAlgorithm(laser.get_t_vec(laser.tau_p / 100));

antenna = PhotoConductiveAntenna(laser, GaAs, march_on, eff_opt = 0.35);

%% PARAMETERS
% LASER
sigma_t = laser.sigma_t;
% PHOTOCONDUCTOR
tau_rec = GaAs.tau_rec;
tau_s = GaAs.tau_s;
% ANTENNA
Vb = 30;
k = antenna.k_const;
Ga = 1 / 70;
% SIMULATION
dt = laser.tau_p / 50;
t_vec = laser.get_t_vec(dt);

%% IMPRESSED CURRENT
i_impr = compute_i_impr(t_vec, k, Vb, tau_s, tau_rec, sigma_t);

%% h_{m}
m_max = length(t_vec);

% h = NaN(1, m_max);
% for m = 1 : 1 : m_max
%     h(m) = compute_hm_sum(t_vec, m, k, tau_rec, sigma_t);
% end

hm = compute_hm(t_vec, k, tau_rec, sigma_t);

% hm_n = compute_hm(t_vec, k, tau_rec, sigma_t);
% hm = sum(hm_n, 2)';

% [v, i_impr, i_int] = compute_v(t_vec, Vb, Ga, k, tau_s, tau_rec, sigma_t);

figure('Position', [250 250 700 400]);
plot(t_vec(2 : end) * 1e12, i_impr, 'LineWidth', 2.0, ...
    'DisplayName', ['i_{impr}, \delta_{t} = ' num2str(dt * 1e15) ' fs']);
% hold on;
% plot(t_vec * 1e12, i_int, 'LineWidth', 2.0, ...
%     'DisplayName', ['i_{int}, \delta_{t} = ' num2str(dt * 1e15) ' fs']);
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('t / ps');
ylabel('i(t) / A');
title(['i @ V_{b} = ' num2str(Vb) ' V, \tau_{rec} = ' num2str(tau_rec * 1e15) ...
    ' fs, \tau_{s} = ' num2str(tau_s * 1e15) ' fs, \sigma_{t} = ' num2str(round(sigma_t * 1e15, 2)) ' fs']);

% figure('Position', [250 250 700 400]);
% plot(t_vec * 1e12, h, 'LineWidth', 2.0, ...
%     'DisplayName', ['h[m], \delta_{t} = ' num2str(dt * 1e15) ' fs']);
% hold on;
% plot(t_vec(2 : end) * 1e12, hm, '--', 'LineWidth', 2.0, ...
%     'DisplayName', ['second method, \delta_{t} = ' num2str(dt * 1e15) ' fs']);
% grid on;
% legend show;
% legend('location', 'bestoutside');
% xlabel('t / ps');
% ylabel('h[m] / S');
% title(['h[m] @ \tau_{rec} = ' num2str(tau_rec * 1e15) ' fs, \tau_{s} = ' num2str(tau_s * 1e15) ...
%     ' fs, \sigma_{t} = ' num2str(round(sigma_t * 1e15, 2)) ' fs']);

% figure('Position', [250 250 700 400]);
% plot(t_vec * 1e12, v, 'LineWidth', 2.0, 'DisplayName', ['v, \delta_{t} = ' num2str(dt * 1e15) ' fs']);
% grid on;
% legend show;
% legend('location', 'bestoutside');
% xlabel('t / ps');
% ylabel('v(t) / V');
% title(['v @ V_{b} = ' num2str(Vb) ' V, \tau_{rec} = ' num2str(tau_rec * 1e15) ...
%     ' fs, \tau_{s} = ' num2str(tau_s * 1e15) ' fs, \sigma_{t} = ' num2str(sigma_t * 1e15) ' fs']);
% 
% figure('Position', [250 250 700 400]);
% plot(t_vec * 1e12, Vb - v, 'LineWidth', 2.0, ...
%     'DisplayName', ['v_{g}, \delta_{t} = ' num2str(dt * 1e15) ' fs']);
% grid on;
% legend show;
% legend('location', 'bestoutside');
% xlabel('t / ps');
% ylabel('v_{g}(t) / V');
% title(['v_{g} @ V_{b} = ' num2str(Vb) ' V, \tau_{rec} = ' num2str(tau_rec * 1e15) ...
%     ' fs, \tau_{s} = ' num2str(tau_s * 1e15) ' fs, \sigma_{t} = ' num2str(sigma_t * 1e15) ' fs']);

%% FUNCTIONS
function i_impr = compute_i_impr(t_vec, vb, k, tau_s, tau_rec, sigma_t)
    m_max = length(t_vec);
    dt = t_vec(2) - t_vec(1);

    [N, M] = meshgrid(1 : 1 : m_max, 2 : 1 : m_max);
    
    hm_n = exp(- 0.5 * (t_vec(N) / sigma_t) .^ 2) .* exp(- (t_vec(M) - t_vec(N)) / tau_rec);

    i_impr = tril( hm_n .* (1 - exp(- (t_vec(M) - t_vec(N)) / tau_s)), 1 );
    i_impr = vb * k * tau_s * dt * sum(i_impr, 2)';
end

function hm_n = compute_hm_n(t_vec, k, tau_rec, sigma_t)
    m_max = length(t_vec);
    dt = t_vec(2) - t_vec(1);

    [N, M] = meshgrid(1 : 1 : m_max, 2 : 1 : m_max);

    hm_n = k * (dt ^ 2) * tril( exp(- 0.5 * (t_vec(N) / sigma_t) .^ 2) ...
        .* exp(- (t_vec(M) - t_vec(N)) / tau_rec), 1 );
end

function hm_sum = compute_hm_sum(t_vec, m, k, tau_rec, sigma_t)
    dt = t_vec(2) - t_vec(1);

    hm_sum = k * (dt ^ 2) * sum( exp(- 0.5 * (t_vec(1 : m) / sigma_t) .^ 2) ...
        .* exp(- (t_vec(m) - t_vec(1 : m)) / tau_rec) );
end

function hm = compute_hm(t_vec, k, tau_rec, sigma_t)
    hm_n = compute_hm_n(t_vec, k, tau_rec, sigma_t);
    hm = sum(hm_n, 2)';
end

function i_int_curr = compute_i_int_step(t_vec, v_curr, hm_curr, i_int_prev, tau_s, tau_rec)
    dt = t_vec(2) - t_vec(1);

    i_int_curr = exp(- dt / tau_rec) * exp(- dt / tau_s) * i_int_prev + v_curr * hm_curr;

%     i_int_curr = vm * compute_hm_sum(t_vec, m, k, tau_rec, sigma_t);
%     i_int_curr = exp(- dt / tau_rec) * exp(- dt / tau_s) * i_int_prev + i_int_curr;
end

function i_int_curr = compute_i_int(t_vec, v, hm_n, tau_s)
    m = find(isnan(v), 1, 'first') - 1;
    
    i_int_curr = 0;
    for n = 1 : 1 : m
        bm_n = sum( exp(- (t_vec(m) - t_vec(n : m)) / tau_s) .* v(n : m) );
        i_int_curr = i_int_curr + sum( hm_n(m - 1, n : m) .* bm_n );
    end
end

function v_curr = compute_v_step(t_vec, i_impr_curr, i_int_prev, hm_curr, gl, tau_rec, tau_s)
    dt = t_vec(2) - t_vec(1);

%     hm_sum = compute_hm_sum(t_vec, m, k, tau_rec, sigma_t);
    v_curr = ( i_impr_curr - exp(- dt / tau_rec) * exp(- dt / tau_s) * i_int_prev ) / (dt * gl + hm_curr);
end

function [v, i_impr, i_int] = compute_v(t_vec, vb, gl, k, tau_s, tau_rec, sigma_t)
    m_max = length(t_vec);

    v = NaN(1, m_max);
    v(1) = 0;

    i_int = NaN(1, m_max);
    i_int(1) = 0;

    i_impr = compute_i_impr(t_vec, vb, k, tau_s, tau_rec, sigma_t);
    hm_n = compute_hm_n(t_vec, k, tau_rec, sigma_t);
    hm = sum(hm_n, 2)';
%     hm = compute_hm(t_vec, k, tau_rec, sigma_t);

    for m = 2 : 1 : m_max
        i_int_prev = i_int(m - 1);
        i_impr_curr = i_impr(m - 1);
        hm_curr = hm(m - 1);

        v(m) = compute_v_step(t_vec, i_impr_curr, i_int_prev, hm_curr, gl, tau_rec, tau_s);
%         i_int(m) = compute_i_int_step(t_vec, v(m), hm_curr, i_int_prev, tau_s, tau_rec);
        i_int(m) = compute_i_int(t_vec, vb - v, hm_n, tau_s);
    end
end

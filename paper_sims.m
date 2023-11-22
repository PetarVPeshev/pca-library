close all;
clear;
clc;

addpath([pwd() '/utils']);

%% PARAMETERS
% LASER
tau_p   = 100 * 1e-15;  % Pulse FWHM
T       = 12.5 * 1e-9;  % Repetition period
wlen    = 780 * 1e-9;   % Wavelength
pwr_avg = 180 * 1e-3;   % Average power
% QO LINK
eta_opt = 0.35;         % Optical efficiency
% PHOTOCONDUCTOR
Wx      = 10 * 1e-6;    % Gap width
Wy      = 10 * 1e-6;    % Length
Wz      = 2 * 1e-6;     % Height
me_GaAs = 0.067;        % Relative carrier effective mass
tau_c   = 300 * 1e-15;  % Carrier recombination time
tau_s   = 8.5 * 1e-15;  % Scattering time
% ANTENNA
Vb      = 30;           % Bias voltage
Gl      = 1 / 70;       % Radiating element admittance
% SIMULATION TIME
t = linspace(- 10 * tau_p, 30 * tau_p, 4001);

%% TEST TIMESTEP OBJECT
[laser, GaAs] = get_default_config();
% default t_vec
N = 4001;
t_vec = linspace(-20 * laser.tau_p, 20 * laser.tau_p, N);
% constants
qe = get_phys_const('ElectronCharge');
h = get_phys_const('PlanckConstant');
% parameters
pwr_opt = 0.35 * laser.P;
K = qe * GaAs.mu_dc * laser.T * pwr_opt ...
    / ( (GaAs.dimensions.Wx ^ 2) * h * laser.freq * GaAs.tau_s * sqrt(2 * pi) * laser.sigma_t);
ga = Gl;
sigma_t = tau_p / sqrt( 8 * log(2) );
% march-on algorithm
march_on = TimeStep(t_vec, K = K, Vb = Vb, ga = ga, tau_c = tau_c, tau_s = tau_s, sigma_t = sigma_t);

tic;
v = NaN(1, N);
i = NaN(1, N);
i_int = NaN(1, N);
for step_idx = 1 : 1 : N
    [v(step_idx), i(step_idx), i_int(step_idx)] = step(march_on);
end
toc

figure('Position', [250 250 1400 700]);
subplot(2,1,1);
yline(Vb, '--', 'LineWidth', 2.0, 'DisplayName', 'V_{b}');
hold on;
plot(t_vec * 1e12, v, 'LineWidth', 2.0, 'DisplayName', 'v');
hold on;
plot(t_vec * 1e12, Vb - v, 'LineWidth', 2.0, 'DisplayName', 'v_{g}');
grid on;
legend('location', 'bestoutside');
ylabel('v / V');
subplot(2,1,2);
plot(t_vec * 1e12, march_on.i_impr, 'LineWidth', 2.0, 'DisplayName', 'i_{impr}');
hold on;
plot(t_vec * 1e12, i, 'LineWidth', 2.0, 'DisplayName', 'i');
hold on;
plot(t_vec * 1e12, i_int, 'LineWidth', 2.0, 'DisplayName', 'i_{int}');
grid on;
legend('location', 'bestoutside');
ylabel('i / A');
xlabel('t / ps');

%% MAIN
pwr_max = pwr_avg * T * sqrt(4 * log(2) / pi) / tau_p;
fprintf('P0 = %.2f kW\n', round(pwr_max * 1e-3, 2));

A = compute_a_const(eta_opt, pwr_avg, T, wlen, tau_p, Wx, Wy, Wz);
fprintf('A = %.2f * 10^38\n', round(A * 1e-38, 2));

i_impr = compute_i_impr(t, Vb, A, me_GaAs, tau_p, tau_c, tau_s, Wx, Wy, Wz);

v        = NaN(1, length(t));
v(1)     = 0;
i        = NaN(1, length(t));
i(1)     = 0;

tic;
for m = 2 : 1 : length(t)
    v(m)     = compute_v_step(t, m, Vb, i(m - 1), Gl, A, me_GaAs, tau_p, tau_c, tau_s, Wx, Wy, Wz);
    i(m)     = compute_i_step(t, m, v, Vb, A, me_GaAs, tau_p, tau_c, tau_s, Wx, Wy, Wz);
end
toc

i_int = i_impr - i; 

%% PLOTS
figure('Position', [250 250 1400 700]);
annotation('textbox', [.91 .60 .3 .3], 'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'String', {'----QO LINK----',['\eta_{opt} = ' num2str(eta_opt)], ...
    ['P_{L} = ' num2str(pwr_avg * 1e3) ' mW'], ['T_{L} = ' num2str(T * 1e9) ' ns'], ...
    ['\tau_{p} = ' num2str(tau_p * 1e15) ' fs'], '------GaAs------', ...
    ['\tau_{s} = ' num2str(tau_s * 1e15) ' fs'], ['\tau_{c} = ' num2str(tau_c * 1e15) ' fs'], ...
    ['W_{x} = ' num2str(Wx * 1e6) ' \mum'], ['W_{y} = ' num2str(Wy * 1e6) ' \mum'], ...
    ['W_{z} = ' num2str(Wz * 1e6) ' \mum']})
sgtitle('Photo-Conductive Antenna Simulations', 'FontWeight', 'bold');

% LASER
subplot(3, 1, 1);
plot(t * 1e12, pwr_max * exp(- 4 * log(2) * (t / tau_p) .^ 2) * 1e-3, 'LineWidth', 2.0, ...
    'DisplayName', ['P_{L}, \delta_{t} = ' num2str((t(2) - t(1)) * 1e15) ' fs']);
grid on;
legend('location', 'bestoutside');
ylabel('P_{L} / kW');

% CURRENTS
subplot(3, 1, 2);
plot(t * 1e12, i_impr, 'LineWidth', 2.0, ...
    'DisplayName', ['i_{impr}, \delta_{t} = ' num2str((t(2) - t(1)) * 1e15) ' fs']);
hold on;
plot(t * 1e12, i_int, 'LineWidth', 2.0, ...
    'DisplayName', ['i_{int}, \delta_{t} = ' num2str((t(2) - t(1)) * 1e15) ' fs']);
hold on;
plot(t * 1e12, i, 'LineWidth', 2.0, ...
    'DisplayName', ['i, \delta_{t} = ' num2str((t(2) - t(1)) * 1e15) ' fs']);
grid on;
ylim([0 3]);
legend('location', 'bestoutside');
xlabel('t / ps');
ylabel('i / A');

% VOLTAGES
subplot(3, 1, 3);
plot(t * 1e12, v, 'LineWidth', 2.0, 'DisplayName', ['v, \delta_t = ' num2str((t(2) - t(1)) * 1e15) ' fs']);
hold on;
yline(Vb, '--', 'LineWidth', 2.0, 'DisplayName', 'V_{b}');
hold on;
plot(t * 1e12, Vb - v, 'LineWidth', 2.0, ...
    'DisplayName', ['v_{g}, \delta_{t} = ' num2str((t(2) - t(1)) * 1e15) ' fs']);
grid on;
ylim([0 35]);
legend('location', 'bestoutside');
xlabel('t / ps');
ylabel('v / V');

%% FUNCTIONS
function a_const = compute_a_const(eta_opt, pwr_avg, period_rep, wlen, tau_p, wx, wy, wz)
    h = get_phys_const('PlanckConstant');
    c = get_phys_const('LightSpeed');

    V = wx * wy * wz;
    pwr_opt = pwr_avg * eta_opt;
    f = c / wlen;

    a_const = pwr_opt * period_rep * sqrt(4 * log(2) / pi) / (h * f * tau_p * V);
end

function i_impr = compute_i_impr(t, vb, a_const, me_pc, tau_p, tau_c, tau_s, wx, wy, wz)
    qe = get_phys_const('ElectronCharge');
    me = get_phys_const('ElectronMass');

    m_max = length(t);
    
    me_pc = me_pc * me;
    dt = t(2) - t(1);

    const = a_const * (qe ^ 2) * wy * wz * tau_s * vb * dt / (me_pc * wx);
    
    i_impr = NaN(1, m_max);
    for m = 1 : 1 : m_max
        i_impr(m) = sum( exp(- 4 * log(2) * (t(1 : m) / tau_p) .^ 2) .* exp(- (t(m) - t(1 : m)) / tau_c) ...
            .* (1 - exp(- (t(m) - t(1 : m)) / tau_s)) );
    end
    i_impr = const * i_impr;
end

function im = compute_i_step(t, m, v, vb, a_const, me_pc, tau_p, tau_c, tau_s, wx, wy, wz)
    qe = get_phys_const('ElectronCharge');
    me = get_phys_const('ElectronMass');
    
    me_pc = me_pc * me;
    dt = t(2) - t(1);

    const = a_const * (qe ^ 2) * wy * wz * (dt ^ 2) / (me_pc * wx);

%     [N, M] = meshgrid(1 : 1 : m, 1 : 1 : m);
%     v_idx = tril(N - M + m);
%     v_idx(v_idx == 0) = 1;
% 
%     int_sum = fliplr( sum(exp(- (t(M) - t(N)) / tau_s) .* (vb - v(v_idx)), 2)' );
%     im = sum( exp(- 4 * log(2) * (t(1 : m) / tau_p) .^ 2) .* exp(- (t(m) - t(1 : m)) / tau_c) .* int_sum );
%     im = const * im;

    im = NaN(1, m);
    for n = 1 : 1 : m
        int_sum = sum( exp(- (t(m) - t(n : m)) / tau_s) .* (vb - v(n : m)) );
        im(n) = exp(- 4 * log(2) * (t(n) / tau_p) ^ 2) * exp(- (t(m) - t(n)) / tau_c) * int_sum;
    end
    im = const * sum(im);
end

function vm = compute_v_step(t, m, vb, i_prev, gl, a_const, me_pc, tau_p, tau_c, tau_s, wx, wy, wz)
    qe = get_phys_const('ElectronCharge');
    me = get_phys_const('ElectronMass');
    
    me_pc = me_pc * me;
    dt = t(2) - t(1);

    const = a_const * (qe ^ 2) * wy * wz * (dt ^ 2) / (me_pc * wx);

    nom = exp(- dt / tau_c) * exp(- dt / tau_s) * i_prev + const * exp(- 4 * log(2) * (t(m) / tau_p) ^ 2) * vb;
    denom = gl + const * exp(- 4 * log(2) * (t(m) / tau_p) ^ 2);

    vm = nom / denom;
end

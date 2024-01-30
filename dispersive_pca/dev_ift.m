close all;
clear;
clc;

addpath('..', '..\slots', '..\utils');

% er_up : medium 2, er_dn : medium 1
config_slot  = struct('er_up', 11.7, 'er_dn', 1);
config_laser = struct('wlen', 780 * 1e-9, 'T', 12.5 * 1e-9, 'tau_p', 100 * 1e-15, 'R_3db', 5 * 1e-6);
config_pcm   = struct('wz', 2 * 1e-6, 'er', 12.96, 'tau_rec', 300 * 1e-15, 'tau_s', 8.5 * 1e-15, ...
                      'me_coef', 0.067, 'alpha', 1 * 1e-6);

ColorLegend = {'I(\omega)', 'V(\omega) / Z(\omega)'};
LineLegend  = {'\Re', '\Im'};
Location    = 'bestoutside';
LineWidth   = 1.5;
Color       = ["#A2142F", "#0072BD"];
XLabel      = 'f [THz]';
YLabel      = '[A/Hz]';
Position    = [680 558 700 420];

%% PARAMETERS
f = (0.05 : 0.005 : 2) * 1e12;
t = (- 0.3 : 0.001 : 2.5) * 1e-12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = [10 2.5] * 1e-6;
% BIAS
Vb = 30;
% OPTICAL POWER
P = [5 10] * 1e-3;
% NUMBER OF POINTS
Nf  = length(f);
Nt  = length(t);
% ADMITTANCE TIME VECTOR
t_ga = (0 : 1 : length(t) - 1) * (t(2) - t(1));

%% CREATE LASER, AND TIME-STEP ALGORITHM OBJECTS
laser     = create_laser_object(config_laser, 'P', P(1));
time_step = TimeStepAlgorithmWeight(t);

%% CURRENTS FREQUENCY DOMAIN
for idx = 1 : 1 : 2
    % Set laser power
    laser.P = P(idx);

    % Create slot, and pcm objects
    slot = SlotInDielectrics(d_gap, ws(idx), config_slot.er_up, config_slot.er_dn);
    pcm  = create_pcm_object(config_pcm, 'd_gap', d_gap, 'ws', ws(idx));

    % Frequency-domain admittance
    Zin = NaN(1, Nf);
    for f_idx = 1 : 1 : Nf
        Zin(f_idx) = slot.compute_zin(f(f_idx));
    end

    % Time-domain admittance
    gin = 2 * real(eval_IFT(t_ga, f, 1 ./ Zin));
    w   = 2 * real(eval_IFT(t_ga, f, (1 ./ Zin) .^ 2));

    % Create photo-conductive antenna object
    pca = PhotoConductiveAntenna(laser, pcm, Vb, gin, 't_vec', t, 'eta_opt', 1);

    % Set time-step algorithm parameters
    time_step.K       = pca.K;
    time_step.Vb      = Vb;
    time_step.ga      = gin;
    time_step.w       = w;
    time_step.tau_c   = pcm.tau_rec;
    time_step.tau_s   = pcm.tau_s;
    time_step.sigma_t = laser.sigma_t;

    % Voltage and currents
    v     = NaN(1, Nt);
    i     = v;
    i_int = v;

    for m = 1 : 1 : Nt
        [v(m), ~, i(m), i_int(m)] = step(time_step);
    end

    i_impr = time_step.i_impr;

    % Frequency domain of transient voltage and current
    V      = eval_FT(t, f, v);
    I      = eval_FT(t, f, i);
    I_int  = eval_FT(t, f, i_int);
    I_impr = eval_FT(t, f, i_impr);
    
    Title = ['@ w_{y} = ' num2str(ws(idx) * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, P_{opt} = ' ...
             num2str(P(idx) * 1e3) ' mW, V_{b} = ' num2str(Vb) ' V'];

    figure('Position', [680 558 700 380]);

    plot(f * 1e-12, 20 * log10(abs(I) / max(abs(I))), 'LineWidth', 1.5, 'Color', "#A2142F", ...
        'DisplayName', 'I(\omega)');
    hold on;
    plot(f * 1e-12, 20 * log10(abs(I_impr) / max(abs(I_impr))), 'LineWidth', 1.5, 'Color', "#0072BD", ...
        'DisplayName', 'I_{impr}(\omega)');
    hold on;
    plot(f * 1e-12, 20 * log10(abs(I_int) / max(abs(I_int))), 'LineWidth', 1.5, 'Color', "#77AC30", ...
        'DisplayName', 'I_{int}(\omega)');

    box off;
    grid on;
    yticks(- 25 : 5 : 0);
    ylim(round([min(20 * log10(abs(I) / max(abs(I)))) max(20 * log10(abs(I) / max(abs(I))))], 0));
    legend('location', 'bestoutside');

    xlabel('f [THz]');
    ylabel('[dBA / Hz]');
    title(Title);

    % Verification
    I_plt = [I; V ./ Zin];
    
    plot_two(f * 1e-12, real(I_plt), imag(I_plt), 'ColorLegend', ColorLegend, 'LineLegend', LineLegend, ...
        'ColorLocation', Location, 'LineLocation', Location, 'LineWidth', LineWidth, 'Color', ...
        Color, 'XLabel', XLabel, 'YLabel', YLabel, 'Title', Title, 'Position', Position);
    grid on;

    % Release and reset algorithm
    release(time_step);
    reset(time_step);
end

%% FUNCTIONS
function Y = eval_FT(t, f, y)
    N = length(f);
    dt = t(2) - t(1);

    Y = NaN(1, N);
    for f_idx = 1 : 1 : N
        Y(f_idx) = sum(y .* exp(- 1j * 2 * pi * f(f_idx) * t)) * dt;
    end
end

function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

% FUNCTIONS USED TO CREATE OBJECTS
function laser = create_laser_object(laser_config, params)
    arguments
        laser_config (1,1) struct
        params.P     (1,1) double
    end

    laser = Laser(laser_config.wlen, laser_config.T, params.P, 'tau_p', laser_config.tau_p, ...
                  'R_3db', laser_config.R_3db);
end

function pcm = create_pcm_object(pcm_config, params)
    arguments
        pcm_config   (1,1) struct
        params.d_gap (1,1) double
        params.ws    (1,1) double
    end

    pcm = PhotoConductor([params.ws params.d_gap pcm_config.wz], pcm_config.er, ...
                         'tau_rec', pcm_config.tau_rec, 'tau_s', pcm_config.tau_s, ...
                         'me_coef', pcm_config.me_coef, 'absorp_len', pcm_config.alpha);
end

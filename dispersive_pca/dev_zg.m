close all;
clear;
clc;

addpath('..');
addpath('..\slots');
addpath('..\utils');

LineLegend = {'\Re', '\Im'};
Location = 'bestoutside';
LineWidth = 1.5;
Color = ["#0072BD", "#EDB120", "#77AC30", "#A2142F"];
Position = [680 558 700 420];

%% PARAMETERS
% f = eps : 1e9 : 2e12;
% t = -0.3e-12 : 0.001e-12 : 3e-12;
% f = 0.05e12 : 5e9 : 2e12;             % Lawrence vectors
% t = -0.3e-12 : 0.0001e-12 : 3e-12;    % Lawrence vectors
t = linspace(-2, 8, 4001) * 1e-12;
f = (0.05 : 0.005 : 2) * 1e12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws_vec = (10 : - 2.5 : 2.5) * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1
% PHOTOCONDUCTOR DIMENSIONS
wx_vec = ws_vec;    % gap between metal
wy = d_gap;         % in PhotoConductiveAntenna the gap between the metal is wx
wz = 2 * 1e-6;
% PHOTOCONDUCTOR RELATIVE PERMITTIVITY
er_pcm = 12.96;
% RECOMBINATION AND SCATTERING TIME
tau_rec = 300 * 1e-15;
tau_s = 8.5 * 1e-15;
% EFFECTIVE ELECTRON MASS COEFFICIENT
me_coef = 0.067;
% ABSORPTION LENGTH
alpha = 1 * 1e-6;
% BIAS
Vb = 30;
% LASER WAVELENGTH
laser_wlen = 780 * 1e-9;
% LASER PERIOD
T = 12.5 * 1e-9;
% OPTICAL POWER
P = 10 * 1e-3;
% LASER FWHM
tau_p = 100 * 1e-15;
% LASER RADIUS HALF WIDTH
R_3db = 5 * 1e-6;

%% IMPEDANCE TIME VECTOR
t_za = t(end) + (1 : 1 : find(t == 0, 1) - 1) * (t(2) - t(1));
t_za = [t t_za];

%% NORTON GENERATOR IMPEDANCE
for wx_idx = 1 : 1 : length(wx_vec)
    slot = SlotInDielectrics(d_gap, ws_vec(wx_idx), er_up, er_dn);

    Zin = NaN(1, length(f));
    for f_idx = 1 : 1 : length(f)
        Zin(f_idx) = slot.compute_zin(f(f_idx));
    end
    zin = 2 * real(eval_IFT(t_za, f, Zin));
    zin_lawrence = IFT_frequency_to_time(Zin, t_za, f);

    % zin(zin < 0) = 0;

    figure('Position', [50 50 700 420]);
    plot(t_za * 1e12, zin, 'LineWidth', 1.5, 'DisplayName', 'IFT, own');
    hold on;
    plot(t_za * 1e12, zin_lawrence, '--', 'LineWidth', 1.5, 'DisplayName', 'IFT, Lawrence');
    grid on;
    legend('location', 'bestoutside');
    xlabel('t [ps]');
    ylabel('z_{in} [\Omega]');
    title(['@ w_{s} = ' num2str(ws_vec(wx_idx) * 1e6) ' \mum']);
    
    pcm = PhotoConductor([wx_vec(wx_idx) wy wz], er_pcm, 'tau_rec', tau_rec, 'tau_s', tau_s, ...
        'me_coef', me_coef, 'absorp_len', alpha);
    laser = Laser(laser_wlen, T, P, 'R_3db', R_3db, 'tau_p', tau_p);

    pca = PhotoConductiveAntenna(laser, pcm, Vb, zin, 'eta_opt', 1, 't_vec', t);
    [v, vg, i_impr, i_int, i] = pca.compute_response();

    figure('Position', [250 250 1400 700]);
    annotation('textbox', [.91 .60 .3 .3], 'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
        'String', {'----QO LINK----',['\eta_{opt} = ' num2str(pca.qo_link.eta_opt)], ...
        ['P_{L} = ' num2str(laser.P * 1e3) ' mW'], ['T_{L} = ' num2str(laser.T * 1e9) ' ns'], ...
        ['\tau_{p} = ' num2str(laser.tau_p * 1e15) ' fs'], '------GaAs------', ...
        ['\tau_{s} = ' num2str(pcm.tau_s * 1e15) ' fs'], ['\tau_{c} = ' num2str(pcm.tau_rec * 1e15) ' fs'], ...
        ['W_{x} = ' num2str(pcm.dimensions.Wx * 1e6) ' \mum'], ...
        ['W_{y} = ' num2str(pcm.dimensions.Wy * 1e6) ' \mum'], ...
        ['W_{z} = ' num2str(pcm.dimensions.Wz * 1e6) ' \mum']})
    sgtitle('Photo-Conductive Antenna Simulations', 'FontWeight', 'bold');
    
    % LASER
    subplot(3, 1, 1);
    
    plot(pca.time_step.t_vec * 1e12, exp(- 0.5 * (pca.time_step.t_vec / laser.sigma_t) .^ 2), ...
        'LineWidth', 2.0, 'DisplayName', ['P_{L}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
    hold on;
    xline(- laser.tau_p * 1e12, '--', 'LineWidth', 2.0, 'Color', [1 0 1], 'DisplayName', '-\tau_{p}');
    hold on;
    xline(laser.tau_p * 1e12, '--', 'LineWidth', 2.0, 'Color', [1 0 1], 'DisplayName', '\tau_{p}');
    
    grid on;
    xlim([-1 2]);
    legend('location', 'bestoutside');
    
    ylabel('P_{L} / kW');
    
    % CURRENTS
    subplot(3, 1, 2);
    
    plot(pca.time_step.t_vec * 1e12, i_impr, 'LineWidth', 2.0, ...
        'DisplayName', ['i_{impr}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
    hold on;
    plot(pca.time_step.t_vec * 1e12, i_int, 'LineWidth', 2.0, ...
        'DisplayName', ['i_{int}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
    hold on;
    plot(pca.time_step.t_vec * 1e12, i, 'LineWidth', 2.0, ...
        'DisplayName', ['i, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
    
    grid on;
    xlim([-1 2]);
    ylim([0 3]);
    legend('location', 'bestoutside');
    
    xlabel('t / ps');
    ylabel('i / A');
    
    % VOLTAGES
    subplot(3, 1, 3);
    
    yline(pca.Vb, '--', 'LineWidth', 2.0, 'DisplayName', 'V_{b}');
    hold on;
    plot(pca.time_step.t_vec * 1e12, v, 'LineWidth', 2.0, ...
        'DisplayName', ['v, \delta_t = ' num2str(pca.time_step.dt * 1e15) ' fs']);
    hold on;
    plot(pca.time_step.t_vec * 1e12, vg, 'LineWidth', 2.0, ...
        'DisplayName', ['v_{g}, \delta_{t} = ' num2str(pca.time_step.dt * 1e15) ' fs']);
    
    grid on;
    xlim([-1 2]);
    ylim([0 35]);
    legend('location', 'bestoutside');
    
    xlabel('t / ps');
    ylabel('v / V');

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

function y = eval_DIFFT(Y)
    N = length(Y);

    Y_temp = NaN(1, 2 * N);
    Y_temp(1 : N) = Y;
    Y_temp(N + 1 : end) = fliplr(Y);

    y = ifft(Y_temp, 'symmetric');
%     y = ifft(Y);
end

function [f] = IFT_frequency_to_time( F, t, freq )
% This function calculates the Fourier transform from frequency to time for
% a signal with a given frequency vector [Hz] and time-domain [s] over which one would
% like to know the Fourier transform. Make sure the frequency points are
% evenly spaced. For this FT only positive frequencies are required which
% is valid for real space-time functions.

    omega = 2 .* pi .* freq;                                        % angular frequency [rad/s]
    d_omega = omega(2) - omega(1);                                  % spacing in time [s]
    if d_omega < 0; disp("Watch definition frequency vector"); end  % if d_omega becomes negative show a message to bring this to the user's attention
    
    f = zeros(1, length(t) );                                       % initialize output vector
    for tt = 1:length(t)                                            % compute FT for each frequency point
        f(tt) = sum( F .* exp( 1i .* omega .* t(tt) ) ) .* d_omega;
        f(tt) = real( f(tt) ) ./ pi;
    end
end

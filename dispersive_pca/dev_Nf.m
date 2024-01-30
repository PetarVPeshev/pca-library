close all;
clear;
clc;

addpath('..\slots', '..\utils');

LineSpec = ["-", "--", ":", "-.", "--+", ":o"];
Color    = {[0 0.4470 0.7410], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], ...
            [0.6350 0.0780 0.1840], [0 0.75 0.75]};

%% PARAMETERS
Nf = [391, 1001, 2001, 4001, 8001, 16001];
t = linspace(-2, 10, 2001) * 1e-12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1
    
%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

gin = NaN(length(Nf), length(t));
w = gin;
for Nf_idx = 1 : 1 : length(Nf)
    % NOTE: THE IMPULSE RESPONSE IS DEPENDENT ON THE MINIMUM FREQUENCY,
    % MAXIMUM FREQUENCY, AND NUMBER OF POINTS
    f = linspace(eps, 2, Nf(Nf_idx)) * 1e12;
%     if Nf_idx == 1  % LAURENCE POINTS
%         f = (0.05 : 0.005 : 2) * 1e12;
%     end

    %% CHECK SAMPLING FREQUENCY
    fs = 1 / (t(2) - t(1));
    if fs <= 2 * f(end)
        error('SamplingFrequency:NyquistFrequencyNotSatisfied', 'The Nyquist frequency is more than 0.5fs.');
    end
    
    %% ADMITTANCE
    Zin = NaN(1, length(f));
    for f_idx = 1 : 1 : length(f)
        Zin(f_idx) = slot.compute_zin(f(f_idx));
    end
    
    Gin = 1 ./ Zin;
    gin(Nf_idx, :) = 2 * real(eval_IFT(t, f, Gin));
    
    %% WEIGHT FUNCTION
    W = Gin .^ 2;
    w(Nf_idx, :) = 2 * real(eval_IFT(t, f, W));
end

figure('Position', [680 558 700 420]);

for Nf_idx = 1 : 1 : length(Nf)
    plot(t * 1e12, gin(Nf_idx, :), LineSpec(Nf_idx), 'LineWidth', 1.5, 'Color', Color{Nf_idx}, ...
        'MarkerIndices', 1 : 100 : length(t), 'DisplayName', ['N_{f} = ' num2str(Nf(Nf_idx))]);
    hold on;
end
hold off;

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('g_{in} [S]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, ' ...
    'f_{i} \approx ' num2str(round(f(1), 2)) ' Hz']);

figure('Position', [680 558 700 420]);

for Nf_idx = 1 : 1 : length(Nf)
    plot(t * 1e12, w(Nf_idx, :), LineSpec(Nf_idx), 'LineWidth', 1.5, 'Color', Color{Nf_idx}, ...
        'MarkerIndices', 1 : 100 : length(t), 'DisplayName', ['N_{f} = ' num2str(Nf(Nf_idx))]);
    hold on;
end
hold off;

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('w [S^{2}]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' \mum, ' ...
    'f_{i} \approx ' num2str(round(f(1), 2)) ' Hz']);

%% FUNCTIONS
function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

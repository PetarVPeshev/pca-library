close all;
clear;
clc;

env;
sim_params;

%% PARAMETERS
f_min = linspace(eps, 50 * 1e9, 100);
f_max = 2 * 1e12;
NFT   = 4001;
% SLOT WIDTH
ws = 10 * 1e-6;
% FEED GAP
d_gap = 4.5 * 1e-6;
% TIME VECTOR
t_min = -2 * 1e-12;
t_max = 10 * 1e-12;
NT    = 2001;

%% TIME POINTS
t = linspace(t_min, t_max, NT);

%% STUDY START FREQUENCY IMPACT
Zin = NaN(length(f_min), NFT);
parfor s_idx = 1 : length(f_min)
    % Slot
    slot = SlotInDielectrics(d_gap, ws, params_slot.er_up, params_slot.er_dn);

    % Frequency points
    f = linspace(f_min(s_idx), f_max, NFT);

    % Frequency-domain impedance
    for f_idx = 1 : NFT
        Zin(s_idx, f_idx) = slot.compute_zin(f(f_idx));
    end
end

% Impulse response admittance
gin = NaN(length(f_min), NT);
parfor s_idx = 1 : length(f_min)
    % Frequency points
    f = linspace(f_min(s_idx), f_max, NFT);

    % Frequency-domain admittance
    Gin = 1 ./ Zin(s_idx, :);

    % Time-domain admittance
    gin(s_idx, :) = 2 * real(eval_IFT(t, f, Gin));
end

% Figures
figure('Position', [680 558 700 420]);

plot(t * 1e12, gin(1, :), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(1) * 1e-9, 2)) ' GHz']);
hold on;
plot(t * 1e12, gin(50, :), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(50) * 1e-9, 2)) ' GHz']);
hold on;
plot(t * 1e12, gin(100, :), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(100) * 1e-9, 2)) ' GHz']);

grid on;
legend('location', 'bestoutside');

xlabel('t [ps]');
ylabel('g_{in} [S/s]');
title(['@ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' ...
      num2str(d_gap * 1e6) ' \mum, NFT = ' num2str(NFT)]);

figure('Position', [680 558 700 420]);

plot(linspace(f_min(1), f_max, NFT) * 1e-12, real(Zin(1, :)), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(1) * 1e-9, 2)) ' GHz']);
hold on;
plot(linspace(f_min(50), f_max, NFT) * 1e-12, real(Zin(50, :)), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(50) * 1e-9, 2)) ' GHz']);
hold on;
plot(linspace(f_min(100), f_max, NFT) * 1e-12, real(Zin(100, :)), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(100) * 1e-9, 2)) ' GHz']);

grid on;
legend('location', 'bestoutside');

xlabel('f [THz]');
ylabel('Z_{in} [S/s]');
title(['\Re @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' ...
      num2str(d_gap * 1e6) ' \mum, NFT = ' num2str(NFT)]);

figure('Position', [680 558 700 420]);

plot(linspace(f_min(1), f_max, NFT) * 1e-12, imag(Zin(1, :)), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(1) * 1e-9, 2)) ' GHz']);
hold on;
plot(linspace(f_min(50), f_max, NFT) * 1e-12, imag(Zin(50, :)), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(50) * 1e-9, 2)) ' GHz']);
hold on;
plot(linspace(f_min(100), f_max, NFT) * 1e-12, imag(Zin(100, :)), 'LineWidth', 1.5, ...
     'DisplayName', ['f_{m} = ' num2str(round(f_min(100) * 1e-9, 2)) ' GHz']);

grid on;
legend('location', 'bestoutside');

xlabel('f [THz]');
ylabel('Z_{in} [S/s]');
title(['\Im @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' ...
      num2str(d_gap * 1e6) ' \mum, NFT = ' num2str(NFT)]);

%% FUNCTIONS
function y = eval_IFT(t, f, Y)
    Nt = length(t);
    df = f(2) - f(1);

    [F, T] = meshgrid(f, t);
    Y = repmat(Y, Nt, 1);
    y = sum(Y .* exp(1j .* 2 .* pi .* F .* T), 2)' * df;
end

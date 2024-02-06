close all;
clear;
clc;

addpath('..\..\slots', '..\..\utils');

%% PARAMETERS
f = linspace(eps, 2, 391) * 1e12;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1
% SLOT POSITIONS
dx = 0.5 * 1e-6;
x  = linspace(0, 10, 101) * 1e-6;
xq = [0.5 1 2] * 1e-6;
% INTEGRATION DOMAIN
C = 100;
% NUMBER OF POINTS
Nf = length(f);
Nx = length(x);

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% FREQUENCY DOMAIN
Zx_1000k0  = NaN(Nx, Nf);
Zx_Ck0 = Zx_1000k0;
for f_idx = 1 : 1 : Nf
    Zx_1000k0(:, f_idx) = evaluate_Zx(x, f(f_idx), slot, dx);
    Zx_Ck0(:, f_idx)    = evaluate_Zx_Ck0(x, f(f_idx), slot, dx, C);
end

xq_idx   = find(ismember(x, xq));
Position = [680 558 750 400];
YLim   = [0 80];

% REAL PART
YLabel = '\Re\{Z_{x}\} [\Omega]';
Title  = ['\Re @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' mum, d_{x} = ' ...
          num2str(dx * 1e6) ' \mum'];
plot_comparison(f, real(Zx_1000k0), real(Zx_Ck0), x, xq_idx, Position, YLim, YLabel, Title, true, C);

% IMAG PART
YLabel = '\Im\{Z_{x}\} [\Omega]';
Title  = ['\Im @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' mum, d_{x} = ' ...
          num2str(dx * 1e6) ' \mum'];
plot_comparison(f, imag(Zx_1000k0), imag(Zx_Ck0), x, xq_idx, Position, YLim, YLabel, Title, false, C);

%% FUNCTIONS
function Zx = evaluate_Zx(x, f, slot, dx)
    c0 = get_phys_const('LightSpeed');
    k0 = 2 * pi * f / c0;

    Nx = length(x);

    D  = @(kx) slot.compute_D(kx, f, 'TopSheet');
    F  = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx = @(kx) sinc(- kx * dx / (2 * pi));

    Zx = NaN(1, Nx);
    for x_idx = 1 : 1 : Nx
        start_pt = - (1000 * k0 + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        integrand = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D(kx);
        Zx(x_idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
    end
end

function Zx = evaluate_Zx_Ck0(x, f, slot, dx, C)
    c0 = get_phys_const('LightSpeed');
    k0 = 2 * pi * f / c0;

    Nx = length(x);

    D  = @(kx) slot.compute_D(kx, f, 'TopSheet');
    F  = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx = @(kx) sinc(- kx * dx / (2 * pi));

    Zx = NaN(1, Nx);
    for x_idx = 1 : 1 : Nx
        start_pt = - (C * k0 + 1j * 0.01);
        end_pt   = - start_pt;
        waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;

        integrand = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D(kx);
        Zx(x_idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
    end
end

function fig = plot_comparison(f, y1, y2, x, xq_idx, Position, YLim, YLabel, Title, zoom_axes, C)
    Color = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"];
    Nf    = length(f);

    fig = figure('Position', Position);
    
    AX1       = axes(fig);
    L         = gobjects(1, length(xq_idx));
    Legend    = cell(1, length(xq_idx));
    color_idx = 1;
    for x_idx = xq_idx
        L(color_idx) = plot(f * 1e-12, y1(x_idx, :), 'LineWidth', 1.5, ...
                            'Color', Color(color_idx), 'Parent', AX1);
        hold on;
                       plot(f * 1e-12, y2(x_idx, :), '--+', 'LineWidth', 1.5, ...
                            'Color', Color(color_idx), 'MarkerIndices', 1 : 20 : Nf, 'Parent', AX1);
        hold on;
    
        Legend(color_idx) = cellstr(['x = ' num2str(x(x_idx) * 1e6) ' \mum']);
        color_idx         = color_idx + 1;
    end
    
    grid on;
    ylim(YLim);
    set(AX1, 'Box', 'off');
    LGND1 = legend(AX1, L, Legend, 'location', 'bestoutside');
    hold on;
    
    AX2 = copyobj(AX1, gcf);
    delete(get(AX2, 'Children'));
    H1 = plot(NaN, 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
    H2 = plot(NaN, '--+', 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
    
    set(AX2, 'Color', 'none', 'XTick', [], 'Visible', 'off');
    LGND2 = legend([H1 H2], {'1000k_{0}', [num2str(C) 'k_{0}']}, 'location', 'bestoutside');
    LGND2.Position(1) = LGND1.Position(1);
    LGND2.Position(2) = LGND1.Position(2) - LGND2.Position(4) - 0.02;
    LGND2.Box   = 'off';
    LGND2.Color = 'none';
    hold off;

    if zoom_axes
        AX3       = axes('Position', [0.38 0.22 0.35 0.35], 'Box', 'off');
        color_idx = 1;
        for x_idx = xq_idx
            plot(f(f >= 1.4e12 & f <= 2e12) * 1e-12, y1(x_idx, f >= 1.4e12 & f <= 2e12), ...
                 'LineWidth', 1.5, 'Color', Color(color_idx), 'Parent', AX3);
            hold on;
            plot(f(f >= 1.4e12 & f <= 2e12) * 1e-12, y2(x_idx, f >= 1.4e12 & f <= 2e12), '--+', ...
                 'LineWidth', 1.5, 'Color', Color(color_idx), ...
                 'MarkerIndices', 7 : 20 : length(f(f >= 1.4e12 & f <= 2e12)), 'Parent', AX3);
            hold on;
        
            color_idx         = color_idx + 1;
        end

        grid on;
        xlim([1.4 2]);
        annotation('rectangle', [0.555 0.75 0.195 0.17], 'LineWidth', 0.75);
    end
    
    xlabel('f [THz]', 'Parent', AX1);
    ylabel(YLabel, 'Parent', AX1);
    title(Title, 'Parent', AX1);
end

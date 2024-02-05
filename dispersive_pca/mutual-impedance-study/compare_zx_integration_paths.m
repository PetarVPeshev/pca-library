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
% xq = [0 0.5 1 1.5 2 2.5 3] * 1e-6;
xq = [0.5 1 2] * 1e-6;
% NUMBER OF POINTS
Nf = length(f);
Nx = length(x);

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% FREQUENCY DOMAIN
Zx_full  = NaN(Nx, Nf);
Zx_path3 = Zx_full;
for f_idx = 1 : 1 : Nf
    Zx_full(:, f_idx)  = evaluate_Zx(x, f(f_idx), slot, dx);
    Zx_path3(:, f_idx) = evaluate_Zx_path3(x, f(f_idx), slot, dx);
end

xq_idx   = find(ismember(x, xq));
Position = [680 558 750 400];

% REAL PART
YLim   = [0 80];
YLabel = '\Re\{Z_{x}\} [\Omega]';
Title  = ['\Re @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' mum, d_{x} = ' ...
          num2str(dx * 1e6) ' \mum'];
plot_comparison(f, real(Zx_full), real(Zx_path3), x, xq_idx, Position, YLim, YLabel, Title, true);

% IMAG PART
YLim   = [0 120];
YLabel = '\Im\{Z_{x}\} [\Omega]';
Title  = ['\Im @ w_{s} = ' num2str(ws * 1e6) ' \mum, \Delta = ' num2str(d_gap * 1e6) ' mum, d_{x} = ' ...
          num2str(dx * 1e6) ' \mum'];
plot_comparison(f, imag(Zx_full), imag(Zx_path3), x, xq_idx, Position, YLim, YLabel, Title, false);

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

function Zx = evaluate_Zx_path3(x, f, slot, dx)
    c0 = get_phys_const('LightSpeed');
    k0 = 2 * pi * f / c0;

    Nx = length(x);

    D_tt = @(kx) slot.compute_D(kx, f, 'TopSheet');
    D_tb = @(kx) slot.compute_D(kx, f, 'BottomSheet');
    D_bb = @(kx) compute_Dbb(slot, kx, f);
    F    = @(kx) sinc(kx * slot.d_gap / (2 * pi));
    Fx   = @(kx) sinc(- kx * dx / (2 * pi));

    Zx = NaN(1, Nx);
    for x_idx = 1 : 1 : Nx
        if x(x_idx) == 0
            start_pt = - (1000 * k0 + 1j * 0.01);
            end_pt = - start_pt;
            waypts = [-(1 + 1j) (1 + 1j)] * 0.01;

            integrand = @(kx) F(kx) .* Fx(kx) ./ D_tt(kx);
            Zx(x_idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        else
            k2 = k0 * sqrt(slot.er_up);
            kxp = slot.find_kxp(f);

            integrand_tt = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D_tt(kx);
            integrand_tb = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D_tb(kx);
            integrand_bb = @(kx) F(kx) .* Fx(kx) .* exp(- 1j * kx * abs(x(x_idx))) ./ D_bb(kx);

            %% integration path 3
            % Bottom-Bottom sheet integral
            start_pt = k0 - 0.01 - 1j * 50 * k0;
            end_pt = k0 - 0.01;
            int_p1 = integral(integrand_bb, start_pt, end_pt) / (2 * pi);

            % Top-Top sheet integral 1
            start_pt = k0 - 0.01;
            end_pt = k0 + 0.01;
            waypts = [(k0 - 0.01 + 1j * 0.01) (k0 + 0.01 + 1j * 0.01)];
            int_p2 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);

            % Top-Bottom sheet integral
            start_pt = k0 + 0.01;
            end_pt = k2 - 0.01;
            waypts = [(k0 + 0.01 - 1j * 50 * k0) (k2 - 0.01 - 1j * 50 * k0) ...
                      (k2 - 0.01 + 1j * (imag(kxp) - 1)) (real(kxp) - 1 + 1j * (imag(kxp) - 1)) ...
                      (real(kxp) - 1 + 1j * (imag(kxp) + 1)) (k2 - 0.01 + 1j * (imag(kxp) + 1))];
            int_p3 = integral(integrand_tb, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);

            % Top-Top sheet integral 2
            start_pt = k2 - 0.01;
            end_pt = k2 + 0.01 - 1j * 50 * k0;
            waypts = [(k2 - 0.01 + 1j * 0.01) (k2 + 0.01 + 1j * 0.01)];
            int_p4 = integral(integrand_tt, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);

            Zx(x_idx) = int_p1 + int_p2 + int_p3 + int_p4;
        end
    end
end

function D = compute_Dbb(slot, kx, f)
    c0 = get_phys_const('LightSpeed');
    eta0 = get_phys_const('VacuumImpedance');
    
    k0 = 2 * pi * f / c0;
    k1 = k0 * sqrt(slot.er_dn);
    k2 = k0 * sqrt(slot.er_up);

    % Solution in bottom-bottom Riemann sheet
    K1 = 1j * sqrt(- k1 ^ 2 + kx .^ 2);
    K2 = 1j * sqrt(- k2 ^ 2 + kx .^ 2);

    D = (K1 .^ 2) .* compute_J0H02(K1 * slot.ws / 4) + (K2 .^ 2) .* compute_J0H02(K2 * slot.ws / 4);
    D = D / (2 * k0 * eta0);
end

function fig = plot_comparison(f, y1, y2, x, xq_idx, Position, YLim, YLabel, Title, zoom_axes)
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
    H2 = plot(NaN, '--', 'LineWidth', 1.5, 'Color', [0 0 0], 'Parent', AX2);
    
    set(AX2, 'Color', 'none', 'XTick', [], 'Visible', 'off');
    LGND2 = legend([H1 H2], {'FULL', 'AROUND CUTS'}, 'location', 'bestoutside');
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
        annotation('line', [0.73 0.75], [0.57 0.75], 'LineWidth', 0.75);
        annotation('line', [0.38 0.555], [0.57 0.75], 'LineWidth', 0.75);
        annotation('rectangle', [0.555 0.75 0.195 0.17], 'LineWidth', 0.75);
    end
    
    xlabel('f [THz]', 'Parent', AX1);
    ylabel(YLabel, 'Parent', AX1);
    title(Title, 'Parent', AX1);
end

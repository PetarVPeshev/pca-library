close all;
clear;
clc;

addpath('..\..\slots', '..\..\utils');

c0 = get_phys_const('LightSpeed');

%% PARAMETERS
f  = [1e-12 2] * 1e12;
k0 = 2 * pi * f / c0;
% FEED GAP
d_gap = 4.5 * 1e-6;
% SLOT WIDTH
ws = 10 * 1e-6;
% DIELECTRICS
er_up = 11.7;  % medium 2
er_dn = 1;     % medium 1
% SLOT POSITIONS
dx = 0.5 * 1e-6;
x  = 0 * 1e-6;
% REAL PART
kx_r = linspace(-8, 0, 2001);
kx_r = [kx_r fliplr(- kx_r)];
% IMAG PART
kx_i = linspace(-1, 0, 4001) * 1j;
% INTEGRATION PATH
int_path_r = [min(kx_r) -0.01 0.01 max(kx_r)];
int_path_i = [-0.01 -0.01 0.01 0.01];
% NUMBER OF POINTS
Nkx_r = length(kx_r);
Nkx_i = length(kx_i);
Nf    = length(f);
Nx    = length(x);

%% SLOT OBJECT
slot = SlotInDielectrics(d_gap, ws, er_up, er_dn);

%% KX
mid_idx = ceil(length(kx_r) / 2) + 1;

[KX_R, KX_I]           = meshgrid(kx_r, kx_i);
KX_I(:, mid_idx : end) = - KX_I(:, mid_idx : end);

%% POLE & INTEGRAND
kxp = NaN(1, Nf);
ZX  =  NaN(Nkx_i, Nkx_r, Nx, Nf);
for f_idx = 1 : 1 : Nf
    KX = (KX_R + KX_I) * k0(f_idx);
    
    % KXP
    kxp(f_idx) = slot.find_kxp(f(f_idx));
    
    % INTEGRAND
    for x_idx = 1 : 1 : Nx
        ZX(:, :, x_idx, f_idx)  = integrand_ZX(KX, x(x_idx), f(f_idx), slot, dx);

        ZX_max = max(abs(ZX(:, :, x_idx, f_idx)), [], 'all');
        
        figure('Position', [680 558 700 420]);
        
        SURF = surface(KX_R, imag(KX_I), abs(ZX(:, :, x_idx, f_idx)) / ZX_max, 'LineStyle', 'none');
        colormap('hot');
        colorbar;
        
        hold on;
        
        ZData = get(SURF, 'ZData');
        set(SURF, 'ZData', ZData - 10);
        
        INT_PATH = plot(int_path_r, int_path_i, 'LineWidth', 3, 'Color', '#00FFFF', 'DisplayName', 'PATH');
        hold on;
        POLE     = plot(real(kxp(f_idx) / k0(f_idx)), imag(kxp(f_idx) / k0(f_idx)), 'x', ...
                        'LineWidth', 3, 'Color', '#D95319', 'MarkerSize', 10, 'DisplayName', 'POLE');
        hold on;
                   plot(real(- kxp(f_idx) / k0(f_idx)), imag(- kxp(f_idx) / k0(f_idx)), 'x', ...
                       'LineWidth', 3, 'Color', '#D95319', 'MarkerSize', 10);
        
        legend([INT_PATH POLE], 'location', 'bestoutside');
        
        xlabel('\Re');
        ylabel('\Im');
    
        [magn_symb, magn_order] = find_magnitude(f(f_idx));
        title("Normalized |Z_{x}'(k_{x})| @ x = " + num2str(x(x_idx) * 1e6) + " \mum, f = " ...
              + num2str(f(f_idx) / magn_order) + " " + magn_symb + "Hz");
    end
end

%% FUNCTIONS
function ZX = integrand_ZX(kx, x, f, slot, dx)
    D  = slot.compute_D(kx, f, 'TopSheet');
    F  = sinc(kx * slot.d_gap / (2 * pi));
    Fx = sinc(- kx * dx / (2 * pi));

    ZX = F .* Fx .* exp(- 1j * kx * abs(x)) ./ D;
end

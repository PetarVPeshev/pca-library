classdef SlotSemiFinite < SlotABC
    %SLOTSEMIFINITE Summary of this class goes here
    %   Detailed explanation goes here
    % FIXME: This class requires refactoring; the class should not work following the SlotABC refactoring 
    
    properties (SetAccess = immutable)
        d_gap   (1,1) double
        ws      (1,1) double
        ds      (1,1) double
        er_up   (1,1) double
        er_dn   (1,1) double
    end
    
    methods
        function obj = SlotSemiFinite(d_gap, ws, ds, er_up, er_dn)
            %SLOTSEMIFINITE Construct an instance of this class
            %   Detailed explanation goes here

            arguments
                d_gap   (1,1) double {mustBeReal, mustBeNonnegative}
                ws      (1,1) double {mustBeReal, mustBePositive}
                ds      (1,1) double {mustBeReal, mustBePositive}
                er_up   (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(er_up, 1)}
                er_dn   (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(er_dn, 1)}
            end
            
            % TODO: Find better solution to adding the utility path in every class.
            add_dir_path('utils');

            obj.d_gap = d_gap;
            obj.ws = ws;
            obj.ds = ds;
            obj.er_up = er_up;
            obj.er_dn = er_dn;
        end
        
        function kxp = get_kxp_init_guess(obj, f)
            %KXP_INIT_GUESS Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');
    
            k0 = 2 * pi * f / c0;
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);
    
            kxp = sqrt( (k1 ^ 2 + k2 ^ 2) / 2 );
        end

        function D = compute_D(obj, kx, f, sheet)
            %COMPUTE_D Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');
            eta0 = get_phys_const('VacuumImpedance');
            
            k0 = 2 * pi * f / c0;
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);

            % Solution in top Riemann sheet
            K1 = - 1j * sqrt(- k1 ^ 2 + kx .^ 2);
            K2 = - 1j * sqrt(- k2 ^ 2 + kx .^ 2);

            if strcmp(sheet, 'BottomSheet')
                K2 = - K2;
            end

            D = (K1 .^ 2) .* compute_J0H02(K1 * obj.ws / 4) + (K2 .^ 2) .* compute_J0H02(K2 * obj.ws / 4);
            D = D / (2 * k0 * eta0);
        end

        function Fd = compute_Fd(obj, kx)
            %COMPUTE_FD Summary of this method goes here
            %   Detailed explanation goes here
            Fd = sinc(kx * obj.d_gap / (2 * pi));
        end

        function Fg = compute_Fg(obj, kx, f)
            %COMPUTE_FG Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');

            wlen = c0 / f;
            g = (5 / 3) * sqrt(obj.ws * wlen);

            Fg = exp(1j * kx * g / 2) .* (besselj(0, kx * g / 2) - 1j * StruveH0(kx * g / 2) ...
                - (2 / pi) * sinc(kx * g / (4 * pi)) .* exp(- 1j * kx * g / 4));
        end

        function z = compute_z_int(obj, kx, f)
            %COMPUTE_Z_INT Summary of this method goes here
            %   Detailed explanation goes here
            D = obj.compute_D(kx, f, 'TopSheet');
            z = obj.compute_Fd(kx) .* obj.compute_Fd(- kx) ./ D;  % sinc is an even function
        end

        function Zdd = compute_Zdd(obj, f)
            %COMPUTE_ZDD Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');
            k0 = 2 * pi * f / c0;

            waypts = obj.waypts_int;
            if isnan(waypts)
                waypts = [-(1 + 1j) (1 + 1j)] * 0.01;
            end

            if isnan(obj.limits_int)
                start_pt = - (100 * k0 + 1j * 0.01);
                end_pt = - start_pt;
            else
                start_pt = obj.limits_int(1);
                end_pt = obj.limits_int(2);
            end

            integrand = @(kx) obj.compute_z_int(kx, f);
            Zdd = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        end

        function zdg = compute_zdg_int(obj, kx, f)
            %COMPUTE_ZDG_INT Summary of this method goes here
            %   Detailed explanation goes here
            D = obj.compute_D(kx, f, 'TopSheet');
            zdg = obj.compute_Fg(kx, f) .* obj.compute_Fd(- kx) .* exp(1j * kx * obj.ds) ./ D;
        end

        function Zdg = compute_Zdg(obj, f)
            %COMPUTE_ZDG Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');
            k0 = 2 * pi * f / c0;

            waypts = obj.waypts_int;
            if isnan(waypts)
                waypts = [-(1 + 1j) (1 + 1j)] * 0.01;
            end

            if isnan(obj.limits_int)
                start_pt = - (100 * k0 + 1j * 0.01);
                end_pt = - start_pt;
            else
                start_pt = obj.limits_int(1);
                end_pt = obj.limits_int(2);
            end

            integrand = @(kx) obj.compute_zdg_int(kx, f);
            Zdg = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        end

        function zgd = compute_zgd_int(obj, kx, f)
            %COMPUTE_ZDG_INT Summary of this method goes here
            %   Detailed explanation goes here
            D = obj.compute_D(kx, f, 'TopSheet');
            zgd = obj.compute_Fd(kx) .* obj.compute_Fg(- kx, f) .* exp(- 1j * kx * obj.ds) ./ D;
        end

        function Zgd = compute_Zgd(obj, f)
            %COMPUTE_ZGD Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');
            k0 = 2 * pi * f / c0;

            waypts = obj.waypts_int;
            if isnan(waypts)
                waypts = [-(1 + 1j) (1 + 1j)] * 0.01;
            end

            if isnan(obj.limits_int)
                start_pt = - (100 * k0 + 1j * 0.01);
                end_pt = - start_pt;
            else
                start_pt = obj.limits_int(1);
                end_pt = obj.limits_int(2);
            end

            integrand = @(kx) obj.compute_zgd_int(kx, f);
            Zgd = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        end

        function zgg = compute_zgg_int(obj, kx, f)
            %COMPUTE_ZGG_INT Summary of this method goes here
            %   Detailed explanation goes here
            D = obj.compute_D(kx, f, 'TopSheet');
            zgg = obj.compute_Fg(kx, f) .* obj.compute_Fg(- kx, f) ./ D;
        end

        function Zgg = compute_Zgg(obj, f)
            %COMPUTE_ZGG Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');
            k0 = 2 * pi * f / c0;

            waypts = obj.waypts_int;
            if isnan(waypts)
                waypts = [-(1 + 1j) (1 + 1j)] * 0.01;
            end

            if isnan(obj.limits_int)
                start_pt = - (100 * k0 + 1j * 0.01);
                end_pt = - start_pt;
            else
                start_pt = obj.limits_int(1);
                end_pt = obj.limits_int(2);
            end

            integrand = @(kx) obj.compute_zgg_int(kx, f);
            Zgg = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        end

        function zin = compute_zin(obj, f)
            %COMPUTE_ZIN Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f (1,1) double
            end

            Zdd = obj.compute_Zdd(f);
            Zdg = obj.compute_Zdg(f);
            Zgd = obj.compute_Zgd(f);
            Zgg = obj.compute_Zgg(f);

            zin = Zdd - Zdg * Zgd / Zgg;
        end

        function ig = compute_ig(obj, f)
            %COMPUTE_IG Summary of this method goes here
            %   Detailed explanation goes here
            Zgd = obj.compute_Zgd(f);
            Zgg = obj.compute_Zgg(f);

            ig = - Zgd / Zgg;
        end

        function v = compute_v_int(obj, kx, x, f)
            %COMPUTE_V_INT Summary of this method goes here
            %   Detailed explanation goes here
            ig = obj.compute_ig(f);

            D = obj.compute_D(kx, f, 'TopSheet');
            v = (obj.compute_Fd(kx) + ig .* obj.compute_Fg(kx, f)) .* exp(- 1j * kx * x) ./ D;
        end
    end
end

function add_dir_path(sub_dir)
%ADD_FOLDER_PATH Summary of this function goes here
%   sub_dir Sub-directory name
    dir_path = pwd();
    num_pca_lib = count(dir_path, '\pca-library');

    lib_path = extractBefore(dir_path, '\pca-library');
    for dir_num = 1 : 1 : num_pca_lib
        lib_path = append(lib_path, '\pca-library');
    end

    sub_dir_path = append(lib_path, '\', sub_dir);
    addpath(sub_dir_path);
end


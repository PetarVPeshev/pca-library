classdef (Abstract) SlotABC < handle
    %SLOT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        delta_deriv (1,1) double {mustBeInRange(delta_deriv, 0, 1, 'exclude-upper')} = 0.002
        kxp_toler   (1,1) double {mustBeInRange(kxp_toler, 0, 1, 'exclude-upper')} = 1e-6
        limits_int  (1,:) double {mustBeLengthOrNan(limits_int, 2)} = NaN
        waypts_int  (1,:) double = NaN
    end

    methods (Abstract)
        kxp = get_kxp_init_guess(obj, f)
        D = compute_D(obj, kx, f, sheet)
        v = compute_v_int(obj, kx, x, f)
        z = compute_z_int(obj, kx, f)
    end
    
    methods
        function Dp = compute_Dp(obj, kx, f, sheet)
            %COMPUTE_DP Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');

            k0 = 2 * pi * f / c0;
            delta_k = k0 * obj.delta_deriv;
    
            D1 = obj.compute_D(kx + delta_k / 2, f, sheet);
            D2 = obj.compute_D(kx - delta_k / 2, f, sheet);
            Dp = (D1 - D2) / delta_k;
        end

        function kxp = find_kxp(obj, f)
            %FIND_KXP Summary of this method goes here
            %   Detailed explanation goes here
            kxp_prev = 0;
            kxp = obj.get_kxp_init_guess(f);

            while abs(kxp - kxp_prev) > obj.kxp_toler
                kxp_prev = kxp;
                
                D = obj.compute_D(kxp_prev, f, 'BottomSheet');
                Dp = obj.compute_Dp(kxp_prev, f, 'BottomSheet');
        
                kxp = kxp_prev - D / Dp;
            end
        end

        function v = compute_v(obj, x, f)
            %COMPUTE_V Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                x (1,:) double
                f (1,1) double
            end

            c0 = get_phys_const('LightSpeed');
            k0 = 2 * pi * f / c0;

            waypts = obj.waypts_int;
            if isnan(waypts)
                waypts = [-(1 + 1j) (1 + 1j)] * 0.01;
            end

            if isnan(obj.limits_int)
                start_pt = - (50 * k0 + 1j * 0.01);
                end_pt = - start_pt;
            else
                start_pt = obj.limits_int(1);
                end_pt = obj.limits_int(2);
            end

            v = NaN(1, length(x));
            for x_idx = 1 : 1 : length(x)
                integrand = @(kx) obj.compute_v_int(kx, x(x_idx), f);
                v(x_idx) = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
            end
        end

        function zin = compute_zin(obj, f)
            %COMPUTE_ZIN Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f (1,1) double
            end

            c0 = get_phys_const('LightSpeed');
            k0 = 2 * pi * f / c0;

            waypts = obj.waypts_int;
            if isnan(waypts)
                waypts = [-(1 + 1j) (1 + 1j)] * 0.01;
            end

            if isnan(obj.limits_int)
                start_pt = - (1000 * k0 + 1j * 0.01);
                end_pt = - start_pt;
            else
                start_pt = obj.limits_int(1);
                end_pt = obj.limits_int(2);
            end

            integrand = @(kx) obj.compute_z_int(kx, f);
            zin = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        end
    end
end


classdef (Abstract) SlotABC < handle
    %SLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        kxp_init_guess
    end

    properties
        delta_deriv (1,1) double {mustBeInRange(delta_deriv, 0, 1, 'exclude-upper')} = 0.002
        delta_int   (1,1) double {mustBeInRange(delta_int, 0, 1, 'exclude-upper')} = 0.01
        kxp_toler   (1,1) double {mustBeInRange(kxp_toler, 0, 1, 'exclude-upper')} = 1e-6
        limits_int  (1,2) double = NaN
        waypts_int  (1,:) double = NaN
    end

    methods (Abstract)
        D = compute_D(obj, kx, f)
        v = compute_v_int(obj, kx, x, f)
        z = compute_z_int(obj, kx, f)
    end
    
    methods
        function Dp = compute_Dp(obj, kx, f)
            %COMPUTE_DP Summary of this method goes here
            %   Detailed explanation goes here
            c0 = get_phys_const('LightSpeed');

            k0 = 2 * pi * f / c0;
            delta_k = k0 * obj.delta_deriv;

    
            D1 = obj.compute_D(kx + delta_k / 2, f);
            D2 = obj.compute_D(kx - delta_k / 2, f);
            Dp = (D1 - D2) / delta_k;
        end

        function kxp = find_kxp(obj, f)
            %FIND_KXP Summary of this method goes here
            %   Detailed explanation goes here
            kxp_prev = 0;
            kxp = obj.kxp_init_guess;

            while abs(kxp - kxp_prev) > obj.kxp_toler
                kxp_prev = kxp;
                
                D = obj.compute_D(kxp_prev, f);
                Dp = obj.compute_Dp(kxp_prev, f);
        
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

            [start_pt, end_pt, waypts] ...
                = obj.verify_int_pts(k0, obj.limits_int(1), obj.limits_int(2), obj.waypts_int, 'v');

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

            [start_pt, end_pt, waypts] ...
                = obj.verify_int_pts(k0, obj.limits_int(1), obj.limits_int(2), obj.waypts_int, 'zin');

            integrand = @(kx) obj.compute_z_int(kx, f);
            zin = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        end
    end

    methods (Access = private, Sealed, Hidden)
        function [start_pt, end_pt, waypts] = verify_int_pts(obj, k0, start_pt, end_pt, waypts, fun_type)
            %SET_INT_PTS Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                k0       (1,1) double
                start_pt (1,1) double
                end_pt   (1,1) double
                waypts   (1,:) double
                fun_type (1,:) char {mustBeMember(fun_type, {'v', 'zin'})} 
            end

            m = 1000;
            if fun_type == 'v'
                m = 50;
            end
            
            if isnan(start_pt)
                start_pt = - m * k0 - 1j * obj.delta_int;
            end

            if isnan(end_pt)
                end_pt = m * k0 + 1j * obj.delta_int;
            end

            if isnan(waypts)
                waypts = [-(1 + 1j) (1 + 1j)] * obj.delta_int;
            end
        end
    end
end


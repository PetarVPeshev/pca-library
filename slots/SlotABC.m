classdef (Abstract) SlotABC < handle
    %SLOT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        delta_deriv (1,1) double {mustBeInRange(delta_deriv, 0, 1, 'exclude-upper')} = 0.002
        kxp_toler   (1,1) double {mustBeInRange(kxp_toler, 0, 1, 'exclude-upper')} = 1e-6
        N           (1,1) double {mustBePositive, mustBeInteger} = 15
    end

    methods (Abstract)
        kxp    = get_kxp_init_guess(obj, k0)
        D      = compute_D(obj, kx, k0, sheet_k1, sheet_k2)
        v      = compute_v_int(obj, x, k0)
        z      = compute_z_int(obj, k0)
        kx_max = get_integration_domain(obj, f, N)
    end
    
    methods
        function Dp = compute_Dp(obj, kx, k0, sheet_k1, sheet_k2)
            %COMPUTE_DP Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                kx       (1,:) double
                k0       (1,1) double {mustBePositive}
                sheet_k1 (1,:) char   {mustBeMember(sheet_k1, {'Bottom', 'Top'})}
                sheet_k2 (1,:) char   {mustBeMember(sheet_k2, {'Bottom', 'Top'})}
            end
            % TODO: update this function to support k0 as a vector, check if parfor can be used

            delta_k = k0 * obj.delta_deriv;
    
            D1 = obj.compute_D(kx + delta_k / 2, k0, sheet_k1, sheet_k2);
            D2 = obj.compute_D(kx - delta_k / 2, k0, sheet_k1, sheet_k2);
            Dp = (D1 - D2) / delta_k;
        end

        function kxp = find_kxp(obj, f)
            %FIND_KXP Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f   (1,:) double {mustBePositive}
            end

            c0 = get_phys_const('LightSpeed');
            k0 = 2 * pi * f / c0;

            % TODO: implement a parfor loop to improve calculation speed
            kxp = obj.get_kxp_init_guess(k0);
            for idx = 1 : Nf
                kxp_prev = 0;
    
                while abs(kxp(idx) - kxp_prev) > obj.kxp_toler
                    kxp_prev = kxp(idx);
                    
                    D  = obj.compute_D(kxp_prev, k0(idx), 'Top', 'Bottom');
                    Dp = obj.compute_Dp(kxp_prev, k0(idx), 'Top', 'Bottom');
            
                    kxp(idx) = kxp_prev - D / Dp;
                end
            end
        end

        function v = compute_v(obj, x, f)
            %COMPUTE_V Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                x (1,:) double {mustBeReal}
                f (1,:) double {mustBePositive}
            end

            c0 = get_phys_const('LightSpeed');
            Nf = length(f);
            Nx = length(x);

            k0         = 2 * pi * f / c0;
            kx_max     = obj.get_integration_domain(f, obj.N);
            integrands = obj.compute_v_int(x, k0);

            v = NaN(Nf, Nx);
            for f_idx = 1 : Nf
                start_pt = - (kx_max(f_idx) * k0(f_idx) + 1j * 0.01);
                end_pt   = - start_pt;
                waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;
    
                parfor idx = 1 : Nx
                    v(f_idx, idx) = integral(integrands{f_idx, idx}, start_pt, end_pt, 'Waypoints', waypts);
                end
            end

            v = v / (2 * pi);
        end

        function zin = compute_zin(obj, f)
            %COMPUTE_ZIN Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f   (1,:) double {mustBePositive}
            end

            c0 = get_phys_const('LightSpeed');
            Nf = length(f);

            k0         = 2 * pi * f / c0;
            kx_max     = obj.get_integration_domain(f, obj.N);
            integrands = obj.compute_z_int(k0);

            zin = NaN(1, Nf);
            parfor idx = 1 : Nf
                start_pt = - (kx_max(idx) * k0(idx) + 1j * 0.01);
                end_pt   = - start_pt;
                waypts   = [-(1 + 1j) (1 + 1j)] * 0.01;
    
                zin(idx) = integral(integrands{idx}, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
            end
        end
    end
end


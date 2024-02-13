classdef SlotInDielectrics < SlotABC
    %SLOTBETWEENDIELECTRICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        er_up (1,1) double
        er_dn (1,1) double
    end

    properties
        d_gap (1,1) double {mustBeReal, mustBeNonnegative}
        ws    (1,1) double {mustBeReal, mustBeNonnegative}
    end
    
    methods
        function obj = SlotInDielectrics(d_gap, ws, er_up, er_dn)
            %SLOTINDIELECTRICS Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                d_gap   (1,1) double
                ws      (1,1) double
                er_up   (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(er_up, 1)}
                er_dn   (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(er_dn, 1)}
            end

            obj.d_gap = d_gap;
            obj.ws    = ws;
            obj.er_up = er_up;
            obj.er_dn = er_dn;
        end

        function kxp = get_kxp_init_guess(obj, k0)
            %KXP_INIT_GUESS Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                k0  (1,:) double {mustBePositive}
            end
            
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);
    
            kxp = sqrt( (k1 .^ 2 + k2 .^ 2) / 2 );
        end

        function D = compute_D(obj, kx, k0, sheet_k1, sheet_k2)
            %COMPUTE_D Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                kx       (1,:) double
                k0       (1,1) double {mustBePositive}
                sheet_k1 (1,:) char   {mustBeMember(sheet_k1, {'Bottom', 'Top'})}
                sheet_k2 (1,:) char   {mustBeMember(sheet_k2, {'Bottom', 'Top'})}
            end
            % TODO: update this function to support k0 as a vector, check if parfor can be used

            eta0 = get_phys_const('VacuumImpedance');
            
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);

            % Solution in top Riemann sheet
            K1 = - 1j * sqrt(- k1 ^ 2 + kx .^ 2);
            K2 = - 1j * sqrt(- k2 ^ 2 + kx .^ 2);

            if strcmp(sheet_k1, 'Bottom')
                K1 = - K1;
            end
            if strcmp(sheet_k2, 'Bottom')
                K2 = - K2;
            end

            D = (K1 .^ 2) .* compute_J0H02(K1 * obj.ws / 4) + (K2 .^ 2) .* compute_J0H02(K2 * obj.ws / 4);
            D = D / (2 * k0 * eta0);
        end

        function v_int = compute_v_int(obj, kx, x, k0)
            %COMPUTE_V_INT Summary of this method goes here
            %   Detailed explanation goes here
            feed  = sinc(kx * obj.d_gap / (2 * pi));
            D     = obj.compute_D(kx, k0, 'Top', 'Top');
            v_int = feed .* exp(- 1j * kx * x) ./ D;
        end

        function z_int = compute_z_int(obj, kx, k0)
            %COMPUTE_Z_INT Summary of this method goes here
            %   Detailed explanation goes here
            feed  = sinc(kx * obj.d_gap / (2 * pi));
            D     = obj.compute_D(kx, k0, 'Top', 'Top');
            z_int = (feed .^ 2) ./ D;
        end

        function kx_max = get_integration_domain(obj, f, N)
            %GET_INTEGRATION_DOMAIN Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f   (1,:) double {mustBePositive}
                N   (1,1) double {mustBePositive} = 15
            end

            c0   = get_phys_const('LightSpeed');
            wlen = c0 ./ f;

            kx_max = N * wlen / obj.d_gap;
        end

        function kxp = evaluate_kxp_analyt(obj, f)
            %EVALUATE_KXP_ANALYT Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f   (1,:) double {mustBePositive}
            end
            
            c0 = get_phys_const('LightSpeed');
            
            k0 = 2 * pi * f / c0;
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);
            
            beta = sqrt( (k1 .^ 2 + k2 .^ 2) / 2 );
            kd   = sqrt( (- k1 .^ 2 + k2 .^ 2) / 2 );

            B1 = compute_J0H02(- 1j * kd * obj.ws / 4);
            B2 = compute_J0H02(kd * obj.ws / 4);

            A1 = obj.compute_A(- 1j * kd * obj.ws / 4);
            A2 = obj.compute_A(kd * obj.ws / 4);

            kxp = beta + ( (kd .^ 2) ./ beta ) .* (B1 - B2) ./ (A1 + A2);
        end

        function Dp = evaluate_Dp_analyt(obj, kx, k0, sheet_k1, sheet_k2)
            %EVALUATE_DP_ANALYT Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                kx       (1,:) double
                k0       (1,1) double {mustBePositive}
                sheet_k1 (1,:) char   {mustBeMember(sheet_k1, {'Bottom', 'Top'})}
                sheet_k2 (1,:) char   {mustBeMember(sheet_k2, {'Bottom', 'Top'})}
            end

            eta0 = get_phys_const('VacuumImpedance');
            
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);

            % Solution in top Riemann sheet
            K1 = - 1j * sqrt(- k1 ^ 2 + kx(idx) .^ 2);
            K2 = - 1j * sqrt(- k2 ^ 2 + kx(idx) .^ 2);

            if strcmp(sheet_k1, 'Bottom')
                K1 = - K1;
            end
            if strcmp(sheet_k2, 'Bottom')
                K2 = - K2;
            end

            Dp = kx .* (obj.compute_A(K1 * obj.ws / 4) + obj.compute_A(K2 * obj.ws / 4));
            Dp = Dp / (2 * k0 * eta0);
        end

        function vlw = compute_vlw(obj, x, f)
            %COMPUTE_VLW Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                x   (1,:) double {mustBeReal}
                f   (1,:) double {mustBePositive}
            end

            c0 = get_phys_const('LightSpeed');
            Nx = length(x);
            Nf = length(f);

            k0  = 2 * pi * f / c0;
            kxp = obj.find_kxp(f);

            feed = sinc(kxp * obj.d_gap / (2 * pi));

            % TODO: implement a parfor loop to improve calculation speed
            vlw = NaN(Nf, Nx);
            for idx = 1 : Nf
                Dp          = obj.compute_Dp(kxp(idx), k0(idx), 'Top', 'Bottom');
                vlw(idx, :) = - 1j * feed(idx) * exp(- 1j * kxp(idx) * abs(x)) / Dp;
            end
        end

        function Zdyn = compute_Zdyn(obj, f)
            %COMPUTE_ZDYN Summary of this method goes here
            %   Detailed explanation goes here
            % FIXME: This function does not work produce correct result; refactoring is required
            c0 = get_phys_const('LightSpeed');
            eta0 = get_phys_const('VacuumImpedance');
            
            k0 = 2 * pi * f / c0;
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);

            start_pt = - (50 * k0 + 1j * 0.01);
            end_pt = - start_pt;
            waypts = [-(1 + 1j) (1 + 1j)] * 0.01;
            
            const = 1 / (2 * k0 * eta0);
            const_qs = pi * obj.ws / (4 * const);

            qs_part = @(kx) const_qs * (sqrt(k1 ^ 2 - kx .^ 2) - sqrt(k2 ^ 2 - kx .^ 2)) / (k1 ^ 2 - k2 ^ 2);
            integrand = @(kx) (1 ./ obj.compute_D(kx, k0, 'Top', 'Top') - qs_part(kx)) ...
                .* (sinc(kx * obj.d_gap / (2 * pi)) .^ 2);

            Zdyn = integral(integrand, start_pt, end_pt, 'Waypoints', waypts) / (2 * pi);
        end

        function Zqs = compute_Zqs(obj, f)
            %COMPUTE_ZQS Summary of this method goes here
            %   Detailed explanation goes here
            % FIXME: This function does not work produce correct result; refactoring is required
            c0 = get_phys_const('LightSpeed');
            eta0 = get_phys_const('VacuumImpedance');
            eug = exp(0.57721566490153286060651209008240243);
            
            k0 = 2 * pi * f / c0;
            k1 = k0 * sqrt(obj.er_dn);
            k2 = k0 * sqrt(obj.er_up);
            
            const = eta0 * k0 * pi * obj.ws / 2;

            Zqs_p = - 1j * ( (k1 ^ 2) * log(k1 / 2) - (k2 ^ 2) * log(k2 / 2) ) / (pi * (k1 ^ 2 - k2 ^ 2));
            Zqs_p = Zqs_p + (1 - 1j * (2 / pi) * (log(eug * obj.d_gap) - 3 / 2)) / 2;
        
            Zqs_pp = besselh(0, 2, k1 * obj.d_gap) - besselh(0, 2, k2 * obj.d_gap) ...
                + 1j * (2 / pi) * log(k1 / k2);
            Zqs_pp = Zqs_pp ./ (2 * (obj.d_gap .^ 2) * (k1 ^ 2 - k2 ^ 2));

            Zqs = const * (Zqs_p + Zqs_pp);
        end

        function [nd, Z0] = evaluate_TL(obj, f)
            %EVALUATE_TL Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f   (1,:) double {mustBePositive}
            end

            c0 = get_phys_const('LightSpeed');
            Nf = length(f);

            k0  = 2 * pi * f / c0;
            kxp = obj.find_kxp(f);

            nd = sinc(kxp * obj.d_gap / (2 * pi));

            % TODO: implement a parfor loop to improve calculation speed
            Z0 = NaN(1, Nf);
            for idx = 1 : Nf
                Dp         = obj.compute_Dp(kxp(idx), k0(idx), 'Top', 'Bottom');
                Z0(1, idx) = - 2j / Dp;
            end
        end

        function [ZTL, Zrem] = evaluate_EC(obj, f)
            %EVALUATE_EC Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                f   (1,:) double {mustBePositive}
            end

            [nd, Z0] = obj.evaluate_TL(f);
            ZTL = (nd .^ 2) .* Z0 / 2;

            Zin = obj.compute_zin(f);
            Zrem = Zin - ZTL;
        end
    end

    methods (Static)
        function A = compute_A(z)
            %COMPUTE_B Summary of this method goes here
            %   Detailed explanation goes here
            B = compute_J0H02(z);
            A = - 2 * B + z .* besselj(1, z) .* besselh(0, 2, z) + z * besselh(1, 2, z) .* besselj(0, z);
        end
    end
end

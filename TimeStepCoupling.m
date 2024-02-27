classdef TimeStepCoupling < handle
    % TimeStepAlgorithmCoupling Add summary here
    %   Detailed explanation goes here

    properties (SetAccess = immutable)
        t           (1,:) double
        dt          (1,1) double
        m_max       (1,1) double
    end

    properties
        K           (1,1)   double
        Vb          (2,1)   double
        h           (2,2,:) double
        tau_c       (1,1)   double
        tau_s       (1,1)   double
        tau_d       (1,2)   double
        sigma_t     (1,1)   double
    end

    properties (SetAccess = protected)
        Y           (2,2) double
        H           (2,2) double
        Fsm         (2,:) double
        alpha       (1,1) double
    end

    methods
        function obj = TimeStepCoupling(t, params)
            % Support name-value pair arguments when constructing object
            %   Detailed explanation goes here

            arguments
                t               (1,:)   double
                params.K        (1,1)   double = 0
                params.Vb       (2,1)   double = [0; 0]
                params.h        (2,2,:) double = zeros(2, 2, 1)
                params.tau_c    (1,1)   double = 0
                params.tau_s    (1,1)   double = 0
                params.tau_d    (1,2)   double = [0 0]
                params.sigma_t  (1,1)   double = 0
            end

            obj.t     = t;
            obj.dt    = t(2) - t(1);
            obj.m_max = length(t);

            obj.K       = params.K;
            obj.Vb      = params.Vb;
            obj.h       = params.h;
            obj.tau_c   = params.tau_c;
            obj.tau_s   = params.tau_s;
            obj.tau_d   = params.tau_d;
            obj.sigma_t = params.sigma_t;
        end
    end

    methods(Access = protected)
        function setup(obj)

            obj.h     = obj.h * obj.dt;
            obj.Y     = [obj.h(1,1,1) 0; 0 obj.h(2,2,1)];
            obj.H     = [1 obj.h(1,2,1); obj.h(2,1,1) 1];
            
            obj.Fsm   = obj.compute_Fsm(obj.t, obj.K, obj.sigma_t, obj.tau_c, obj.tau_d);
            obj.alpha = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function reset(obj)
            obj.Y      = zeros(2, 2);
            obj.H      = zeros(2, 2);
            obj.Fsm    = zeros(2, obj.m_max);
            obj.alpha  = 0;
        end

        function validate_properties(obj)
            % Validate related or interdependent property values

            mustBeValidTimeStepScalar(obj.K);
            mustBeValidTimeStepScalar(obj.tau_c);
            mustBeValidTimeStepScalar(obj.tau_s);
            mustBeValidTimeStepScalar(obj.sigma_t);

            mustBeValidTimeStepMatrix(obj.Vb,    [2 1]);
            mustBeValidTimeStepMatrix(obj.h,     [2 2 obj.m_max]);
            mustBeValidTimeStepMatrix(obj.tau_d, [1 2]);
        end
    end

    methods
        function [v, i] = simulate(obj)
            obj.validate_properties();
            obj.reset();
            obj.setup();

            v = NaN(2, obj.m_max);
            i = v;

            % Initial conditions
            Fsm_m = diag(obj.Fsm(:, 1));
            v(:, 1) = (obj.Y + obj.H * Fsm_m) \ (obj.H * Fsm_m * obj.Vb);
            i(:, 1) = Fsm_m * (obj.Vb - v(:, 1));

            for m = 2 : obj.m_max
                Fsm_m = diag(obj.Fsm(:, m));
                C     = obj.compute_C(i(:, 1 : m - 1), v(:, 1 : m - 1), obj.h(:, :, 2 : m));

                v(:, m) = obj.compute_vm(obj.alpha, obj.Y, obj.H, Fsm_m, i(:, m - 1), obj.Vb, C);
                i(:, m) = obj.compute_im(obj.alpha, i(:, m - 1), obj.Vb, v(:, m), Fsm_m);
            end
        end
    end

    methods (Static)
        function Fsm = compute_Fsm(t, K, sigma_t, tau_c, tau_d)
            %COMPUTE_FSM Summary of this method goes here
            %   Detailed explanation goes here

            m_max = length(t);
            dt    = t(2) - t(1);

            tau_d1 = tau_d(1);
            tau_d2 = tau_d(2);

            Fsm_1  = NaN(1, m_max);
            Fsm_2  = Fsm_1;

            parfor m = 1 : m_max
                pcm_response = exp(- (t(m) - t(1 : m)) / tau_c);

                Fsm_1(m) = sum(exp(- 0.5 * ( (t(1 : m) - tau_d1) / sigma_t ) .^ 2) .* pcm_response);
                Fsm_2(m) = sum(exp(- 0.5 * ( (t(1 : m) - tau_d2) / sigma_t ) .^ 2) .* pcm_response);
            end

            Fsm = (dt ^ 2) * K * [Fsm_1; Fsm_2];
        end

        function i_impr = compute_i_impr(t, Vb, K, sigma_t, tau_c, tau_d, tau_s)
            %COMPUTE_I_IMPR Summary of this method goes here
            %   Detailed explanation goes here

            m_max  = length(t);
            dt     = t(2) - t(1);

            tau_d1 = tau_d(1);
            tau_d2 = tau_d(2);

            i_impr_1 = NaN(1, m_max);
            i_impr_2 = i_impr_1;

            parfor m = 1 : m_max
                pcm_response  = exp(- (t(m) - t(1 : m)) / tau_c);
                scat_response = (1 - exp(- (t(m) - t(1 : m)) / tau_s));

                i_impr_1(m) = sum(exp(- 0.5 * ( (t(1 : m) - tau_d1) / sigma_t ) ^ 2) ...
                                  .* pcm_response .* scat_response);
                i_impr_2(m) = sum(exp(- 0.5 * ( (t(1 : m) - tau_d2) / sigma_t ) ^ 2) ...
                                  .* pcm_response .* scat_response);
            end

            i_impr = dt * tau_s * Vb * K * [i_impr_1; i_impr_2];
        end

        function vm = compute_vm(alpha, Y, H, Fsm_m, ip, Vb, C)
            %COMPUTE_VM Summary of this method goes here
            %   Detailed explanation goes here

            vm = (Y + H * Fsm_m) \ (alpha * H * ip + H * Fsm_m * Vb + C);
        end

        function im = compute_im(alpha, ip, Vb, vm, Fsm_m)
            %COMPUTE_IM Summary of this method goes here
            %   Detailed explanation goes here
            
            im = alpha * ip + Fsm_m * (Vb - vm);
        end

        function C = compute_C(i, v, h)
            %COMPUTE_C Summary of this method goes here
            %   Detailed explanation goes here

            Ci    = NaN(2, 1);
            Ci(1) = sum(i(1, :) .* fliplr(permute(h(2, 1, :), [1 3 2])));
            Ci(2) = sum(i(2, :) .* fliplr(permute(h(1, 2, :), [1 3 2])));

            Cv    = NaN(2, 1);
            Cv(1) = sum(v(1, :) .* fliplr(permute(h(1, 1, :), [1 3 2])));
            Cv(2) = sum(v(2, :) .* fliplr(permute(h(2, 2, :), [1 3 2])));

            C = Ci - Cv;
        end
    end
end
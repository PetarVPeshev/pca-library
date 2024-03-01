classdef TimeStepCoupling1 < handle
    % TimeStepAlgorithmCoupling Add summary here
    %   Detailed explanation goes here

    properties (SetAccess = immutable)
        t           (1,:) double
        dt          (1,1) double
        m_max       (1,1) double
    end

    properties
        K           (1,1)   double
        Vb          (1,1)   double
        w           (1,:)   double
        h           (2,2,:) double
        tau_c       (1,1)   double
        tau_s       (1,1)   double
        tau_d       (1,2)   double
        sigma_t     (1,1)   double
    end

    properties (SetAccess = protected)
        W0          (2,2) double
        H0          (2,2) double
        Fs          (2,:) double
        i_impr      (1,:) double
        C_impr      (2,2,:) double
        alpha       (1,1) double
    end

    methods
        function obj = TimeStepCoupling1(t, params)
            % Support name-value pair arguments when constructing object
            %   Detailed explanation goes here

            arguments
                t               (1,:)   double
                params.K        (1,1)   double = 0
                params.Vb       (1,1)   double = 0
                params.w        (1,:)   double = double.empty(1, 0)
                params.h        (2,2,:) double = double.empty(2, 2, 0)
                params.tau_c    (1,1)   double = 0
                params.tau_s    (1,1)   double = 0
                params.tau_d    (1,2)   double = zeros(1, 2)
                params.sigma_t  (1,1)   double = 0
            end

            obj.t     = t;
            obj.dt    = t(2) - t(1);
            obj.m_max = length(t);

            obj.K       = params.K;
            obj.Vb      = params.Vb;
            obj.h       = params.h;
            obj.w       = params.w;
            obj.tau_c   = params.tau_c;
            obj.tau_s   = params.tau_s;
            obj.tau_d   = params.tau_d;
            obj.sigma_t = params.sigma_t;
        end
    end

    methods(Access = protected)
        function setup(obj)

            obj.W0     = obj.w(1) * eye(2);
            obj.H0     = obj.h(:, :, 1);
            
            obj.Fs     = obj.compute_Fs(obj.t, obj.K, obj.sigma_t, obj.tau_c, obj.tau_d);
            obj.i_impr = obj.compute_i_impr(obj.t, obj.Vb, obj.K, obj.sigma_t, obj.tau_c, obj.tau_d(1), obj.tau_s);
            obj.C_impr = obj.compute_C_impr(obj.i_impr, obj.h);
            obj.alpha  = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function reset(obj)

            obj.W0     = zeros(2, 2);
            obj.H0     = zeros(2, 2);

            obj.Fs     = zeros(2, obj.m_max);
            obj.i_impr = zeros(1, obj.m_max);
            obj.C_impr = zeros(2, 2, obj.m_max);
            obj.alpha  = 0;
        end

        function validate_properties(obj)
            % Validate related or interdependent property values

            mustBeValidTimeStepScalar(obj.K);
            mustBeValidTimeStepScalar(obj.Vb);
            mustBeValidTimeStepScalar(obj.tau_c);
            mustBeValidTimeStepScalar(obj.tau_s);
            mustBeValidTimeStepScalar(obj.sigma_t);

            mustBeValidTimeStepMatrix(obj.w,     [1 obj.m_max]);
            mustBeValidTimeStepMatrix(obj.h,     [2 2 obj.m_max]);
            mustBeValidTimeStepMatrix(obj.tau_d, [1 2]);
        end
    end

    methods
        function [v, i] = simulate(obj)
            obj.validate_properties();
            obj.reset();
            obj.setup();

            v     = NaN(2, obj.m_max);
            i_int = v;

            % Initial conditions
            Fsm         = diag(obj.Fs(:, 1));
            v(:, 1)     = (obj.W0 + obj.H0 * Fsm) \ (obj.C_impr(:, :, 1) * [1; 1]);
            i_int(:, 1) = Fsm * v(:, 1);

            for m = 2 : obj.m_max
                Fsm     = diag(obj.Fs(:, m));
                Cm_int  = obj.compute_Cm_int(i_int(:, 1 : m - 1), obj.h(:, :, 2 : m));
                Cm_v    = obj.compute_Cm_v(v(:, 1 : m - 1), obj.w(1, 2 : m));

                v(:, m)     = obj.compute_vm(obj.alpha, obj.W0, obj.H0, Fsm, i_int(:, m - 1), ...
                                             obj.C_impr(:, :, m), Cm_int, Cm_v);
                i_int(:, m) = obj.compute_im_int(obj.alpha, i_int(:, m - 1), v(:, m), Fsm);
            end

            i = obj.i_impr - i_int;
        end
    end

    methods (Static)
        function Fs = compute_Fs(t, K, sigma_t, tau_c, tau_d)
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

            Fs = (dt ^ 2) * K * [Fsm_1; Fsm_2];
        end

        function i_impr = compute_i_impr(t, Vb, K, sigma_t, tau_c, tau_d, tau_s)
            %COMPUTE_I_IMPR Summary of this method goes here
            %   Detailed explanation goes here

            m_max  = length(t);
            dt     = t(2) - t(1);

            i_impr = NaN(1, m_max);

            parfor m = 1 : m_max
                pcm_response  = exp(- (t(m) - t(1 : m)) / tau_c);
                scat_response = (1 - exp(- (t(m) - t(1 : m)) / tau_s));

                i_impr(m) = sum(exp(- 0.5 * ( (t(1 : m) - tau_d) / sigma_t ) .^ 2) ...
                                  .* pcm_response .* scat_response);
            end

            i_impr = dt * tau_s * K * Vb * i_impr;
        end

        function C_impr = compute_C_impr(i_impr, h)
            %COMPUTE_C_IMPR Summary of this method goes here
            %   Detailed explanation goes here

            m_max  = size(i_impr, 2);
            C_impr = NaN(2, m_max);

            h      = squeeze(h(:, 1, :));
            i_impr = repmat(i_impr, [2 1]);

            parfor m = 1 : m_max
                C_impr(:, m) = sum(i_impr(:, 1 : m) .* fliplr(h(:, 1 : m)), 2);
            end
        end

        function Cm_int = compute_Cm_int(i_int, h)
            %COMPUTE_CM_INT Summary of this method goes here
            %   Detailed explanation goes here
            % TODO: left @ this point in modifying the algorithm

            h = squeeze(h(:, 1, :));
            i_int = repmat(i_int, );

            h     = permute(h, [2 3 1]);
            i_int = repmat(i_int, [1 1 2]);

            Cm_int = permute(sum(i_int .* fliplr(h), 2), [3 1 2]);
        end

        function Cm_v = compute_Cm_v(v, w)
            %COMPUTE_CM_V Summary of this method goes here
            %   Detailed explanation goes here

            w    = repmat(w, [2 1]);
            Cm_v = sum(v .* fliplr(w), 2);
        end

        function vm = compute_vm(alpha, W0, H0, Fsm, ip_int, Cm_impr, Cm_int, Cm_v)
            %COMPUTE_VM Summary of this method goes here
            %   Detailed explanation goes here

            vm = (W0 + H0 * Fsm) \ ((Cm_impr - Cm_int) * [1; 1] - alpha * H0 * ip_int - Cm_v);
        end

        function im_int = compute_im_int(alpha, ip, vm, Fsm)
            %COMPUTE_IM Summary of this method goes here
            %   Detailed explanation goes here
            
            im_int = alpha * ip + Fsm * vm;
        end
    end
end
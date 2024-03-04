classdef TimeStepCoupling2 < handle
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
        alpha       (1,1) double
    end

    methods
        function obj = TimeStepCoupling2(t, params)
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
            obj.alpha  = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function reset(obj)

            obj.W0     = zeros(2, 2);
            obj.H0     = zeros(2, 2);

            obj.Fs     = zeros(2, obj.m_max);
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

            v = NaN(2, obj.m_max);
            i = v;

            % Initial conditions
            Fsm     = diag(obj.Fs(:, 1));
            v(:, 1) = obj.compute_vm(obj.alpha, obj.Vb, obj.W0, obj.H0, Fsm, [0; 0], [0 0; 0 0], [0; 0]);
            i(:, 1) = obj.compute_im(obj.alpha, obj.Vb, [0; 0], v(:, 1), Fsm);

            for m = 2 : obj.m_max
                Fsm = diag(obj.Fs(:, m));

                Cm_i = obj.compute_Cm_i(i(:, 1 : m - 1), obj.h(:, :, 2 : m));
                Cm_v = obj.compute_Cm_v(v(:, 1 : m - 1), obj.w(:, 2 : m));  % REMOVE ROW SPECIFIER

                v(:, m) = obj.compute_vm(obj.alpha, obj.Vb, obj.W0, obj.H0, Fsm, i(:, m - 1), Cm_i, Cm_v);
                i(:, m) = obj.compute_im(obj.alpha, obj.Vb, i(:, m - 1), v(:, m), Fsm);
            end
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

        function Cm_i = compute_Cm_i(i, h)
            %COMPUTE_CM_I Summary of this method goes here
            %   Detailed explanation goes here

            h = permute(h, [2 3 1]);
            i = repmat(i, [1 1 2]);

            Cm_i = permute(sum(i .* fliplr(h), 2), [3 1 2]);
        end

        function Cm_v = compute_Cm_v(v, w)
            %COMPUTE_CM_V Summary of this method goes here
            %   Detailed explanation goes here

            w    = repmat(w, [2 1]);

            Cm_v = sum(v .* fliplr(w), 2);
        end

        function vm = compute_vm(alpha, Vb, W0, H0, Fsm, ip, Cm_i, Cm_v)
            %COMPUTE_VM Summary of this method goes here
            %   Detailed explanation goes here

            T  = [-1 0; 0 1];
            R  = [1; 1];
            Vb = [Vb; 0];

            vm = (W0 - H0 * Fsm * T) \ (alpha * H0 * ip + H0 * Fsm * Vb + Cm_i * R - Cm_v);
        end

        function im = compute_im(alpha, Vb, ip, vm, Fsm)
            %COMPUTE_IM Summary of this method goes here
            %   Detailed explanation goes here

            T  = [-1 0; 0 1];
            Vb = [Vb; 0];
            
            im = alpha * ip + Fsm * Vb + Fsm * T * vm;
        end
    end
end
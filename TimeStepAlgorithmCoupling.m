classdef TimeStepAlgorithmCoupling < matlab.System
    % TimeStepAlgorithmCoupling Add summary here
    %   Detailed explanation goes here

    properties (SetAccess = immutable)
        t           (1,:) double
        dt          (1,1) double
        m_max       (1,1) double
    end

    properties (Nontunable)
        K           (1,1)   double
        Vb          (2,1)   double
        h           (2,2,:) double
        tau_c       (1,1)   double
        tau_s       (1,1)   double
        tau_d       (1,2)   double
        sigma_t     (1,1)   double
    end

    properties (DiscreteState)
        m
        Ci
        Cv
        ip
        i
        v
        is_done
    end

    properties (SetAccess = protected)
        Y           (2,2) double
        H           (2,2) double
        Fsm         (2,:) double
        alpha       (1,1) double
    end

    methods
        function obj = TimeStepAlgorithmCoupling(t, params)
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
        function setupImpl(obj)

            obj.Y     = [obj.h(1,1,1) * obj.dt,               0; 
                                   0,             obj.h(2,2,1) * obj.dt];

            obj.H     = [          1             obj.h(1,2,1) * obj.dt; 
                         obj.h(2,1,1) * obj.dt              1          ];
            
            obj.Fsm   = obj.compute_Fsm(obj.t, obj.K, obj.sigma_t, obj.tau_c, obj.tau_d);
            obj.alpha = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function [vm, im] = stepImpl(obj)
%             % Current pcm response
%             Fsm_m = diag(obj.Fsm(obj.m));
% 
%             % Current transient voltage step
%             vm = obj.compute_vm(obj.dt, obj.alpha, obj.Y, obj.H, Fsm_m, obj.ip, obj.Vb, obj.Ci, obj.Cv);
% 
%             % Current transient current step
%             im = obj.compute_im(obj.alpha, obj.ip, obj.Vb, vm, Fsm_m);
% 
%             % Update discrete states
%             obj.v(:, obj.m) = vm;
%             obj.i(:, obj.m) = im;
% 
%             if ~obj.is_done
%                 % Update discrete states
%                 obj.m = obj.m + 1;
%     
%                 % Next Ci and Cv
%                 ind = 1 : obj.m - 1;
%                 [obj.Ci, obj.Cv] = obj.compute_C(obj.dt, obj.i(:, ind), obj.v(:, ind), obj.h(:, :, ind + 1));
%     
%                 % Flag
%                 obj.is_done = obj.m == obj.m_max;
%             end

            % Update discrete state time index
            mm = obj.m + 1;

            if ~obj.is_done
                Fsm_m = diag(obj.Fsm(:, mm));
                
                % Compute transient voltage at m
                vm           = obj.compute_vm(obj.alpha, obj.Y, obj.H, Fsm_m, obj.ip, obj.Vb, obj.Ci, obj.Cv);
                obj.v(:, mm) = vm;

                % Compute radiating current at m
                im           = obj.compute_im(obj.alpha, obj.ip, obj.Vb, vm, Fsm_m);
                obj.i(:, mm) = im;

                % Set is_done flag
                obj.is_done = mm == obj.m_max;

                % Compute next step parameters
                if ~obj.is_done
                    % Next Ci and Cv
                    ind = 1 : mm;
                    [obj.Ci, obj.Cv] = obj.compute_C(obj.dt, obj.i(:, ind), obj.v(:, ind), obj.h(:, :, ind + 1));
                end

                obj.m = mm;
            else
                error('TimeStepWeight:isDone:true', 'Time step algorithm finished in last step.');
            end
        end

        function resetImpl(obj)

            obj.i       = NaN(2, obj.m_max);
            obj.v       = NaN(2, obj.m_max);
            obj.is_done = false;
            
%             obj.m       = 1;
            obj.m       = 0;
            obj.Ci      = zeros(2, 1);
            obj.Cv      = zeros(2, 1);
            obj.ip      = zeros(2, 1);
        end

        function releaseImpl(obj)

            obj.Y      = zeros(2, 2);
            obj.H      = zeros(2, 2);
            obj.Fsm    = double.empty(1, 0);
            obj.alpha  = 0;
        end

        function validatePropertiesImpl(obj)
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

        function vm = compute_vm(alpha, Y, H, Fsm_m, ip, Vb, Ci, Cv)
            %COMPUTE_VM Summary of this method goes here
            %   Detailed explanation goes here

            vm = (Y + H * Fsm_m) \ (alpha * H * ip + (H * Fsm_m) * Vb + Ci - Cv);
        end

        function im = compute_im(alpha, ip, Vb, vm, Fsm_m)
            %COMPUTE_IM Summary of this method goes here
            %   Detailed explanation goes here
            
            im = alpha * ip + Fsm_m * (Vb - vm);
        end

        function [Ci, Cv] = compute_C(dt, i, v, h)
            %COMPUTE_C Summary of this method goes here
            %   Detailed explanation goes here

            Ci    = NaN(2, 1);
            Ci(1) = dt * sum(i(1, :) .* fliplr(permute(h(2, 1, :), [1 3 2])));
            Ci(2) = dt * sum(i(2, :) .* fliplr(permute(h(1, 2, :), [1 3 2])));

            Cv    = NaN(2, 1);
            Cv(1) = dt * sum(v(1, :) .* fliplr(permute(h(1, 1, :), [1 3 2])));
            Cv(2) = dt * sum(v(2, :) .* fliplr(permute(h(2, 2, :), [1 3 2])));
        end
    end
end
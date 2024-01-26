classdef TimeStepAlgorithm < matlab.System
    % TimeStep Add summary here
    % Algorithm is defined as closely related to simulation time, for
    % different simulation time, new object must be constructed

    properties (SetAccess = immutable)
        t_vec       (1,:) double
        dt          (1,1) double
        m_max       (1,1) double
        t0_idx      (1,:) double
    end

    properties (Nontunable)
        K           (1,1) double
        Vb          (1,1) double
        za          (1,:) double
        tau_c       (1,1) double
        tau_s       (1,1) double
        sigma_t     (1,1) double
    end

    properties (DiscreteState)
        m
        i_int
        is_done
    end

    properties (SetAccess = protected)
        Fsm         (1,:) double
        i_impr      (1,:) double
        v_impr      (1,:) double
        alpha       (1,1) double  % attenuation of current for one time step
    end

    methods
        function obj = TimeStepAlgorithm(t_vec, params)
            % Support name-value pair arguments when constructing object
            %   t_vec [t] Simulation time [s]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna impedance [Ohm]
            %   tau_c [\tau_{c}] Recombination time [s]
            %   tau_s [\tau_{s}] Scattering time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]

            arguments
                t_vec           (1,:) double
                params.K        (1,1) double = 0
                params.Vb       (1,1) double = 0
                params.za       (1,:) double = 0
                params.tau_c    (1,1) double = 0
                params.tau_s    (1,1) double = 0
                params.sigma_t  (1,1) double = 0
            end

            obj.dt    = t_vec(2) - t_vec(1);
            obj.t_vec = t_vec;
            obj.m_max = length(t_vec);
            obj.t0_idx = find(t_vec == 0);
            if isempty(obj.t0_idx)
                error('TimeStep:t_vec:noSampleAtT0', 'Time vector does not contain sample at t = 0 s.');
            elseif ~isequal(size(obj.t0_idx), [1 1])
                error('TimeStep:t_vec:manySamplesAtT0', 'Time vector contains several samples at t = 0 s.');
            end

            obj.K = params.K;
            obj.Vb = params.Vb;
            obj.za = params.za;
            obj.tau_c = params.tau_c;
            obj.tau_s = params.tau_s;
            obj.sigma_t = params.sigma_t;
        end
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculation of impressed current (can be computed before time step run-time)

            hm_n = obj.compute_hm_n(obj.t_vec, obj.m_max, obj.tau_c, obj.sigma_t);
            obj.Fsm = obj.compute_Fsm(obj.dt, hm_n, obj.K);
            obj.i_impr = obj.compute_i_impr(obj.t_vec, obj.m_max, hm_n, obj.Vb, obj.K, obj.tau_s);
            if length(obj.za) == 1
                obj.v_impr = obj.i_impr * obj.za;
            else
                obj.v_impr = conv(obj.i_impr, obj.za, 'same') * obj.dt;
%                 obj.v_impr = zeros(1, obj.m_max);
%                 for k = 1 : 1 : obj.m_max
%                     for j = 1 : 1 : obj.m_max
%                         obj.v_impr(k) = obj.v_impr(k) + obj.i_impr(j) * obj.za(k - j + 1); 
%                     end
%                 end
%                 obj.v_impr = obj.v_impr * obj.dt;
            end
            obj.alpha = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function [vm, vgm, im, im_int] = stepImpl(obj)
            % Calculate transient voltage, radiating and internal current at time step m as a function of 
            % the antenna parameters and discrete states.

            % Update discrete state time index
            obj.m = obj.m + 1;

            if ~obj.is_done
                % Compute transient voltage at time instance m
                vm = obj.compute_v_step(obj.m, obj.t0_idx, obj.dt, obj.v_impr, obj.i_int, obj.Vb, obj.za, ...
                                        obj.Fsm, obj.alpha);

                % Compute gap voltage at time instance m
                vgm = obj.Vb - vm;

                % Compute internal current at time instance m
                im_int = obj.compute_iint_step(obj.m, obj.i_int, vm, obj.Fsm, obj.alpha);
                obj.i_int(obj.m) = im_int;

                % Compute radiating current at time instance m
                im = obj.i_impr(obj.m) - im_int;

                % Set is_done flag to true
                obj.is_done = obj.m == obj.m_max;
            else
                error('TimeStep:isDone:true', 'Time step algorithm finished in last step.');
            end
        end

        function resetImpl(obj)
            % Initialize / reset time idex to 0, and initial internal current to 0 A

            obj.is_done = false;
            obj.m = 1;
            obj.i_int = NaN(1, obj.m_max);
            obj.i_int(1) = 0;
        end

        function releaseImpl(obj)
            % Release impressed current (set to double empty row vector) and impulse response (set to empty 
            % matrix )

            obj.i_impr = double.empty(1, 0);
            obj.v_impr = double.empty(1, 0);
            obj.Fsm = double.empty(1, 0);
            obj.alpha = 0;
        end

        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values

            obj.check_general_property(obj.K, 'K');
            obj.check_general_property(obj.Vb, 'Vb');
            obj.check_general_property(obj.tau_c, 'tau_c');
            obj.check_general_property(obj.tau_s, 'tau_s');
            obj.check_general_property(obj.sigma_t, 'sigma_t');

            if length(obj.za) ~= obj.m_max && length(obj.za) ~= 1
                error('TimeStep:ga:invalidSize', 'Property ga must be equal in size to t_vec, or be scalar.');
            end

            if any(isnan(obj.za))
                error('TimeStep:ga:NaN', 'Property ga contains NaN values.');
            end
        end

        function check_general_property(~, property, property_name)
            %CHECK_GENERAL_PROPERTY Summary of this method goes here

            if ~(property > 0)
                error(['TimeStep:' property_name ':notPositive'], ...
                    ['Property ' property_name ' is not positive.']);
            end

            if ~isreal(property)
                error(['TimeStep:' property_name ':notReal'], ...
                    ['Property ' property_name ' is not real.']);
            end

            if ~isnumeric(property)
                error(['TimeStep:' property_name ':notNumeric'], ...
                    ['Property ' property_name ' is not numeric.']);
            end

        end
    end

    methods (Static)
        function hm_n = compute_hm_n(t_vec, m_max, tau_c, sigma_t)
            %COMPUTE_HM_N Summary of this method goes here
            %   t_vec [t] Simulation time [s]
            %   m_max [m_{max}] Maximum step
            %   tau_c [\tau_{c}] Recombination time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]

            [N, M] = meshgrid(1 : 1 : m_max, 1 : 1 : m_max);
            
            hm_n = exp(- 0.5 * (t_vec(N) / sigma_t) .^ 2) .* exp(- (t_vec(M) - t_vec(N)) / tau_c);
            hm_n = tril(hm_n);
        end

        function Fsm = compute_Fsm(dt, hm_n, K)
            %COMPUTE_FSM Summary of this method goes here
            %   dt [\delta_{t}] Time step [s]
            %   hm_n [h_{m}[n]] Impulse response (m - rows, n - cols)
            %   K [K] Photo-conductive antenna constant [S/s^2]

            Fsm_n = (dt ^ 2) * K * hm_n;
            Fsm = sum(Fsm_n, 2)';
        end

        function i_impr = compute_i_impr(t_vec, m_max, hm_n, Vb, K, tau_s)
            %COMPUTE_I_IMPR Summary of this method goes here
            %   t_vec [t] Simulation time [s]
            %   m_max [m_{max}] Maximum step
            %   hm_n [h_{m}[n]] Impulse response (m - rows, n - cols)
            %   Vb [V_{b}] Bias voltage [V]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   tau_s [\tau_{s}] Scattering time [s]

            if ~isequal([m_max m_max], size(hm_n))
                error('TimeStep:compute_i_impr:notEqualSize', ...
                    'Argument hm_n has incompatible size to t_vec.');
            end

            dt = t_vec(2) - t_vec(1);
            [N, M] = meshgrid(1 : 1 : m_max, 1 : 1 : m_max);
        
            i_impr = dt * Vb * tau_s * K * sum(hm_n .* (1 - exp(- (t_vec(M) - t_vec(N)) / tau_s)), 2)';
        end

        function vm = compute_v_step(m, t0_idx, dt, v_impr, i_int, Vb, za, Fsm, alpha)
            %COMPUTE_V_STEP Summary of this method goes here
            %   m [m] Step index [-]
            %   dt [\delta_{t}] Time step [s]
            %   v [v[m]] Transient voltage v[m] [V]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna impedance [S]
            %   delta_atten [d\alpha] Attenuation in previous current for one time step
            %   broadband [-] Logical parameter specifying if the antenna is considered broadband [-]

            if length(za) == 1
                vm = (v_impr(m) - za * alpha * i_int(m - 1)) / (1 + za * Fsm(m));
            else
                conv_part = sum( i_int(1 : m - 1) .* fliplr(za(t0_idx + (1 : m - 1))) ) * dt;
                conv_part_2 = 0;
                for k = 1 : 1 : m - 1
                    conv_part_2 = conv_part_2 + i_int(k) * za(t0_idx + m - k);
                end
                conv_part_2 = conv_part_2 * dt;
                if conv_part ~= conv_part_2
                    error('TimeStep:compute_v_step:conv_part', ...
                        ['Convolution part is not equal at m = ' num2str(m)]);
                end

                nom = v_impr(m) - conv_part - dt * za(t0_idx) * alpha * i_int(m - 1);
                den = 1 + dt * za(t0_idx) * Fsm(m);
                vm = nom / den;
            end

            if vm > Vb
                vm = Vb;
            end
        end

        function im_int = compute_iint_step(m, i_int, vm, Fsm, alpha)
            %COMPUTE_IINT_STEP Summary of this method goes here
            %   m [m] Sample index [-]
            %   i_int [i_{int}(m)] Internal current [A]
            %   vm [v(m)] Transient voltage at time instance m [V]
            %   Fsm [F^{s}_{m}] Response at time instance m
            %   alpha [\alpha] Attenuation in previous current for one time step [-]

            im_int = alpha * i_int(m - 1) + vm * Fsm(m);
        end

        function im = compute_i_step(m, i, vgm, Fsm, alpha)
            %COMPUTE_I_STEP Summary of this method goes here
            %COMPUTE_IINT_STEP Summary of this method goes here
            %   m [m] Sample index [-]
            %   i [i(m)] Transient current [A]
            %   vm [v(m)] Transient voltage at time instance m [V]
            %   Fsm [F^{s}_{m}] Response at time instance m
            %   alpha [\alpha] Attenuation in previous current for one time step [-]
            
            im = alpha * i(m - 1) + vgm * Fsm(m);
        end
    end
end

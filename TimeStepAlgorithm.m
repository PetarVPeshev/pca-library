classdef TimeStepAlgorithm < matlab.System
    % TimeStep Add summary here
    % Algorithm is defined as closely related to simulation time, for
    % different simulation time, new object must be constructed

    properties (SetAccess = immutable)
        t_vec       (1,:) double
        dt          (1,1) double
    end

    properties (Hidden, SetAccess = immutable)
        m_max       (1,1) double
    end

    properties (Dependent, Hidden, SetAccess = protected)
        broadband   (1,1) logical
    end

    properties (Nontunable)
        K           (1,1) double
        Vb          (1,1) double
        ga          (1,:) double
        tau_c       (1,1) double
        tau_s       (1,1) double
        sigma_t     (1,1) double
    end

    properties (DiscreteState)
        m
        v
        i
        is_done
    end

    properties (SetAccess = protected)
        Fsm         (1,:) double
        i_impr      (1,:) double
        delta_atten (1,1) double  % attenuation of current for one time step
    end

    methods
        function obj = TimeStepAlgorithm(t_vec, params)
            % Support name-value pair arguments when constructing object
            %   t_vec [t] Simulation time [s]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna impedance [S]
            %   tau_c [\tau_{c}] Recombination time [s]
            %   tau_s [\tau_{s}] Scattering time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]

            arguments
                t_vec           (1,:) double
                params.K        (1,1) double = 0
                params.Vb       (1,1) double = 0
                params.ga       (1,:) double = 0
                params.tau_c    (1,1) double = 0
                params.tau_s    (1,1) double = 0
                params.sigma_t  (1,1) double = 0
            end

            obj.dt    = t_vec(2) - t_vec(1);
            obj.t_vec = t_vec;
            obj.m_max = length(t_vec);

            obj.K = params.K;
            obj.Vb = params.Vb;
            obj.ga = params.ga;
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
            obj.delta_atten = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function [vm, im, im_int] = stepImpl(obj)
            % Calculate transient voltage, radiating and internal current at time step m as a function of 
            % the antenna parameters and discrete states.

            % Update discrete state time index
            obj.m = obj.m + 1;

            if ~obj.is_done
                % Previous time step current
                i_prev = 0;
                if obj.m ~= 1
                    i_prev = obj.i(obj.m - 1);
                end

                % Compute transient voltage at time instance m
                vm = obj.compute_v_step(obj.m, obj.dt, obj.v, i_prev, obj.Fsm(obj.m), obj.Vb, obj.ga, ...
                                        obj.delta_atten, obj.broadband);
                obj.v(obj.m) = vm;

                % Compute radiating current at time instance m
                im = obj.compute_i_step(obj.m, obj.dt, obj.v, obj.ga, obj.broadband);
                obj.i(obj.m) = im;

                % Compute internal current at time instance m
                im_int = obj.i_impr(obj.m) - im;

                % Set is_done flag to true
                obj.is_done = obj.m == obj.m_max;
            else
                error('TimeStep:isDone:true', 'Time step algorithm finished in last step.');
            end
        end

        function resetImpl(obj)
            % Initialize / reset time idex to 0, and initial internal current to 0 A

            obj.m = 0;
            obj.v = NaN(1, obj.m_max);
            obj.i = NaN(1, obj.m_max);
            obj.is_done = false;
        end

        function releaseImpl(obj)
            % Release impressed current (set to double empty row vector) and impulse response (set to empty 
            % matrix )

            obj.i_impr = double.empty(1, 0);
            obj.Fsm = double.empty(1, 0);
            obj.delta_atten = 0;
        end

        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values

            obj.check_general_property(obj.K, 'K');
            obj.check_general_property(obj.Vb, 'Vb');
            obj.check_general_property(obj.tau_c, 'tau_c');
            obj.check_general_property(obj.tau_s, 'tau_s');
            obj.check_general_property(obj.sigma_t, 'sigma_t');

            if length(obj.ga) ~= obj.m_max && length(obj.ga) ~= 1
                error('TimeStep:ga:invalidSize', 'Property ga must be equal in size to t_vec, or be scalar.');
            end

            if any(isnan(obj.ga))
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

    methods
        function broadband = get.broadband(obj)
            %GET.BROADBAND Summary of this method goes here

            broadband = length(obj.ga) == 1;
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

        function vm = compute_v_step(m, dt, v, i_prev, Fsm, Vb, ga, delta_atten, broadband)
            %COMPUTE_V_STEP Summary of this method goes here
            %   m [m] Step index [-]
            %   dt [\delta_{t}] Time step [s]
            %   v [v[m]] Transient voltage v[m] [V]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna impedance [S]
            %   delta_atten [d\alpha] Attenuation in previous current for one time step
            %   broadband [-] Logical parameter specifying if the antenna is considered broadband [-]

%             assert(broadband == true, 'Not implemented non-broadband pca antenna algorithm.');

            if broadband
                vm = (delta_atten * i_prev + Vb * Fsm) / (ga + Fsm);
            else
%                 conv_part = 0;
%                 for n = 1 : 1 : m - 1
%                     conv_part = conv_part + v(n) * ga(m + 1 - n);
%                 end

                if m ~= 1
                    v = v(1 : m - 1);
                    ga_r = fliplr( ga(2 : m) );
                    conv_part = sum(v .* ga_r);
                else
                    conv_part = 0;
                end

                nom = delta_atten * i_prev + Vb * Fsm - dt * conv_part;
                den = (dt * ga(1) + Fsm);
                vm = nom / den;
%                 conv_part = 0;
%                 for n = 1 : 1 : m - 1
%                     conv_part = conv_part + v(n) * ga(m - n + 1);
%                 end

%                 if m == 1
%                     conv_part = 0;
%                 else
%                     conv_part = sum( v(1 : m - 1) .* flip(ga(m : 2)) );
%                 end
% 
%                 nom = dt * delta_atten * i_prev + dt * Vb * hm - conv_part;
%                 den = ga(1) + dt * hm;
%                 vm = nom / den;
            end

%             if vm > Vb
%                 vm = Vb;
%             end
        end

        function im = compute_i_step(m, dt, v, ga, broadband)
            %COMPUTE_I_STEP Summary of this method goes here
            %   dt [\delta_{t}] Time step [s]
            %   vm [v(m)] Transient voltage at time instance m [V]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   hm [h_{m}] Impulse response h_{m} at time instance m
            %   Vb [V_{b}] Bias voltage [V]
            %   delta_atten [d\alpha] Attenuation in previous current for one time step
            
            if broadband
%                 im = delta_atten * i(m - 1) + vg(m) * Fsm;
                im = v(m) * ga;
            else
                im = dt * sum( v(1 : m) .* fliplr(ga(1 : m)) );
            end
        end
    end
end

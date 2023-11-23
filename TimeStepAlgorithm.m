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
        i_prev
        is_done
    end

    properties (SetAccess = protected)
        hm_n        (:,:) double
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

            obj.hm_n = obj.compute_hm_n(obj.t_vec, obj.m_max, obj.K, obj.tau_c, obj.sigma_t);
            obj.i_impr = obj.compute_i_impr(obj.t_vec, obj.m_max, obj.hm_n, obj.Vb, obj.tau_s);
            obj.delta_atten = exp(- obj.dt * (1 / obj.tau_c + 1 / obj.tau_s));
        end

        function [vm, im, im_int] = stepImpl(obj)
            % Calculate transient voltage, radiating and internal current at time step m as a function of 
            % the antenna parameters and discrete states.

            % Update discrete state time index
            obj.m = obj.m + 1;

            if ~obj.is_done
                % Compute transient voltage at time instance m
                hm_m = obj.hm_n(obj.m, obj.m);
                vm = obj.compute_v_step(obj.dt, obj.i_prev, hm_m, obj.Vb, obj.ga, ...
                                        obj.delta_atten, obj.broadband);

                % Compute radiating current at time instance m
                hm = sum(obj.hm_n(obj.m, :));
                im = obj.compute_i_step(obj.dt, vm, obj.i_prev, hm, obj.Vb, obj.delta_atten);

                % Compute internal current at time instance m
                im_int = obj.i_impr(obj.m) - im;

                % Update discrete state previous internal current
                obj.i_prev = im;

                % Set is_done flag to true
                obj.is_done = obj.m == obj.m_max;
            else
                error('TimeStep:isDone:true', 'Time step algorithm finished in last step.');
            end
        end

        function resetImpl(obj)
            % Initialize / reset time idex to 0, and initial internal current to 0 A

            obj.m = 0;
            obj.i_prev = 0;
            obj.is_done = false;
        end

        function releaseImpl(obj)
            % Release impressed current (set to double empty row vector) and impulse response (set to empty 
            % matrix )

            obj.i_impr = double.empty(1, 0);
            obj.h_mn = double.empty(0, 0);
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
        function hm_n = compute_hm_n(t_vec, m_max, K, tau_c, sigma_t)
            %COMPUTE_HM_N Summary of this method goes here
            %   t_vec [t] Simulation time [s]
            %   m_max [m_{max}] Maximum step
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   tau_c [\tau_{c}] Recombination time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]

            [N, M] = meshgrid(1 : 1 : m_max, 1 : 1 : m_max);
            
            hm_n = K * exp(- 0.5 * (t_vec(N) / sigma_t) .^ 2) .* exp(- (t_vec(M) - t_vec(N)) / tau_c);
            hm_n = tril(hm_n);
        end

        function i_impr = compute_i_impr(t_vec, m_max, hm_n, Vb, tau_s)
            %COMPUTE_I_IMPR Summary of this method goes here
            %   t_vec [t] Simulation time [s]
            %   m_max [m_{max}] Maximum step
            %   hm_n [h_{m}[n]] Impulse response (m - rows, n - cols)
            %   Vb [V_{b}] Bias voltage [V]
            %   tau_s [\tau_{s}] Scattering time [s]

            if ~isequal([m_max m_max], size(hm_n))
                error('TimeStep:compute_i_impr:notEqualSize', ...
                    'Argument hm_n has incompatible size to t_vec.');
            end

            dt = t_vec(2) - t_vec(1);
            [N, M] = meshgrid(1 : 1 : m_max, 1 : 1 : m_max);
        
            i_impr = hm_n .* (1 - exp(- (t_vec(M) - t_vec(N)) / tau_s));
            i_impr = Vb * tau_s * dt * sum(i_impr, 2)';
        end

        function vm = compute_v_step(dt, i_prev, hm_m, Vb, ga, delta_atten, broadband)
            %COMPUTE_V_STEP Summary of this method goes here
            %   dt [\delta_{t}] Time step [s]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   hm_m [h_{m}[m]] Impulse response h_{m}[m] at time instance m
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna impedance [S]
            %   delta_atten [d\alpha] Attenuation in previous current for one time step
            %   broadband [-] Logical parameter specifying if the antenna is considered broadband [-]

            assert(broadband == true, 'Not implemented non-broadband pca antenna algorithm.');

            vm = (delta_atten * i_prev + (dt ^ 2) * hm_m * Vb) / (ga + (dt ^ 2) * hm_m);
        end

        function im = compute_i_step(dt, vm, i_prev, hm, Vb, delta_atten)
            %COMPUTE_I_STEP Summary of this method goes here
            %   dt [\delta_{t}] Time step [s]
            %   vm [v(m)] Transient voltage at time instance m [V]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   hm [h_{m}] Impulse response h_{m} at time instance m
            %   Vb [V_{b}] Bias voltage [V]
            %   delta_atten [d\alpha] Attenuation in previous current for one time step
            
            im = delta_atten * i_prev + (Vb - vm) * (dt ^ 2) * hm;
        end
    end
end

classdef TimeStepAlgorithmParallel < matlab.System
    % TimeStep Add summary here
    % Algorithm is defined as closely related to simulation time, for
    % different simulation time, new object must be constructed

    properties (SetAccess = immutable)
        t           (1,:) double
        dt          (1,1) double
        m_max       (1,1) double
    end

    properties (Nontunable)
        K           (1,1) double
        Vb          (1,1) double
        ga          (1,:) double
        w           (1,:) double
        tau_c       (1,1) double
        tau_s       (1,1) double
        sigma_t     (1,1) double
    end

    properties (DiscreteState)
        m
        i_int_prev
        i_int
        v
        is_done
    end

    properties (SetAccess = protected)
        Fsm         (1,:) double
        i_impr      (1,:) double
        v_impr      (1,:) double
        alpha       (1,1) double  % attenuation of current for one time step
    end

    methods
        function obj = TimeStepAlgorithmParallel(t, params)
            % Support name-value pair arguments when constructing object
            %   t [t] Simulation time [s]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna admittance [S]
            %   w [w] Weight function [-]
            %   tau_c [\tau_{c}] Recombination time [s]
            %   tau_s [\tau_{s}] Scattering time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]

            arguments
                t               (1,:) double
                params.K        (1,1) double = 0
                params.Vb       (1,1) double = 0
                params.ga       (1,:) double = 0
                params.w        (1,:) double = 0
                params.tau_c    (1,1) double = 0
                params.tau_s    (1,1) double = 0
                params.sigma_t  (1,1) double = 0
            end

            obj.dt    = t(2) - t(1);
            obj.t     = t;
            obj.m_max = length(t);

            obj.K       = params.K;
            obj.Vb      = params.Vb;
            obj.ga      = params.ga;
            obj.w       = params.w;
            obj.tau_c   = params.tau_c;
            obj.tau_s   = params.tau_s;
            obj.sigma_t = params.sigma_t;
        end
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculation of Fsm, and alpha

            obj.Fsm    = obj.compute_Fsm(obj.t, obj.K, obj.sigma_t, obj.tau_c);
            obj.i_impr = obj.compute_i_impr(obj.t, obj.Vb, obj.K, obj.sigma_t, obj.tau_c, obj.tau_s);
            obj.v_impr = obj.compute_v_impr(obj.i_impr, obj.ga);
            obj.alpha  = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function [vm, vgm, im, im_int] = stepImpl(obj)
            % Calculate transient voltage, radiating and internal current at time step m as a function of 
            % the antenna parameters and discrete states.

            % Update discrete state time index
            obj.m = obj.m + 1;

            if ~obj.is_done
                Fsm_m = obj.Fsm(obj.m);

                % Compute transient voltage at m
                vm = obj.compute_vm(obj.m, obj.v, obj.v_impr(obj.m), obj.i_int, obj.i_int_prev, obj.ga, ...
                                    obj.w, Fsm_m, obj.alpha);
                obj.v(obj.m) = vm;

                % Compute total voltage at m
                vgm = obj.Vb - vm;

                % Compute transient current at m
                im_int = obj.compute_im_int(obj.i_int_prev, vm, Fsm_m, obj.alpha);
                obj.i_int_prev = im_int;
                obj.i_int(obj.m) = im_int;

                % Compute total current at m
                im = obj.i_impr(obj.m) - im_int;

                % Set is_done flag
                obj.is_done = obj.m == obj.m_max;
            else
                error('TimeStepWeight:isDone:true', 'Time step algorithm finished in last step.');
            end
        end

        function resetImpl(obj)
            % Initialize / reset time idex to 0, initial current at 0 A, allocated space for current vector, 
            % and set is_done flag to false

            obj.m           = 0;
            obj.i_int_prev  = 0;
            obj.i_int       = NaN(1, obj.m_max);
            obj.v           = NaN(1, obj.m_max);
            obj.is_done     = false;
        end

        function releaseImpl(obj)
            % Release Fsm, impressed current (set to double empty row vector), and step attenuation (set to 
            % empty matrix)

            obj.Fsm    = double.empty(1, 0);
            obj.i_impr = double.empty(1, 0);
            obj.v_impr = double.empty(1, 0);
            obj.alpha  = 0;
        end

        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values

            obj.check_general_property(obj.K, 'K');
            obj.check_general_property(obj.Vb, 'Vb');
            obj.check_general_property(obj.tau_c, 'tau_c');
            obj.check_general_property(obj.tau_s, 'tau_s');
            obj.check_general_property(obj.sigma_t, 'sigma_t');

            if length(obj.ga) ~= obj.m_max
                error('TimeStepWeight:ga:invalidSize', ['Antenna impulse response admittance must be ' ...
                    'same length as t.']);
            end

            if any(isnan(obj.ga))
                error('TimeStepWeight:ga:NaN', 'Antenna impulse response admittance contains NaN values.');
            end

            if length(obj.w) ~= obj.m_max
                error('TimeStepWeight:w:invalidSize', 'Weight function must be same length as t.');
            end

            if any(isnan(obj.w))
                error('TimeStepWeight:w:NaN', 'Weight function contains NaN values.');
            end
        end

        function check_general_property(~, property, property_name)
            %CHECK_GENERAL_PROPERTY Summary of this method goes here

            if ~(property > 0)
                error(['TimeStepWeight:' property_name ':notPositive'], ...
                    ['Property ' property_name ' is not positive.']);
            end

            if ~isreal(property)
                error(['TimeStepWeight:' property_name ':notReal'], ...
                    ['Property ' property_name ' is not real.']);
            end

            if ~isnumeric(property)
                error(['TimeStepWeight:' property_name ':notNumeric'], ...
                    ['Property ' property_name ' is not numeric.']);
            end

        end
    end

    methods (Static)
        function Fsm = compute_Fsm(t, K, sigma_t, tau_c)
            %COMPUTE_FSM Summary of this method goes here
            %   dt [\delta_{t}] Time step [s]
            %   hm_n [h_{m}[n]] Impulse response (m - rows, n - cols)
            %   K [K] Photo-conductive antenna constant [S/s^2]

            m_max = length(t);
            dt    = t(2) - t(1);
            Fsm   = NaN(1, m_max);

            parfor m = 1 : m_max
                hm = 0;
                for n = 1 : m
                    hm = hm + exp(- 0.5 * (t(n) / sigma_t) ^ 2) * exp(- (t(m) - t(n)) / tau_c);
                end

                Fsm(m) = hm;
            end

            Fsm = (dt ^ 2) * K * Fsm;
        end

        function i_impr = compute_i_impr(t, Vb, K, sigma_t, tau_c, tau_s)
            %COMPUTE_I_IMPR Summary of this method goes here
            %   t [t] Simulation time [s]
            %   hm_n [h_{m}[n]] Impulse response (m - rows, n - cols)
            %   Vb [V_{b}] Bias voltage [V]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   tau_s [\tau_{s}] Scattering time [s]

            m_max = length(t);
            dt     = t(2) - t(1);
            i_impr = NaN(1, m_max);

            parfor m = 1 : m_max
                Hm = 0;
                for n = 1 : m
                    Hm = Hm + exp(- 0.5 * (t(n) / sigma_t) ^ 2) * exp(- (t(m) - t(n)) / tau_c) ...
                            * (1 - exp(- (t(m) - t(n)) / tau_s));
                end

                i_impr(m) = Hm;
            end

            i_impr = dt * tau_s * Vb * K * i_impr;
        end

        function v_impr = compute_v_impr(i_impr, ga)
            %COMPUTE_V_IMPR Summary of this method goes here
            %   i_impr [i_{impr}] Impressed current [A]
            %   ga [g_{a}] Antenna impulse response admittance [S]

            Nga = length(ga);
            if length(i_impr) ~= Nga
                error('TimeStepWeight:compute_v_impr:notEqualSize', ...
                      'Arguments i_impr and ga must be of equal length');
            end

            v_impr = NaN(1, Nga);

            parfor m = 1 : Nga
                v_impr(m) = sum(i_impr(1 : m) .* fliplr(ga(1 : m)));
            end
        end

        function vm = compute_vm(m, v, vm_impr, i_int, i_int_prev, gl, w, Fsm_m, alpha)
            %COMPUTE_V_STEP Summary of this method goes here
            %   m [m] Step index [-]
            %   dt [\delta_{t}] Time step [s]
            %   v [v[m]] Transient voltage v[m] [V]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Antenna impulse response admittance [S]
            %   alpha [\alpha] Attenuation in previous current for one time step [-]

            if m == 1
                vm = (vm_impr - alpha * i_int_prev * gl(1)) / (Fsm_m * gl(1) + w(1));
            else
                conv_1 = sum(v(1 : m - 1) .* fliplr(w(2 : m)));
                conv_2 = sum(i_int(1 : m - 1) .* fliplr(gl(2 : m)));
                vm = (vm_impr - alpha * i_int_prev * gl(1) - conv_1 - conv_2) / (Fsm_m * gl(1) + w(1));
            end
        end

        function im_int = compute_im_int(i_int_prev, vm, Fsm_m, alpha)
            %COMPUTE_I_STEP Summary of this method goes here
            %   i_int_prev [i(m-1)] Internal current at time instance m - 1 [A]
            %   vm [v(m)] Transient voltage at time instance m [V]
            %   Fsm_m [F_{m}^{s}[m]] Impulse response [S]
            %   alpha [\alpha] Attenuation in previous current for one time step [-]
            
            im_int = alpha * i_int_prev + vm * Fsm_m;
        end
    end
end
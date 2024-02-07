classdef TimeStepAlgorithm < matlab.System
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
        tau_c       (1,1) double
        tau_s       (1,1) double
        sigma_t     (1,1) double
    end

    properties (DiscreteState)
        m
        i_prev
        v
        is_done
    end

    properties (SetAccess = protected)
        Fsm         (1,:) double
        i_impr      (1,:) double
        alpha       (1,1) double  % attenuation of current for one time step
    end

    methods
        function obj = TimeStepAlgorithm(t, params)
            % Support name-value pair arguments when constructing object
            %   t [t] Simulation time [s]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna admittance [S]
            %   tau_c [\tau_{c}] Recombination time [s]
            %   tau_s [\tau_{s}] Scattering time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]

            arguments
                t               (1,:) double
                params.K        (1,1) double = 0
                params.Vb       (1,1) double = 0
                params.ga       (1,:) double = 0
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
            obj.tau_c   = params.tau_c;
            obj.tau_s   = params.tau_s;
            obj.sigma_t = params.sigma_t;
        end
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculation of hm_n, Fsm, and alpha

            hm_n       = obj.compute_hm_n(obj.t, obj.tau_c, obj.sigma_t);
            obj.Fsm    = obj.compute_Fsm(obj.dt, hm_n, obj.K);
            obj.i_impr = obj.compute_i_impr(obj.t, hm_n, obj.Vb, obj.K, obj.tau_s);
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
                if length(obj.ga) == 1
                    vm = obj.compute_vm_bband(obj.i_prev, obj.Vb, obj.ga, Fsm_m, obj.alpha);
                else
                    vm = obj.compute_vm(obj.m, obj.dt, obj.v, obj.i_prev, obj.Vb, obj.ga, Fsm_m, obj.alpha);
                end
                obj.v(obj.m) = vm;

                % Compute total voltage at m
                vgm = obj.Vb - vm;

                % Compute gap current at m
                im = obj.compute_im(obj.i_prev, vgm, Fsm_m, obj.alpha);
                if im < 0
                    im = 0;
                end
                if im > obj.i_impr(obj.m)
                    im = obj.i_impr(obj.m);
                end
                obj.i_prev = im;

                % Compute internal current at m
                im_int = obj.i_impr(obj.m) - im;

                % Set is_done flag
                obj.is_done = obj.m == obj.m_max;
            else
                error('TimeStep:isDone:true', 'Time step algorithm finished in last step.');
            end
        end

        function resetImpl(obj)
            % Initialize / reset time idex to 0, initial current at 0 A, allocated space for current vector, 
            % and set is_done flag to false

            obj.m       = 0;
            obj.i_prev  = 0;
            obj.v       = NaN(1, obj.m_max);
            obj.is_done = false;
        end

        function releaseImpl(obj)
            % Release Fsm, impressed current (set to double empty row vector), and step attenuation (set to 
            % empty matrix)

            obj.Fsm    = double.empty(1, 0);
            obj.i_impr = double.empty(1, 0);
            obj.alpha  = 0;
        end

        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values

            obj.check_general_property(obj.K, 'K');
            obj.check_general_property(obj.Vb, 'Vb');
            obj.check_general_property(obj.tau_c, 'tau_c');
            obj.check_general_property(obj.tau_s, 'tau_s');
            obj.check_general_property(obj.sigma_t, 'sigma_t');

            if length(obj.ga) ~= obj.m_max && length(obj.ga) ~= 1
                error('TimeStep:ga:invalidSize', ['Antenna impulse response admittance must be same ' ...
                    'length as t, or scalar.']);
            end

            if any(isnan(obj.ga))
                error('TimeStep:ga:NaN', 'Antenna impulse response admittance contains NaN values.');
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
        function hm_n = compute_hm_n(t, tau_c, sigma_t)
            %COMPUTE_HM_N Summary of this method goes here
            %   t [t] Simulation time [s]
            %   m_max [m_{max}] Maximum step
            %   tau_c [\tau_{c}] Recombination time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]
            
            m_max = length(t);
            [N, M] = meshgrid(1 : m_max, 1 : m_max);
            
            hm_n = exp(- 0.5 * (t(N) / sigma_t) .^ 2) .* exp(- (t(M) - t(N)) / tau_c);
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

        function i_impr = compute_i_impr(t, hm_n, Vb, K, tau_s)
            %COMPUTE_I_IMPR Summary of this method goes here
            %   t [t] Simulation time [s]
            %   hm_n [h_{m}[n]] Impulse response (m - rows, n - cols)
            %   Vb [V_{b}] Bias voltage [V]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   tau_s [\tau_{s}] Scattering time [s]

            m_max = length(t);
            if ~isequal([m_max m_max], size(hm_n))
                error('TimeStep:compute_i_impr:notEqualSize', ...
                    'Argument hm_n has incompatible size to t.');
            end

            dt = t(2) - t(1);
            [N, M] = meshgrid(1 : m_max, 1 : m_max);
        
            bm_n = tril(1 - exp(- (t(M) - t(N)) / tau_s));
            i_impr = dt * Vb * K * tau_s * sum(hm_n .* bm_n, 2)';
        end

        function vm = compute_vm_bband(i_prev, Vb, Gl, Fsm_m, alpha)
            %COMPUTE_V_STEP Summary of this method goes here
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   Vb [V_{b}] Bias voltage [V]
            %   Gl [G_{l}] Radiating element admittance [S]
            %   Fsm_m [F_{m}^{s}[m]] Impulse response [S]
            %   alpha [\alpha] Attenuation in previous current for one time step [-]

            vm = (alpha * i_prev + Vb * Fsm_m) / (Gl + Fsm_m);

            % Saturation condition
            if vm > Vb
                vm = Vb;
            end
        end

        function vm = compute_vm(m, dt, v, i_prev, Vb, ga, Fsm_m, alpha)
            %COMPUTE_V_STEP Summary of this method goes here
            %   m [m] Step index [-]
            %   dt [\delta_{t}] Time step [s]
            %   v [v[m]] Transient voltage v[m] [V]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Antenna impulse response admittance [S]
            %   alpha [\alpha] Attenuation in previous current for one time step [-]

            if m == 1
                vm = (alpha * i_prev + Vb * Fsm_m) / (dt * ga(1) + Fsm_m);
            else
                conv_part = dt * sum(v(1 : m - 1) .* fliplr(ga(2 : m)));
                vm = (alpha * i_prev + Vb * Fsm_m - conv_part) / (dt * ga(1) + Fsm_m);
            end

            % Saturation condition
            if vm > Vb
                vm = Vb;
            end

            % No negative transient voltage
            % Caused by oscillatory behaviour in the admittance due to the
            % numerically computed IFT
            if vm < 0
                vm = 0;
            end
        end

        function im = compute_im(i_prev, vgm, Fsm_m, alpha)
            %COMPUTE_I_STEP Summary of this method goes here
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   vgm [v_{g}(m)] Gap voltage at time instance m [V]
            %   Fsm_m [F_{m}^{s}[m]] Impulse response [S]
            %   alpha [\alpha] Attenuation in previous current for one time step [-]
            
            im = alpha * i_prev + vgm * Fsm_m;
        end
    end
end
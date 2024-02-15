classdef TimeStepAlgorithmTwoFeeds2 < matlab.System
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
        h           (2,2) cell
        w           (1,:) double
        tau_c       (1,1) double
        tau_s       (1,1) double
        sigma_t     (1,1) double
    end

    properties (DiscreteState)
        m
        i_int
        v
        T_prev
        C_int
        Cv
        is_done
    end

    properties (SetAccess = protected)
        w0          (2,2) double
        h0          (2,2) double
        Fsm         (1,:) double
        i_impr      (1,:) double
        C_impr      (2,2) cell
        alpha       (1,1) double  % attenuation of current for one time step
    end

    methods
        function obj = TimeStepAlgorithmTwoFeeds2(t, params)
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
                params.h        (2,2) cell   = cell(2,2)
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
            obj.h       = params.h;
            obj.w       = params.w;
            obj.tau_c   = params.tau_c;
            obj.tau_s   = params.tau_s;
            obj.sigma_t = params.sigma_t;
        end
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculation of Fsm, i_impr, and alpha

            obj.w0     = eye(2) * obj.w(1);
            obj.h0     = [obj.h{1,1}(1) obj.h{1,2}(1); obj.h{2,1}(1) obj.h{2,2}(1)];
            obj.Fsm    = obj.compute_Fsm(obj.t, obj.K, obj.sigma_t, obj.tau_c);
            obj.i_impr = obj.compute_i_impr(obj.t, obj.Vb, obj.K, obj.sigma_t, obj.tau_c, obj.tau_s);
            obj.C_impr = obj.compute_Cimpr(obj.i_impr, obj.h);
            obj.alpha  = exp(- obj.dt / obj.tau_c) * exp(- obj.dt / obj.tau_s);
        end

        function [vm, vgm, im, im_int] = stepImpl(obj)
            % Calculate transient voltage, radiating and internal current at time step m as a function of 
            % the antenna parameters and discrete states.

            % Update discrete state time index
            mm = obj.m + 1;

            if ~obj.is_done
                Fsm_m   = obj.Fsm(mm);
                Cimpr_m = [obj.C_impr{1,1}(mm) obj.C_impr{1,2}(mm); ...
                           obj.C_impr{2,1}(mm) obj.C_impr{2,2}(mm)];
                
                % Compute transient voltage at m
                vm           = obj.compute_vm(obj.w0, Fsm_m, obj.h0, Cimpr_m, obj.T_prev, obj.C_int, obj.Cv);
                obj.v(:, mm) = vm;

                % Compute gap voltage at m
                vgm = obj.Vb - vm;

                % Compute radiating current at m
                im_int           = obj.compute_im_int(vm, obj.T_prev, Fsm_m);
                obj.i_int(:, mm) = im_int;

                % Compute internal current at m
                im = obj.i_impr(mm) - im_int;

                % Set is_done flag
                obj.is_done = mm == obj.m_max;

                % Compute next step parameters
                if ~obj.is_done
                    obj.T_prev = obj.alpha * im_int;

                    Cint_11   = sum(obj.i_int(1, 1 : mm) .* fliplr(obj.h{1,1}(2 : mm + 1)));
                    Cint_12   = sum(obj.i_int(2, 1 : mm) .* fliplr(obj.h{1,2}(2 : mm + 1)));
                    Cint_21   = sum(obj.i_int(1, 1 : mm) .* fliplr(obj.h{2,1}(2 : mm + 1)));
                    Cint_22   = sum(obj.i_int(2, 1 : mm) .* fliplr(obj.h{2,2}(2 : mm + 1)));
                    obj.C_int = [Cint_11 Cint_12; Cint_21 Cint_22];

                    Cv1    = sum(obj.v(1, 1 : mm) .* fliplr(obj.w(2 : mm + 1)));
                    Cv2    = sum(obj.v(2, 1 : mm) .* fliplr(obj.w(2 : mm + 1)));
                    obj.Cv = [Cv1; Cv2];
                end

                obj.m = mm;
            else
                error('TimeStepWeight:isDone:true', 'Time step algorithm finished in last step.');
            end
        end

        function resetImpl(obj)
            % Initialize / reset time idex to 0, initial current at 0 A, allocated space for current vector, 
            % and set is_done flag to false

            obj.m       = 0;
            obj.i_int       = NaN(2, obj.m_max);
            obj.v       = NaN(2, obj.m_max);
            obj.T_prev  = [0; 0];
            obj.C_int   = [0 0; 0 0];
            obj.Cv      = [0; 0];
            obj.is_done = false;
        end

        function releaseImpl(obj)
            % Release Fsm, impressed current (set to double empty row vector), and step attenuation (set to 
            % empty matrix)

            obj.w0     = 0;
            obj.h0     = 0;
            obj.Fsm    = double.empty(1, 0);
            obj.i_impr = double.empty(1, 0);
            obj.C_impr = cell(2, 2);
            obj.alpha  = 0;
        end

        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values

            obj.check_general_property(obj.K, 'K');
            obj.check_general_property(obj.Vb, 'Vb');
            obj.check_general_property(obj.tau_c, 'tau_c');
            obj.check_general_property(obj.tau_s, 'tau_s');
            obj.check_general_property(obj.sigma_t, 'sigma_t');

%             if length(obj.h) ~= obj.m_max
%                 error('TimeStepWeight:ga:invalidSize', ['Antenna impulse response admittance must be ' ...
%                     'same length as t.']);
%             end
% 
%             if any(isnan(obj.h))
%                 error('TimeStepWeight:ga:NaN', 'Antenna impulse response admittance contains NaN values.');
%             end
% 
%             if length(obj.w) ~= obj.m_max
%                 error('TimeStepWeight:w:invalidSize', 'Weight function must be same length as t.');
%             end
% 
%             if any(isnan(obj.w))
%                 error('TimeStepWeight:w:NaN', 'Weight function contains NaN values.');
%             end
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

        function Cimpr = compute_Cimpr(i_impr, h)
            %COMPUTE_CIMPR Summary of this method goes here
            % Detailed explanation goes here
            m_max = length(i_impr);

            h_11 = h{1,1};
            h_12 = h{1,2};
            h_21 = h{2,1};
            h_22 = h{2,2};

            parfor m = 1 : m_max
                Cimpr_11(m) = sum(i_impr(1 : m) .* fliplr(h_11(1 : m)));
                Cimpr_12(m) = sum(i_impr(1 : m) .* fliplr(h_12(1 : m)));
                Cimpr_21(m) = sum(i_impr(1 : m) .* fliplr(h_21(1 : m)));
                Cimpr_22(m) = sum(i_impr(1 : m) .* fliplr(h_22(1 : m)));
            end

            Cimpr = {Cimpr_11, Cimpr_12; Cimpr_21, Cimpr_22};
        end

        function vm = compute_vm(w0, Fsm_m, h0, Cimpr_m, T_prev, C_int, Cv)
            %COMPUTE_V_STEP Summary of this method goes here
            %   m [m] Step index [-]
            %   dt [\delta_{t}] Time step [s]
            %   v [v[m]] Transient voltage v[m] [V]
            %   i_prev [i(m-1)] Radiating current at time instance m - 1 [A]
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Antenna impulse response admittance [S]
            %   alpha [\alpha] Attenuation in previous current for one time step [-]

%             vm = inv(w(1) + Fsm_m * h0) * (h0 * Tm_prev + C_prev * [1; 1] - Cv);
            vm = (w0 + Fsm_m * h0) \ (Cimpr_m * [1; 1] - h0 * T_prev - C_int * [1; 1] - Cv);
        end

        function im_int = compute_im_int(vm, T_prev, Fsm_m)
            %COMPUTE_I_STEP Summary of this method goes here
            %   i_int_prev [i(m-1)] Internal current at time instance m - 1 [A]
            %   vm [v(m)] Transient voltage at time instance m [V]
            %   Fsm_m [F_{m}^{s}[m]] Impulse response [S]
            %   alpha [\alpha] Attenuation in previous current for one time step [-]
            
            im_int = T_prev + vm * Fsm_m;
        end
    end
end
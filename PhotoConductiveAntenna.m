classdef PhotoConductiveAntenna < handle
    %PHOTOCONDUCTIVEANTENNA Summary of this class goes here
    %   Detailed explanation goes here

    properties (Hidden, SetAccess = immutable)
        laser           Laser
        pc_material     PhotoConductor
    end
    
    properties (SetAccess = immutable)
        qo_link         QuasiOpticalLink
    end

    properties (SetAccess = protected)
        time_step
    end

    properties
        Vb              (1,1) double {mustBeReal, mustBeNonNan}
        ga              (1,1) double {mustBeNonNan}
    end

    properties (Dependent, SetAccess = protected)
        P_opt           (1,1) double
        K               (1,1) double
    end
    
    methods
        function obj = PhotoConductiveAntenna(laser, pc_material, Vb, ga, options)
            %PHOTOCONDUCTIVEANTENNA Construct an instance of this class
            %   laser Pulsed laser object
            %   pc_material Photo-conductor object
            %   Vb [V_{b}] Bias voltage [V]
            %   ga [g_{a}] Radiating element impulse response / antenna impedance [S]
            %   t_vec [t] Simulation time [s]
            %   eta_opt [\eta_{opt}] Manually defined optical efficiency (not calculated from parameters) [-]
            %   theta_inc [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   er_inc [\epsilon_{r,inc}] Relative permittivity of medium one (incidence medium) [-]
            
            arguments
                laser               (1,1) Laser
                pc_material         (1,1) PhotoConductor
                Vb                  (1,1) double {mustBePositive, mustBeReal, mustBeNonNan}
                ga                  (1,:) double {mustBeNonNan}
                options.t_vec       (1,2) double {mustBeReal} = NaN
                options.eta_opt     (1,1) double {mustBeInRange(options.eta_opt, 0, 1)} = 0
                options.theta_inc   (1,1) double {mustBeInRange(options.theta_inc, 0, 90)} = 0
                options.er_inc      (1,1) double {mustBeGreaterThanOrEqual(options.er_inc, 1), mustBeReal} = 1
            end
            
            % TODO: Find better solution to adding the utility path in every class.
            addpath([pwd() '\utils']);

            obj.laser = laser;
            obj.pc_material = pc_material;

            obj.Vb = Vb;
            obj.ga = ga;

            obj.qo_link = QuasiOpticalLink(laser, pc_material, eta_opt = options.eta_opt, ...
                                           theta_inc = options.theta_inc, er_inc = options.er_inc);
            
            if isnan(options.t_vec)
                options.t_vec = (-20 : 1e-2 : 20) * laser.tau_p;
            end
            obj.time_step = TimeStepAlgorithm(options.t_vec);
        end

        function P_opt = get.P_opt(obj)
            %GET.PWR_OPT Summary of this method goes here

            P_opt = obj.qo_link.eta_opt * obj.laser.P;
        end
        
        function K = get.K(obj)
            %GET.K Summary of this method goes here
            
            qe = get_phys_const('ElectronCharge');
            h = get_phys_const('PlanckConstant');

            nom = qe * obj.pc_material.mu_dc * obj.laser.T * obj.P_opt;
            den = (obj.pc_material.dimensions.Wx ^ 2) * h * obj.laser.freq * obj.pc_material.tau_s ...
                * sqrt(2 * pi) * obj.laser.sigma_t;
            K = nom / den;
        end

        function change_sim_time(obj, t_vec)
            %CHANGE_SIM_TIME Summary of this method goes here
            %   t_vec [t] New simulation time [s]

            delete(obj.time_step);
            obj.time_step = TimeStepAlgorithm(t_vec);
        end

        function [v, vg, i_impr, i_int, i] = compute_response(obj)
            %COMPUTE_RESPONSE Summary of this method goes here
            
            % Reset time step algorithm
            release(obj.time_step);
            reset(obj.time_step);

            % Update parameters
            obj.time_step.K = obj.K;
            obj.time_step.Vb = obj.Vb;
            obj.time_step.ga = obj.ga;
            obj.time_step.tau_c = obj.pc_material.tau_rec;
            obj.time_step.tau_s = obj.pc_material.tau_s;
            obj.time_step.sigma_t = obj.laser.sigma_t;

            % Compute response
            m_max = length(obj.time_step.t_vec);
            
            v = NaN(1, m_max);
            i = NaN(1, m_max);
            i_int = NaN(1, m_max);
            
            for m = 1 : 1 : m_max
                [v(m), i(m), i_int(m)] = step(obj.time_step);
            end

            vg = obj.Vb - v;
            i_impr = obj.time_step.i_impr;
        end
    end

    methods (Static)
        function i_impr = compute_i_impr_analytical(t_vec, K, Vb, tau_c, tau_s, sigma_t)
            %COMPUTE_I_IMPR_ANALYTICAL Summary of this function goes here
            %   t_vec [t] Simulation time [s]
            %   K [K] Photo-conductive antenna constant [S/s^2]
            %   Vb [V_{b}] Bias voltage [V]
            %   tau_c [\tau_{c}] Recombination time [s]
            %   tau_s [\tau_{s}] Scattering time [s]
            %   sigma_t [\sigma_{t}] Temporal standard deviation of laser [s]

            %   TODO: the analytical solution does produce correct results, check implementation and 
            % mathematical formulation
                i_impr = exp(- t_vec / tau_c) .* exp(0.5 * (sigma_t / tau_c) ^ 2) ...
                    .* ( erf(t_vec / (sqrt(2) * sigma_t) - sigma_t / (sqrt(2) * tau_c)) + 1 - exp(- t_vec / tau_s) ...
                    .* exp((sigma_t ^ 2) / (tau_c * tau_s)) .* exp(0.5 * (sigma_t / tau_s) ^ 2) ...
                    .* (erf(t_vec / (sqrt(2) * sigma_t) - sigma_t  * (tau_s + tau_c) / (sqrt(2) * tau_s * tau_c) + 1)) );
                i_impr = i_impr * K * Vb * tau_s * sigma_t * sqrt(pi / 2);
        end
    end
end

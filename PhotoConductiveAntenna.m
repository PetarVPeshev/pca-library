classdef PhotoConductiveAntenna < handle
    %PHOTOCONDUCTIVEANTENNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = protected, SetAccess = immutable)
        laser           Laser
        pc_material     PhotoConductor
        qo_link         QuasiOpticalLink
        marching_on     TimeStepAlgorithm
    end

    properties (SetAccess = immutable)
        inc_angle       (1,1) double
        inc_rel_permit  (1,1) double
    end

    properties (Dependent, SetAccess = protected)
        pwr_opt (1,1) double
        k_const (1,1) double
    end
    
    methods
        function obj = PhotoConductiveAntenna(laser, photoconductor, marching_on, options)
            %PHOTOCONDUCTIVEANTENNA Construct an instance of this class
            %   laser Pulsed laser object
            %   photoconductor Photo-conductor object
            %   eff_opt [\eta_{opt}] Manually defined optical efficiency (not calculated from parameters) [-]
            %   inc_angle [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   inc_rel_permit [\epsilon_{r}^{inc}] Relative permittivity of medium one (incidence medium) [-]
            
            arguments
                laser                   (1,1) Laser
                photoconductor          (1,1) PhotoConductor
                marching_on             (1,1) TimeStepAlgorithm
                options.eff_opt         (1,1) double {mustBeInRange(options.eff_opt, 0, 1)} = 0
                options.inc_angle       (1,1) double {mustBeInRange(options.inc_angle, 0, 90)} = 0
                options.inc_rel_permit  (1,1) double {mustBeGreaterThanOrEqual(options.inc_rel_permit, 1), ...
                                                     mustBeReal} = 1
            end
            
            % TODO: Find better solution to adding the utility path in every class.
            addpath([pwd() '\utils']);

            obj.laser = laser;
            obj.pc_material = photoconductor;
            obj.marching_on = marching_on;

            obj.qo_link = QuasiOpticalLink(laser, photoconductor, eta_opt = options.eff_opt, ...
                                        theta_inc = options.inc_angle, er_inc = options.inc_rel_permit);

            obj.inc_angle = options.inc_angle;
            obj.inc_rel_permit = options.inc_rel_permit;
        end

        function pwr_opt = get.pwr_opt(obj)
            %GET.PWR_OPT Summary of this method goes here

            pwr_opt = obj.qo_link.eta_opt * obj.laser.P;
        end
        
        function k_const = get.k_const(obj)
            %METHOD1 Summary of this method goes here

            k_const = get_phys_const('ElectronCharge') * obj.pc_material.mu_dc ...
                * obj.laser.T * obj.pwr_opt / ( (obj.pc_material.dimensions.Wx ^ 2) ...
                * get_phys_const('PlanckConstant') * obj.laser.freq * obj.pc_material.tau_s ...
                * sqrt(2 * pi) * obj.laser.sigma_t);
        end
    end
end


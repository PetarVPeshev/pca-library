classdef QuasiOpticalLink < handle
    %QUASIOPTICALLINK Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess = protected, SetAccess = immutable)
        % TODO: Add a default size (1,1) for laser and photoconductor properties, the object constructors must
        % be changed to have zero argument inputs
        laser           Laser
        pc_material     PhotoConductor
        default_eta_opt (1,1) double = NaN
    end

    properties
        theta_inc       (1,1) double
        er_inc          (1,1) double
    end
    
    properties (Dependent, SetAccess = protected)
        eta_tx      (1,1) double
        eta_absorp  (1,1) double
        eta_so      (1,1) double
        eta_opt     (1,1) double
    end
    
    methods
        function obj = QuasiOpticalLink(laser, pc_material, options)
            %EFFICIENCY Construct an instance of this class
            %   laser Pulsed laser
            %   pc_material Photo-conductor
            %   eta_opt [\eta_{opt}] Manually defined optical efficiency (not calculated from parameters) [-]
            %   theta_inc [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   er_inc [\epsilon_{r,inc}] Relative permittivity of medium one (incidence medium) [-]

            arguments
                laser               (1,1) Laser
                pc_material         (1,1) PhotoConductor
                options.eta_opt     (1,1) double {mustBeInRange(options.eta_opt, 0, 1)} = 0
                options.theta_inc   (1,1) double {mustBeInRange(options.theta_inc, 0, 90)} = 0
                options.er_inc      (1,1) double {mustBeGreaterThanOrEqual(options.er_inc, 1), mustBeReal} = 1
            end
            
            obj.laser = laser;
            obj.pc_material = pc_material;
            
            if options.eta_opt ~= 0
                obj.default_eta_opt = options.eta_opt;
            end

            assert(options.theta_inc == 0, 'Non-normal incidence is not implemented.');

            obj.theta_inc = options.theta_inc;
            obj.er_inc = options.er_inc;
        end

        function set.theta_inc(obj, theta_inc)
            %SET.THETA_INC Summary of this method goes here
            %   theta_inc [\theta_{inc}] Incidence angle of laser beam at interface [deg]

            arguments
                obj
                theta_inc (1,1) double {mustBeInRange(theta_inc, 0, 90)}
            end

            assert(theta_inc == 0, 'Non-normal incidence is not implemented.');

            obj.theta_inc = theta_inc;
        end

        function set.er_inc(obj, er_inc)
            %SET.ER_INC Summary of this method goes here
            %   er_inc [\epsilon_{r,inc}] Relative permittivity of medium one (incidence medium) [-]

            arguments
                obj
                er_inc (1,1) double {mustBeGreaterThanOrEqual(er_inc, 1), mustBeReal}
            end

            obj.er_inc = er_inc;
        end
        
        function eta_tx = get.eta_tx(obj)
            %GET.ETA_TX Summary of this method goes here

            if ~isnan(obj.default_eta_opt)
                eta_tx = NaN;
                return;
            end

            eta_tx = obj.compute_transmittance();
        end
        
        function eta_absorp = get.eta_absorp(obj)
            %GET.ETA_ABSORP Summary of this method goes here

            if ~isnan(obj.default_eta_opt)
                eta_absorp = NaN;
                return;
            end

            eta_absorp = 1 - exp(- obj.pc_material.dimensions.Wz / obj.pc_material.absorp_len);
        end
        
        function eta_so = get.eta_so(obj)
            %GET.ETA_SO Summary of this method goes here

            if ~isnan(obj.default_eta_opt)
                eta_so = NaN;
                return;
            end

            Wx = obj.pc_material.dimensions.Wx;
            Wy = obj.pc_material.dimensions.Wy;
            
            % TODO: According to derivations in "Photoconductive Antenna Basics", the denominator must be 
            % sqrt(8) * sigma_rho; however according to the spatial distribution definition in 
            % "Characterization of an LT GaAs Bow-Tie Antenna", the denominator is computed to be 
            % sqrt(8 * pi) * sigma_rho, which gives closer results to expected efficiency for test 
            % parameters. Check, compare, and fix derivation for spatial distribution of laser.
            distrib_param = sqrt(8) * obj.laser.sigma_rho;
            eta_so = erf(Wx / distrib_param) * erf(Wy / distrib_param);
        end
        
        function eta_opt = get.eta_opt(obj)
            %GET.ETA_OPT Summary of this method goes here

            if ~isnan(obj.default_eta_opt)
                eta_opt = obj.default_eta_opt;
                return;
            end

            eta_opt = obj.eta_tx * obj.eta_absorp * obj.eta_so;
        end

        function R = compute_reflectance(obj)
            %COMPUTE_REFLECTANCE Summary of this method goes here

            n_inc = sqrt(obj.er_inc);
            n_pc  = sqrt(obj.pc_material.er);

            R = abs( (n_inc - n_pc) / (n_inc + n_pc) ) ^ 2;
        end

        function T = compute_transmittance(obj)
            %COMPUTE_TRANSMITTANCE Summary of this method goes here

            R = obj.compute_reflectance();
            T = 1 - R;
        end
    end
end


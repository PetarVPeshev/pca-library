classdef Efficiency < handle
    %EFFICIENCY Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess = protected, SetAccess = immutable)
        % TODO: Add a default size (1,1) for laser and photoconductor properties, the object constructors must
        % be changed to have zero argument inputs
        laser           Laser
        photoconductor  PhotoConductor
        eff_def_value   (1,1) double = NaN
    end

    properties (SetAccess = immutable)
        inc_angle       (1,1) double
        inc_rel_permit  (1,1) double
    end
    
    properties (Dependent, SetAccess = protected)
        eff_tx      (1,1) double
        eff_absorp  (1,1) double
        eff_so      (1,1) double
        eff_opt     (1,1) double
    end
    
    methods
        function obj = Efficiency(laser, photoconductor, options)
            %EFFICIENCY Construct an instance of this class
            %   laser Pulsed laser object
            %   photoconductor Photo-conductor object
            %   eff_opt [\eta_{opt}] Manually defined optical efficiency (not calculated from parameters) [-]
            %   inc_angle [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   inc_rel_permit [\epsilon_{r}^{inc}] Relative permittivity of medium one (incidence medium) [-]

            arguments
                laser                   (1,1) Laser
                photoconductor          (1,1) PhotoConductor
                options.eff_opt         (1,1) double {mustBeInRange(options.eff_opt, 0, 1)} = 0
                options.inc_angle       (1,1) double {mustBeInRange(options.inc_angle, 0, 90)} = 0
                options.inc_rel_permit  (1,1) double {mustBeGreaterThanOrEqual(options.inc_rel_permit, 1), ...
                                                     mustBeReal} = 1
            end
            
            obj.laser = laser;
            obj.photoconductor = photoconductor;
            
            if options.eff_opt ~= 0
                obj.eff_def_value = options.eff_opt;
            end

            obj.inc_angle = options.inc_angle;
            obj.inc_rel_permit = options.inc_rel_permit;
        end
        
        function eff_tx = get.eff_tx(obj)
            %GET.EFF_TX Summary of this method goes here

            if ~isnan(obj.eff_def_value)
                eff_tx = NaN;
                return;
            end

            eff_tx = obj.photoconductor.compute_transmittance(inc_angle = obj.inc_angle, ...
                                                              inc_rel_permit = obj.inc_rel_permit);
        end
        
        function eff_absorp = get.eff_absorp(obj)
            %GET.EFF_ABSORP Summary of this method goes here

            if ~isnan(obj.eff_def_value)
                eff_absorp = NaN;
                return;
            end

            eff_absorp = 1 - exp( - obj.photoconductor.dimens(3) / obj.photoconductor.absorp_coeff );
        end
        
        function eff_so = get.eff_so(obj)
            %GET.EFF_ABSORP Summary of this method goes here

            if ~isnan(obj.eff_def_value)
                eff_so = NaN;
                return;
            end

            eff_so = erf( obj.photoconductor.dimens(1) / (sqrt(8) * obj.laser.spatial_std) ) ...
                * erf( obj.photoconductor.dimens(2) / (sqrt(8) * obj.laser.spatial_std) );
        end
        
        function eff_opt = get.eff_opt(obj)
            %GET.EFF_ABSORP Summary of this method goes here

            if ~isnan(obj.eff_def_value)
                eff_opt = obj.eff_def_value;
                return;
            end

            eff_opt = obj.eff_tx * obj.eff_absorp * obj.eff_so;
        end
    end
end


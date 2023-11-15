classdef PhotoConductor < handle
    %PHOTOCONDUCTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        rec_time            (1,1) double
        scat_time           (1,1) double
        rel_permit          (1,1) double
        electr_mass_coeff   (1,1) double
        absorp_coeff        (1,1) double
        dimens              (1,3) double
    end

    properties (Dependent, SetAccess = protected)
        electr_eff_mass (1,1) double
        permit          (1,1) double
        dc_carr_mobil   (1,1) double
        vol             (1,1) double
    end
    
    methods
        function obj = PhotoConductor(dimens, rel_permit, charact_times, coeffs)
            %PHOTOCONDUCTOR Construct an instance of this class
            %   rec_time [\tau_{rec}] carrier recombination time [s]
            %   scat_time [\tau_{scat}] carrier scattering time [s]
            %   rel_permit [\epsilon_{r}] relative permittivity [-]
            %   electr_eff_mass [m_{eff}] electron effective mass [-]
            %   absorp_coeff [\alpha] absorption coefficient [m]
            %   dimens [w_{x}, w_{y}, w_{z}] photoconductor dimensions [m m m]

            arguments
                dimens                      (1,3) double {mustBePositive, mustBeReal}
                rel_permit                  (1,1) double {mustBeGreaterThanOrEqual(rel_permit, 1), mustBeReal}
                charact_times.rec_time      (1,1) double {mustBePositive, mustBeReal}
                charact_times.scat_time     (1,1) double {mustBePositive, mustBeReal}
                coeffs.electr_mass_coeff    (1,1) double {mustBePositive, mustBeReal}
                coeffs.absorp_coeff         (1,1) double {mustBePositive, mustBeReal}
            end

            obj.dimens = dimens;
            obj.rel_permit = rel_permit;
            obj.rec_time = charact_times.rec_time;
            obj.scat_time = charact_times.scat_time;
            obj.electr_mass_coeff = coeffs.electr_mass_coeff;
            obj.absorp_coeff = coeffs.absorp_coeff;
        end
        
        function electr_eff_mass = get.electr_eff_mass(obj)
            %GET.ELECTR_MASS Summary of this method goes here

            electr_eff_mass = obj.electr_mass_coeff * get_phys_const('ElectronMass');
        end

        function permit = get.permit(obj)
            %GET.PERMIT Summary of this method goes here

            permit = obj.rel_permit * get_phys_const('VacuumPermittivity');
        end

        function dc_carr_mobil = get.dc_carr_mobil(obj)
            %GET.DC_CARR_MOBIL Summary of this method goes here

            dc_carr_mobil = get_phys_const('ElectronCharge') * obj.scat_time / obj.electr_eff_mass;
        end

        function vol = get.vol(obj)
            %GET.VOLUME Summary of this method goes here

            vol = obj.dimens(1) * obj.dimens(2) * obj.dimens(3);
        end

        function reflectance = compute_reflectance(obj, options)
            %COMPUTE_REFL_COEFF Summary of this method goes here
            %   inc_angle [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   inc_rel_permit [\epsilon_{r}^{inc}] Relative permittivity of medium one (incidence medium) [-]

            arguments
                obj
                options.inc_angle       (1,1) double {mustBeInRange(options.inc_angle, 0, 90)} = 0
                options.inc_rel_permit  (1,1) double {mustBeGreaterThanOrEqual(options.inc_rel_permit, 1), ...
                                                      mustBeReal} = 1
            end

            assert(options.inc_angle == 0, 'Non-normal incidence is not implemented.');

            inc_refr_idx = sqrt(options.inc_rel_permit);
            pc_refr_idx = sqrt(obj.rel_permit);

            reflectance = abs( (inc_refr_idx - pc_refr_idx) / (inc_refr_idx + pc_refr_idx) ) ^ 2;
        end

        function transmittance = compute_transmittance(obj, options)
            %COMPUTE_TRANSM_COEFF Summary of this method goes here
            %   inc_angle [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   inc_rel_permit [\epsilon_{r}^{inc}] Relative permittivity of medium one (incidence medium) [-]

            arguments
                obj
                options.inc_angle       (1,1) double {mustBeInRange(options.inc_angle, 0, 90)} = 0
                options.inc_rel_permit  (1,1) double {mustBeGreaterThanOrEqual(options.inc_rel_permit, 1), ...
                                                      mustBeReal} = 1
            end

            reflectance = obj.compute_reflectance(inc_angle = options.inc_angle, ...
                                                  inc_rel_permit = options.inc_rel_permit);

            transmittance = 1 - reflectance;
        end
    end
end


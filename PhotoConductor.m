classdef PhotoConductor < handle
    %PHOTOCONDUCTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        tau_rec     (1,1) double
        tau_s       (1,1) double
        er          (1,1) double
        me_coef     (1,1) double
        absorp_len  (1,1) double
        dimensions  (1,1) struct = struct('Wx', 0, 'Wy', 0, 'Wz', 0)
    end

    properties (Dependent, SetAccess = protected)
        me          (1,1) double
        epsilon     (1,1) double
        mu_dc       (1,1) double
        V           (1,1) double
    end
    
    methods
        function obj = PhotoConductor(dimensions, er, charact_times, coeffs)
            %PHOTOCONDUCTOR Construct an instance of this class
            %   dimensions [w_{x}, w_{y}, w_{z}] photoconductor dimensions [m m m]
            %   er [\epsilon_{r}] relative permittivity [-]
            %   tau_rec [\tau_{rec}] carrier recombination time [s]
            %   tau_s [\tau_{s}] carrier scattering time [s]
            %   me_coeff [m_{e,eff}] electron effective mass [-]
            %   absorp_len [\alpha] absorption coefficient [m]

            arguments
                dimensions              (1,3) double {mustBePositive, mustBeReal}
                er                      (1,1) double {mustBeGreaterThanOrEqual(er, 1), mustBeReal}
                charact_times.tau_rec   (1,1) double {mustBePositive, mustBeReal}
                charact_times.tau_s     (1,1) double {mustBePositive, mustBeReal}
                coeffs.me_coef          (1,1) double {mustBePositive, mustBeReal}
                coeffs.absorp_len       (1,1) double {mustBePositive, mustBeReal}
            end

            obj.dimensions.Wx = dimensions(1);
            obj.dimensions.Wy = dimensions(2);
            obj.dimensions.Wz = dimensions(3);

            obj.er = er;

            obj.tau_rec = charact_times.tau_rec;
            obj.tau_s = charact_times.tau_s;

            obj.me_coef = coeffs.me_coef;
            obj.absorp_len = coeffs.absorp_len;
        end
        
        function me = get.me(obj)
            %GET.ME Summary of this method goes here

            me = obj.me_coef * get_phys_const('ElectronMass');
        end

        function epsilon = get.epsilon(obj)
            %GET.EPSILON Summary of this method goes here

            epsilon = obj.er * get_phys_const('VacuumPermittivity');
        end

        function mu_dc = get.mu_dc(obj)
            %GET.MU_DC Summary of this method goes here

            mu_dc = get_phys_const('ElectronCharge') * obj.tau_s / obj.me;
        end

        function V = get.V(obj)
            %GET.V Summary of this method goes here

            V = obj.dimensions.Wx * obj.dimensions.Wy * obj.dimensions.Wz;
        end

        function R = compute_reflectance(obj, options)
            %COMPUTE_REFLECTANCE Summary of this method goes here
            %   theta_inc [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   er_inc [\epsilon_{r,inc}] Relative permittivity of medium one (incidence medium) [-]

            arguments
                obj
                options.theta_inc   (1,1) double {mustBeInRange(options.theta_inc, 0, 90)} = 0
                options.er_inc      (1,1) double {mustBeGreaterThanOrEqual(options.er_inc, 1), mustBeReal} = 1
            end

            assert(options.theta_inc == 0, 'Non-normal incidence is not implemented.');

            n_inc = sqrt(options.er_inc);
            n_pc  = sqrt(obj.er);

            R = abs( (n_inc - n_pc) / (n_inc + n_pc) ) ^ 2;
        end

        function T = compute_transmittance(obj, options)
            %COMPUTE_TRANSMITTANCE Summary of this method goes here
            %   theta_inc [\theta_{inc}] Incidence angle of laser beam at interface [deg]
            %   er_inc [\epsilon_{r,inc}] Relative permittivity of medium one (incidence medium) [-]

            arguments
                obj
                options.theta_inc   (1,1) double {mustBeInRange(options.theta_inc, 0, 90)} = 0
                options.er_inc      (1,1) double {mustBeGreaterThanOrEqual(options.er_inc, 1), mustBeReal} = 1
            end

            R = obj.compute_reflectance(theta_inc = options.theta_inc, er_inc = options.er_inc);
            T = 1 - R;
        end
    end
end


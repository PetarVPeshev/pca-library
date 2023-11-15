classdef Laser < handle
    %LASERBEAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        wlen                (1,1) double
        freq                (1,1) double
        radius_half_width   (1,1) double
        pulse_half_pwr      (1,1) double
    end

    properties
        rep_period  (1,1) double
        avg_pwr     (1,1) double
    end

    properties (Dependent, SetAccess = protected)
        spatial_std     (1,1) double
        temporal_std    (1,1) double
        peak_pwr        (1,1) double
        peak_pwr_dens   (1,1) double
        pulse_energy    (1,1) double
        peak_e_field    (1,1) double
    end
    
    methods
        function obj = Laser(laser_wlen, rep_period, laser_avg_pwr, beam_distrib)
            %LASERBEAM Construct an instance of this class
            %   laser_wlen [\lambda_{c}] wavelength of the laser beam [m]
            %   rep_period [T_{L}] repetition rate of the laser beam [s]
            %   laser_avg_pwr [P_{avg}] average laser power [W]
            %   diam_half_width [D_{-3dB}] spatial -3dB diameter of the laser beam [m]
            %   pulse_half_pwr [\tau_{p}] temporal half-power pulse width (FWHM) [s]

            arguments
                laser_wlen                      (1,1) double {mustBePositive, mustBeReal}
                rep_period                      (1,1) double {mustBePositive, mustBeReal}
                laser_avg_pwr                   (1,1) double {mustBePositive, mustBeReal}
                beam_distrib.radius_half_width  (1,1) double {mustBePositive, mustBeReal}
                beam_distrib.pulse_half_pwr     (1,1) double {mustBePositive, mustBeReal}
            end
            
            obj.wlen = laser_wlen;
            obj.freq = get_phys_const('LightSpeed') / laser_wlen;
            obj.rep_period = rep_period;
            obj.avg_pwr = laser_avg_pwr;
            obj.radius_half_width = beam_distrib.radius_half_width;
            obj.pulse_half_pwr = beam_distrib.pulse_half_pwr;
        end
        
        function set.avg_pwr(obj, laser_avg_pwr)
            %SET.LASER_AVG_PWR Summary of this method goes here
            %   laser_avg_pwr [P_{avg}] average laser power

            arguments
                obj
                laser_avg_pwr (1,1) double {mustBePositive, mustBeReal}
            end

            obj.avg_pwr = laser_avg_pwr;
        end
        
        function set.rep_period(obj, rep_period)
            %SET.LASER_AVG_PWR Summary of this method goes here
            %   rep_period [T_{L}] repetition rate of the laser beam [s]

            arguments
                obj
                rep_period (1,1) double {mustBePositive, mustBeReal}
            end

            obj.rep_period = rep_period;
        end

        function spatial_std = get.spatial_std(obj)
            %GET.SPATIAL_STD Summary of this method goes here
            
            spatial_std = obj.radius_half_width / sqrt( 8 * log(2) );
        end

        function temporal_std = get.temporal_std(obj)
            %GET.TEMPORAL_STD Summary of this method goes here
            
            temporal_std = obj.pulse_half_pwr / sqrt( 8 * log(2) );
        end

        function peak_pwr = get.peak_pwr(obj)
            %GET.TEMPORAL_STD Summary of this method goes here
            
            peak_pwr = obj.avg_pwr * obj.rep_period / ( obj.temporal_std * sqrt(2 * pi) );
        end

        function peak_pwr_dens = get.peak_pwr_dens(obj)
            %GET.PEAK_PWR_DENS Summary of this method goes here

            peak_pwr_dens = obj.peak_pwr / ( 2 * pi * obj.spatial_std ^ 2 );
        end

        function pulse_energy = get.pulse_energy(obj)
            %GET.PULSE_ENER Summary of this method goes here

            pulse_energy = obj.avg_pwr * obj.rep_period;
        end

        function peak_e_field = get.peak_e_field(obj)
            %GET.PEAK_E_FIELD Summary of this method goes here
            
            peak_e_field = sqrt( 2 * obj.peak_pwr_dens * get_phys_const("VacuumImpedance") );
        end

        function time_vec = get_time_vec(obj, delta_t, time_limits)
            %GET_TIME_VEC Summary of this method goes here
            %   delta_t [\delta_{t}] Discrete time step [s]
            %   max_t [t_{max}] Maximum evaluation time [s]

            arguments
                obj
                delta_t             (1,1) {mustBePositive, mustBeReal}
                time_limits.min_t   (1,1) {mustBeNonpositive, mustBeReal} = 0
                time_limits.max_t   (1,1) {mustBeNonnegative, mustBeReal} = 0
            end

            if time_limits.min_t == 0
                time_limits.min_t = - 5 * obj.pulse_half_pwr;
            end

            if time_limits.max_t == 0
                time_limits.max_t = 5 * obj.pulse_half_pwr;
            end

            time_vec = time_limits.min_t : delta_t : time_limits.max_t;
        end

        function fig = plot_time_envelope(obj, delta_t, max_t, options)
            %PLOT_TIME_ENVELOPE Summary of this method goes here
            %   delta_t [\delta_{t}] Discrete time step [s]
            %   max_t [t_{max}] Maximum evaluation time [s]
            %   z [z] Elevation at which the temporal envelope is evaluated [m]
            %   line_style Plot line style
            %   line_width Plot line width
            %   display_name Plot legend name
            %   time_scale Time scale of x-axis and FWHM parameter

            arguments
                obj
                delta_t                 (1,1) {mustBePositive, mustBeReal}
                max_t                   (1,1) {mustBePositive, mustBeReal} = 5 * obj.pulse_half_pwr
                options.z               (1,1) {mustBeReal} = 0
                options.line_style      (1,:) {mustBeMember(options.line_style, {'-', '--', ':', '-.'})} = '-'
                options.line_width      (1,1) {mustBePositive, mustBeReal} = 1.5
                options.display_name    (1,:) {mustBeText} = 'time envelope'
                options.time_scale      (1,:) {mustBeMember(options.time_scale, {'ns', 'ps', 'fs'})} = 'ps'
            end

            c = get_phys_const('LightSpeed');
            if strcmp(options.time_scale, 'ns')
                time_scale_order = 1e9;
            elseif strcmp(options.time_scale, 'ps')
                time_scale_order = 1e12;
            else
                time_scale_order = 1e15;
            end
            
            time_vec = obj.get_time_vec(delta_t, max_t = max_t);
            temporal_envelope = zeros(1, length(time_vec));
            for rep_idx = 0 : 1 : floor(max_t / obj.rep_period)
                temporal_envelope = temporal_envelope + exp(- 0.5 * ( (time_vec - obj.rep_period * rep_idx ...
                    - options.z / c) / obj.temporal_std ) .^ 2);
            end

            fig = figure();
            plot(time_vec * time_scale_order, temporal_envelope, options.line_style, ...
                'LineWidth', options.line_width, 'DisplayName', options.display_name);
            grid on;
            xlabel(['t / ' options.time_scale]);
            ylabel('P(t) / P_{0}');
            title(['Laser Temporal Envelope @ \tau_{p} = ' num2str(obj.pulse_half_pwr * time_scale_order) ...
                ' ' options.time_scale ', T_{p} = ' num2str(obj.rep_period * 1e9) ' ns']);
        end
    end
end


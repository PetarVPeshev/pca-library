classdef Laser < handle
    %LASERBEAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        wlen        (1,1) double
        freq        (1,1) double
        R_3db       (1,1) double
        tau_p       (1,1) double
    end

    properties
        T           (1,1) double
        P           (1,1) double
    end

    properties (Dependent, SetAccess = protected)
        sigma_rho   (1,1) double
        sigma_t     (1,1) double
        P0          (1,1) double
        S0          (1,1) double
        Ep          (1,1) double
        E0          (1,1) double
    end
    
    methods
        function obj = Laser(wlen, rep_period, pwr_avg, beam_distrib)
            %LASERBEAM Construct an instance of this class
            %   wlen [\lambda_{c}] laser wavelength [m]
            %   rep_period [T] laser repetition rate [s]
            %   pwr_avg [P_{avg}] average laser power [W]
            %   R_3db [R_{-3dB}] spatial -3dB radius of the laser beam [m]
            %   tau_p [\tau_{p}] temporal half-power pulse width (FWHM) [s]

            arguments
                wlen                    (1,1) double {mustBePositive, mustBeReal}
                rep_period              (1,1) double {mustBePositive, mustBeReal}
                pwr_avg                 (1,1) double {mustBePositive, mustBeReal}
                beam_distrib.R_3db      (1,1) double {mustBePositive, mustBeReal}
                beam_distrib.tau_p      (1,1) double {mustBePositive, mustBeReal}
            end
            
            % TODO: Find better solution to adding the utility path in every class.
            add_dir_path('utils');
            
            obj.wlen = wlen;
            obj.freq = get_phys_const('LightSpeed') / wlen;

            obj.T = rep_period;
            obj.P = pwr_avg;

            obj.R_3db = beam_distrib.R_3db;
            obj.tau_p = beam_distrib.tau_p;
        end
        
        function set.P(obj, pwr_avg)
            %SET.P Summary of this method goes here
            %   pwr_avg [P_{avg}] average laser power

            arguments
                obj
                pwr_avg (1,1) double {mustBePositive, mustBeReal}
            end

            obj.P = pwr_avg;
        end
        
        function set.T(obj, rep_period)
            %SET.T Summary of this method goes here
            %   rep_period [T_{L}] repetition rate of the laser beam [s]

            arguments
                obj
                rep_period (1,1) double {mustBePositive, mustBeReal}
            end

            obj.T = rep_period;
        end

        function sigma_rho = get.sigma_rho(obj)
            %GET.SIGMA_RHO Summary of this method goes here
            
            sigma_rho = obj.R_3db / sqrt( 8 * log(2) );
        end

        function sigma_t = get.sigma_t(obj)
            %GET.SIGMA_T Summary of this method goes here
            
            sigma_t = obj.tau_p / sqrt( 8 * log(2) );
        end

        function P0 = get.P0(obj)
            %GET.P0 Summary of this method goes here
            
            P0 = obj.P * obj.T / ( obj.sigma_t * sqrt(2 * pi) );
        end

        function S0 = get.S0(obj)
            %GET.S0 Summary of this method goes here

            S0 = obj.P0 / ( 2 * pi * obj.sigma_rho ^ 2 );
        end

        function Ep = get.Ep(obj)
            %GET.Ep Summary of this method goes here

            Ep = obj.P * obj.T;
        end

        function E0 = get.E0(obj)
            %GET.E0 Summary of this method goes here
            
            E0 = sqrt( 2 * obj.S0 * get_phys_const("VacuumImpedance") );
        end

        function varargout = get_t_vec(obj, delta_t, t_lims)
            %GET_T_VEC Summary of this method goes here
            %   delta_t [\delta_{t}] Discrete time step [s]
            %   t_lims [t_{min}, t_{max}] Minimum and maximum time [s]

            arguments
                obj
                delta_t  (1,1) {mustBePositive, mustBeReal} = obj.tau_p / 100
                t_lims   (1,2) {mustBeReal} = [-10 20] * obj.tau_p
            end

            if t_lims(1) == t_lims(2) || t_lims(1) > t_lims(2)
                error(['Error. Parsed invalid time limits [' num2str(t_lims(1)) ',' num2str(t_lims(2)) '].']);
            end
            
            varargout{1} = t_lims(1) : delta_t : t_lims(2);

            if delta_t == obj.tau_p / 100
                varargout{2} = delta_t;
                varargout{3} = length(varargout{1});
            else
                varargout{2} = length(varargout{1});
            end
        end

        function [fig_obj, plot_obj] = plot_time_envelope(obj, delta_t, t_lims, options)
            %PLOT_TIME_ENVELOPE Summary of this method goes here
            %   delta_t [\delta_{t}] Discrete time step [s]
            %   t_lims [t_{min}, t_{max}] Minimum and maximum time [s]
            %   z [z] Elevation at which the temporal envelope is evaluated [m]
            %   LineStyle Line style
            %   LineWidth Line width
            %   DisplayName Plot name
            %   TimeScale Time scale of x-axis

            arguments
                obj
                delta_t                (1,1) {mustBePositive, mustBeReal} = obj.tau_p / 100
                t_lims                 (1,2) {mustBeReal} = [-10 20] * obj.tau_p
                options.z              (1,1) {mustBeReal} = 0
                options.LineStyle      (1,:) {mustBeMember(options.LineStyle, {'-', '--', ':', '-.'})} = '-'
                options.LineWidth      (1,1) {mustBePositive, mustBeReal} = 1.0
                options.DisplayName    (1,:) {mustBeText} = 'time envelope'
                options.TimeScale      (1,:) {mustBeMember(options.TimeScale, {'ns', 'ps', 'fs'})} = 'ps'
            end

            if strcmp(options.TimeScale, 'ns')
                t_vec_order = 1e9;
            elseif strcmp(options.TimeScale, 'ps')
                t_vec_order = 1e12;
            else
                t_vec_order = 1e15;
            end

            c = get_phys_const('LightSpeed');
            [t_vec, ~, num_pts] = obj.get_t_vec(delta_t, t_lims);
            
            t_envelop = zeros(1, num_pts);
            for idx_rep = 0 : 1 : floor(t_lims(2) / obj.T)
                t_param = t_vec - obj.T * idx_rep - options.z / c;
                t_envelop = t_envelop + exp(- 0.5 * (t_param / obj.sigma_t) .^ 2);
            end
            
            [tau_p_magn, tau_p_magn_order] = find_magnitude(obj.tau_p);
            [T_magn, T_magn_order] = find_magnitude(obj.T);

            fig_obj = figure();
            plot_obj = plot(t_vec * t_vec_order, t_envelop, options.LineStyle, ...
                'LineWidth', options.LineWidth, 'DisplayName', options.DisplayName);
            grid on;
            xlabel(['t / ' options.TimeScale]);
            ylabel('P(t) / P_{0}');
            title(['Laser Temporal Envelope @ \tau_{p} = ' num2str(round(obj.tau_p / tau_p_magn_order, 3)) ...
                ' ' tau_p_magn 's, T_{p} = ' num2str(round(obj.T / T_magn_order, 3)) ' ' T_magn 's']);
        end
    end
end

function add_dir_path(sub_dir)
%ADD_FOLDER_PATH Summary of this function goes here
%   sub_dir Sub-directory name
    dir_path = pwd();
    num_pca_lib = count(dir_path, '\pca-library');

    lib_path = extractBefore(dir_path, '\pca-library');
    for dir_num = 1 : 1 : num_pca_lib
        lib_path = append(lib_path, '\pca-library');
    end

    sub_dir_path = append(lib_path, '\', sub_dir);
    addpath(sub_dir_path);
end

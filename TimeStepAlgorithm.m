classdef TimeStepAlgorithm < handle
    %TIMESTEPALGORITHM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        time_vec (1,:) double
        delta_t (1,1) double
    end
    
    methods
        function obj = TimeStepAlgorithm(time_vec)
            %TIMESTEPALGORITHM Construct an instance of this class
            %   delta_t [\delta_{t}] Discretization time step [s]
            %   max_t [t_{max}] Maximum time of simulation [s] (first laser pulse has a peak power at 0s)

            arguments
                time_vec (1,:) double {mustBeReal}
            end
            
            obj.time_vec = time_vec;
            obj.delta_t = time_vec(2) - time_vec(1);
        end

%         function hm = compute_hm(k, sigma_t, tau_rec, options)
%             %COMPUTE_HM Summary of this function goes here
% 
%             arguments
%                 k           (1,1) double {mustBePositive}
%                 sigma_t     (1,1) double {mustBePositive}
%                 tau_rec     (1,1) double {mustBePositive}
%                 options.m   (1,1) double {mustBeInteger} = 0
%             end
% 
%             hm = exp(- 0.5 () ^ .2)
%         end
        
        function time_vec = get_time_vect(obj, laser)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            arguments
                obj
                laser (1,1) Laser
            end

            t_min = - laser.rep_period / 2;
            time_vec = t_min : obj.delta_t : obj.max_t;
        end
    end
end


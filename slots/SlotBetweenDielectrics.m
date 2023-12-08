classdef SlotBetweenDielectrics < handle
    %SLOTBETWEENDIELECTRICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        er_up   (1,1) double
        er_down (1,1) double
        f0      (1,1) double
        ws      (1,1) double
    end
    
    methods
        function obj = SlotBetweenDielectrics(er_up, er_down, f0, ws)
            %SLOTBETWEENDIELECTRICS Construct an instance of this class
            %   Detailed explanation goes here

            arguments
                er_up   (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(er_up, 1)}
                er_down (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(er_down, 1)}
                f0      (1,1) double {mustBeReal, mustBePositive}
                ws      (1,1) double {mustBeReal, mustBePositive}
            end
            
            % TODO: Find better solution to adding the utility path in every class.
            add_dir_path('utils');

            obj.er_up = er_up;
            obj.er_down = er_down;
            obj.f0 = f0;
            obj.ws = ws;
        end
        
        function D = compute_D(obj, kx)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            arguments
                obj
                kx (1,1) double {mustBeFinite}
            end
            
            k0 = 2 * pi * obj.f0 / get_phys_const('LightSpeed');
            eta_0 = get_phys_const('VacuumImpedance');

            ki = k0 * sqrt([obj.er_up obj.er_down]);
            args = obj.ws * sqrt(ki .^ 2 - kx ^ 2) / 4;

            const = 1 / (2 * k0 * eta_0);
            D = const * sum( (ki .^ 2 - kx ^ 2) .* besselj(0, args) .*  besselh(0, args) );
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


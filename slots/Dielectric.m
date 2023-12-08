classdef Dielectric < handle
    %DIELECTRIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        er (1,1) double
        n  (1,1) double
        c  (1,1) double
    end
    
    methods
        function obj = Dielectric(er)
            %DIELECTRIC Construct an instance of this class
            %   er [\epsilon_{r}] Relative permittivity of dielectric

            arguments
                er (1,1) double {mustBeReal, mustBeGreaterThanOrEqual(er, 1)}
            end
            
            % TODO: Find better solution to adding the utility path in every class.
            add_dir_path('utils');

            obj.er = er;
            obj.n = sqrt(er);
            obj.c = get_phys_const('LightSpeed') / obj.n;
        end
        
        function k = compute_k(obj, f)
            %COMPUTE_K Summary of this method goes here
            %   f [f] Frequency

            arguments
                obj
                f (1,1) double {mustBeReal, mustBePositive}
            end

            k = 2 * pi * f / obj.c;
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

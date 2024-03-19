% ENVIRONMENT SCRIPT
% This scripts adds the library paths to the MATLAB environment path.
% Run this script before executing code requiring this library.

env_dirs = {'', 'utils', 'slots', 'utils\validation-functions', 'dispersive_pca\coupling'};

add_env(env_dirs);
clear('env_dirs');

function lib_path = find_dir_path()
    dir_path = pwd();
    num_pca_lib = count(dir_path, '\pca-library');

    lib_path = extractBefore(dir_path, '\pca-library');
    for dir_num = 1 : 1 : num_pca_lib
        lib_path = append(lib_path, '\pca-library');
    end
end

function add_env(dirs)
    num_dirs = length(dirs);
    lib_path = find_dir_path();

    for dir_idx = 1 : num_dirs
        addpath(append(lib_path, "\", dirs{dir_idx}));
    end
end
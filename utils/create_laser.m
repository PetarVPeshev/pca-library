function laser = create_laser(laser_config, params)
%CREATE_LASER Summary of this function goes here
%   Detailed explanation goes here
    arguments
        laser_config  (1,1) struct {mustBeWithFields(laser_config, {'wlen', 'T', 'P', 'tau_p', 'R_3db'})}
        params.wlen   (1,1) double {mustBePositiveOrNan} = NaN
        params.T      (1,1) double {mustBePositiveOrNan} = NaN
        params.P      (1,1) double {mustBePositiveOrNan} = NaN
        params.tau_p  (1,1) double {mustBePositiveOrNan} = NaN
        params.R_3db  (1,1) double {mustBePositiveOrNan} = NaN
    end

    if isnan(params.wlen)
        params.wlen = laser_config.wlen;
    end
    if isnan(params.T)
        params.T = laser_config.T;
    end
    if isnan(params.P)
        params.P = laser_config.P;
    end
    if isnan(params.tau_p)
        params.tau_p = laser_config.tau_p;
    end
    if isnan(params.R_3db)
        params.R_3db = laser_config.R_3db;
    end

    laser = Laser(params.wlen, params.T, params.P, 'tau_p', params.tau_p, 'R_3db', params.R_3db);
end


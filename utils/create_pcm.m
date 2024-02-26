function pcm = create_pcm(pcm_config, params)
%CREATE_PCM Summary of this function goes here
%   Detailed explanation goes here
    arguments
        pcm_config      (1,1) struct {mustBeWithFields(pcm_config, {'wz', 'er', 'tau_rec', 'tau_s', ...
                                                                    'me_coef', 'alpha'})}
        params.d_gap    (1,1) double {mustBePositive}
        params.ws       (1,1) double {mustBePositive}
        params.wz       (1,1) double {mustBePositiveOrNan} = NaN
        params.er       (1,1) double {mustBeGreaterThanOrEqualOrNan(params.er, 1)} = NaN
        params.tau_rec  (1,1) double {mustBePositiveOrNan} = NaN
        params.tau_s    (1,1) double {mustBePositiveOrNan} = NaN
        params.me_coef  (1,1) double {mustBePositiveOrNan} = NaN
        params.alpha    (1,1) double {mustBePositiveOrNan} = NaN
    end

    if isnan(params.wz)
        params.wz = pcm_config.wz;
    end
    if isnan(params.er)
        params.er = pcm_config.er;
    end
    if isnan(params.tau_rec)
        params.tau_rec = pcm_config.tau_rec;
    end
    if isnan(params.tau_s)
        params.tau_s = pcm_config.tau_s;
    end
    if isnan(params.me_coef)
        params.me_coef = pcm_config.me_coef;
    end
    if isnan(params.alpha)
        params.alpha = pcm_config.alpha;
    end

    pcm = PhotoConductor([params.ws params.d_gap params.wz], params.er, 'tau_rec', params.tau_rec, ...
                         'tau_s', params.tau_s, 'me_coef', params.me_coef, 'absorp_len', params.alpha);
end


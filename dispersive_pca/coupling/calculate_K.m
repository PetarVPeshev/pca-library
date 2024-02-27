function K = calculate_K(laser, pcm, arg)
%CALCULATE_K Summary of this function goes here
%   Detailed explanation goes here
    arguments
        laser           Laser
        pcm             PhotoConductor
        arg.P_opt (1,1) double {mustBePositiveOrNan} = NaN
    end

    if isnan(arg.P_opt)
        arg.P_opt = laser.P;
    end
    
    qe = get_phys_const('ElectronCharge');
    h  = get_phys_const('PlanckConstant');

    nom = qe * pcm.mu_dc * laser.T * arg.P_opt;
    den = (pcm.dimensions.Wx ^ 2) * h * laser.freq * pcm.tau_s * sqrt(2 * pi) * laser.sigma_t;
    K = nom / den;
end


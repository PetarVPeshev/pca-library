function [laser, photoconductor] = get_default_config(args)
%GET_DEFAULT_CONFIG Summary of this function goes here
%   photoconductor - name of photoconductive material

    arguments
        args.photoconductor (1,1) string {mustBeMember(args.photoconductor, ["GaAs", "InAs", "GaInAs"])} = "GaAs"
    end

    run([pwd() '\default_config.m']);

    laser = Laser(laser_config.wlen, laser_config.T, laser_config.P, ...
                  R_3db = laser_config.R_3db, tau_p = laser_config.tau_p);
    
    assert(args.photoconductor == "GaAs", 'InAs and GaInAs are not in default parameter list.')

    photoconductor = PhotoConductor(GaAs_photoconductor_config.dimens, ...
        GaAs_photoconductor_config.rel_permit, ...
        rec_time = GaAs_photoconductor_config.rec_time, ...
        scat_time = GaAs_photoconductor_config.scat_time, ...
        electr_mass_coeff = GaAs_photoconductor_config.electr_eff_mass_coeff, ...
        absorp_coeff = GaAs_photoconductor_config.absorp_coeff);
end


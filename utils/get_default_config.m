function [laser, pc_material] = get_default_config(args)
%GET_DEFAULT_CONFIG Summary of this function goes here
%   photoconductor - name of photoconductive material

    arguments
        args.pc_material (1,1) string {mustBeMember(args.pc_material, ["GaAs", "InAs", "GaInAs"])} = "GaAs"
    end

    run([pwd() '\default_config.m']);

    laser = Laser(laser_config.wlen, laser_config.T, laser_config.P, ...
                  R_3db = laser_config.R_3db, tau_p = laser_config.tau_p);
    
    assert(args.pc_material == "GaAs", 'InAs and GaInAs are not in default parameter list.')

    pc_material = PhotoConductor(GaAs_config.dimensions, GaAs_config.er, tau_rec = GaAs_config.tau_rec, ...
                                 tau_s = GaAs_config.tau_s, me_coef = GaAs_config.me_coef, ...
                                 absorp_len = GaAs_config.absorp_len);
end


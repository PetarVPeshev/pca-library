function phys_const_value = get_phys_const(phys_const)
%GET_PHYS_CONST Summary of this function goes here
%   Detailed explanation goes here

    arguments
        phys_const (1,1) string {mustBeMember(phys_const, ["LightSpeed", ...
            "Boltzmann", "EarthRadius", "VacuumPermittivity", ...
            "VacuumPermeability", "ElectronCharge", "ElectronMass", ...
            "PlanckConstant", "VacuumImpedance"])}
    end

    if any(strcmp(phys_const, ["LightSpeed", "Boltzmann", "EarthRadius"]))
        phys_const_value = physconst(phys_const);
        return;
    elseif strcmp(phys_const, "VacuumPermittivity")
        phys_const_value = 8.8541878128 * 1e-12;
        return;
    elseif strcmp(phys_const, "VacuumPermeability")
        phys_const_value = 1.25663706212 * 1e-6;
        return;
    elseif strcmp(phys_const, "ElectronCharge")
        phys_const_value = 1.602176634 * 1e-19;
        return;
    elseif strcmp(phys_const, "ElectronMass")
        phys_const_value = 9.1093837015 * 1e-31;
        return;
    elseif strcmp(phys_const, "PlanckConstant")
        phys_const_value = 6.62607015 * 1e-34;
        return;
    elseif strcmp(phys_const, "VacuumImpedance")
        phys_const_value = sqrt( 1.25663706212 / 8.8541878128e-6 );
    end

end


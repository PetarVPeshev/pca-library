function slot = create_slot(slot_config, params)
%CREATE_SLOT Summary of this function goes here
%   Detailed explanation goes here
    arguments
        slot_config  (1,1) struct {mustBeWithFields(slot_config, {'er_up', 'er_dn', 'd_gap', 'ws'})}
        params.d_gap (1,1) double {mustBePositiveOrNan} = NaN
        params.ws    (1,1) double {mustBePositiveOrNan} = NaN
        params.er_up (1,1) double {mustBeGreaterThanOrEqualOrNan(params.er_up, 1)} = NaN
        params.er_dn (1,1) double {mustBeGreaterThanOrEqualOrNan(params.er_dn, 1)} = NaN
    end
    
    if isnan(params.d_gap)
        params.d_gap = slot_config.d_gap;
    end
    if isnan(params.ws)
        params.ws = slot_config.ws;
    end
    if isnan(params.er_up)
        params.er_up = slot_config.er_up;
    end
    if isnan(params.er_dn)
        params.er_dn = slot_config.er_dn;
    end

    slot = SlotInDielectrics(params.d_gap, params.ws, params.er_up, params.er_dn);
end

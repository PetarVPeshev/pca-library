function mustBeWithFields(struct_type, fields)
%MUSTBEWITHFIELDS Summary of this function goes here
%   Detailed explanation goes here
    struct_fields = fieldnames(struct_type);

    if ~all(ismember(struct_fields, fields), 'all')
        eid = "Struct:nonValidFields";
        msg = "Struct does not contain or has different fields from : " + string(fields(1));
        for idx = 2 : length(fields)
            msg = append(msg, ", ", fields(idx));
        end

        error(eid, msg);
    end
end

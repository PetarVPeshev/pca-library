function mustBeValidStruct(struct_type, fields)
%MUSTBEVALIDSTRUCT Summary of this function goes here
%   Detailed explanation goes here
    struct_fields = fieldnames(struct_type);
    num_fields    = length(fields);
    valid_fields  = fields(1 : 2 : num_fields);

    field_names = fields(1 : 2 : num_fields);
    for idx = 1 : 2 : num_fields
        if ~isnan(fields{idx + 1})
            field_names{idx - floor(idx / 2)} = [];
        end
    end
    field_names = field_names(~cellfun('isempty', field_names));
    struct_fields = struct_fields(~cellfun())

    if ~all(ismember(struct_fields, field_names), 'all')
        eid = "Struct:nonValidFields";
        msg = "Struct does not contain or has different fields from : " + string(field_names(1));
        for idx = 2 : length(field_names)
            msg = append(msg, ", ", field_names(idx));
        end

        error(eid, msg);
    end
end


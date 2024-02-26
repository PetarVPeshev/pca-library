function mustBeLengthOrNan(var, len)
%MUSTBESIZE Summary of this function goes here
%   Detailed explanation goes here

    if length(var) ~= len && ~isnan(var)
        eid = 'Length:notOfLengthOrNan';
        msg = 'Length of variable is not of specified and variable is not NaN';

        error(eid, msg)
    end
end


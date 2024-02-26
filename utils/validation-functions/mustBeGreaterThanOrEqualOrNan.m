function mustBeGreaterThanOrEqualOrNan(value, c)
%MUSTBEGREATERTHANOREQUALORNAN Summary of this function goes here
%   Detailed explanation goes here
    if ~isnan(value)
        mustBeGreaterThanOrEqual(value, c)
    end
end


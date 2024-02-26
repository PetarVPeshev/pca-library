function mustBePositiveOrNan(value)
%MUSTBEPOSITIVEORNAN Summary of this function goes here
%   Detailed explanation goes here
    if ~isnan(value)
        mustBePositive(value);
    end
end


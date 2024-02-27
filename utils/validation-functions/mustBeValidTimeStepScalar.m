function mustBeValidTimeStepScalar(A)
%MUSTBEVALIDTIMESTEPSCALAR Summary of this function goes here
%   Detailed explanation goes here
    mustBePositive(A);
    mustBeReal(A);
    mustBeNumeric(A);
end


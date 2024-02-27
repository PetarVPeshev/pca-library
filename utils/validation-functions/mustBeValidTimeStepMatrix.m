function mustBeValidTimeStepMatrix(A, B)
%MUSTBEVALIDTIMESTEPMATRIX Summary of this function goes here
%   Detailed explanation goes here
    mustBeNonempty(A);
    mustBeNonNan(A);

    if ~isequal(size(A), B)
        eid = 'Matrix:nonValidSize';
        msg = 'Time step matrix is not of valid size.';

        error(eid, msg);
    end
end


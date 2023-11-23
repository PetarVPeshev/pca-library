function mustBeEqualSize(a, b)
    %MUSTBEEQUALSIZE Summary of this function goes here
    %   Custom validation function for equal size between two matricies.
    %
    %   Source: www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html
    %   © 1994-2023 The MathWorks, Inc.

    if ~isequal(size(a), size(b))
        eid = 'Size:notEqual';
        msg = 'Size of first input must equal size of second input.';

        error(eid, msg)
    end
    
end

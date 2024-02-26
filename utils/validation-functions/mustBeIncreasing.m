function mustBeIncreasing(A)
%MUSTBEINCREASING Summary of this function goes here
%   Detailed explanation goes here
    num_elem = length(A);

    elem = A(1);
    for idx = 2 : num_elem
        if ~(A(idx) > elem)
            eid = 'NumericalArray:notIncreasing';
            msg = 'Elements in numerical array are not increasing';

            error(eid, msg);
        end

        elem = A(idx);
    end
end


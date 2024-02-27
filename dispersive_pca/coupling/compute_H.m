function H = compute_H(Z)
%COMPUTE_H Summary of this function goes here
%   Detailed explanation goes here
    arguments
        Z (2,2,:) double {mustBeNonNan}
    end

    Nf = size(Z, 3);
    H  = NaN(2, 2, Nf);

    H(1, 1, :) = 1 ./ Z(1, 1, :);
    H(1, 2, :) = Z(1, 2, :) ./ Z(1, 1, :);
    H(2, 1, :) = H(1, 2, :);
    H(2, 2, :) = H(1, 1, :);
end


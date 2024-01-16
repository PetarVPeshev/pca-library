function J0H02 = compute_J0H02(z)
%COMPUTE_J0H02 Summary of this method goes here
%   Detailed explanation goes here
    J0H02 = besselj(0, z) .* besselh(0, 2, z);
end

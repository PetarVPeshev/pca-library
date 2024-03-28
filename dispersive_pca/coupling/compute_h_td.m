function h_td = compute_h_td(t, f, H)
%COMPUTE_H_TD Summary of this function goes here
%   Detailed explanation goes here
    Nt = length(t);
    h_td = NaN(2, 2, Nt);

    h_td(1, 1, :) = 2 * real(eval_IFT(t, f, permute(H(1, 1, :), [1 3 2])));
    h_td(2, 2, :) = h_td(1, 1, :);
    h_td(1, 2, :) = 2 * real(eval_IFT(t, f, permute(H(1, 2, :), [1 3 2])));
    h_td(2, 1, :) = h_td(1, 2, :);
end

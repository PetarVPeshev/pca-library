function [H, W] = compute_H(Z, args)
%COMPUTE_H Summary of this function goes here
%   Detailed explanation goes here
    arguments
        Z      (2,2,:) double {mustBeNonNan}
        args.W (1,:)   double {mustBeNonNan} = double.empty(1, 0)
    end

    if isempty(args.W)
        args.W = (1 ./ Z(1, 1, :)) .^ 2;
    else
        args.W = permute(args.W, [1 3 2]);
    end

    H = Z .* args.W;
    W = permute(args.W, [1 3 2]);
end

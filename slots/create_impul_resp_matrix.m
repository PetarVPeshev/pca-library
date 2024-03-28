function [impulResp, weight] = create_impul_resp_matrix(timeArrayForResp, freqArray, impedMatrix, args)
%COMPUTE_H Summary of this function goes here
%   Detailed explanation goes here
    arguments
        timeArrayForResp (1, :) double {mustBeNonnegative}
        freqArray        (1, :) double {mustBeNonnegative}
        impedMatrix      (:, :, :) double {mustBeNonNan}
        args.Weight      (1, :)    double {mustBeNonNan}   = double.empty(1, 0)
    end

    if isempty(args.Weight)
        args.Weight = (1 ./ impedMatrix(1, 1, :)) .^ 2;
    else
        args.Weight = permute(args.Weight, [1 3 2]);
    end

    impulRespFreqDom = impedMatrix .* args.Weight;
    weightFreqDom    = permute(args.Weight, [1 3 2]);

    [numFeeds, ~, ~] = size(impedMatrix);
    numTime          = length(timeArrayForResp);

    impulRespTimeDomArray = NaN(1, numFeeds, numTime);
    for idxFeed = 1 : numFeeds
        feedImpulResp = permute(impulRespFreqDom(1, idxFeed, :), [1 3 2]);
        impulRespTimeDomArray(1, idxFeed, :) = 2 * real(eval_IFT(timeArrayForResp, freqArray, feedImpulResp));
    end

    impulRespTimeDom = NaN(numFeeds, numFeeds, numTime);
    for idxTime = 1 : numTime
        impulRespTimeDom(:, :, idxTime) = toeplitz(impulRespTimeDomArray(:, :, idxTime));
    end

    weightTimeDom = 2 * real(eval_IFT(timeArrayForResp, freqArray, weightFreqDom));

    impulResp = struct('freq_domain', impulRespFreqDom, 'time_domain', impulRespTimeDom);
    weight = struct('freq_domain', weightFreqDom, 'time_domain', weightTimeDom);
end


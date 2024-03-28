function impedMatrix = create_imped_matrix(freqPts, slotObj, args)
%CREATE_IMPED_MATRIX Create impedance matrix for uniformly spaced and equal dimension
%feeds along a slot between two semi-infinite homogeneous dielectrics
%   The function returns a 3D matrix with impedance values between
%   uniformly spaced elements with equal dimensions along a slot between 
%   two semi-infinite homogeneous dielectrics. Third dimension
%   corresponds to the provided frequency points, while first and second
%   dimensions correspond to the impedance between ith and jth elements.
%   The function takes the frequency points as first argument, a slot
%   between dielectrics object as a second argument, and three optional
%   name-value arguments specifying distance between feeds (default: 0),
%   number of feeds (default: 1), and maximum evaluation using total
%   integration (default: 0, use total integration path for all). In case
%   the feed is one, the function returns a 3D matrix with first two
%   dimensions equal to one. In case more than one feed are specified the
%   distance between feeds must be larger than the feed's length, else the
%   function returns an error.
    arguments
        freqPts           (1,:) double             {mustBePositive}
        slotObj                 SlotInDielectrics
        args.DistFeeds    (1,1) double             {mustBeReal, mustBeNonnegative} = 0
        args.NumFeeds     (1,1) double             {mustBeReal, mustBePositive}    = 1
        args.MaxTotalEval (1,1) double             {mustBeNonnegative}             = 0
    end

    if args.NumFeeds > 1 && args.DistFeeds < slotObj.d_gap
        error('DistFeeds:SmallerThanFeedLength', 'Distance between feeds is smaller than feed length.');
    end
    
    numFreq    = length(freqPts);
    impedArray = NaN(1, args.NumFeeds, numFreq);

    impedArray(1, 1, :) = slotObj.compute_zin(freqPts);
    for idxFeed = 2 : args.NumFeeds
        % TODO: Combine the total integration & pole capture functions into
        % one, use argument EvalType = ['Total', 'Residue'] with default
        % 'Total' to determine the type of evaluation; use only one
        % function call and define the string with one if statement

        distAtIdxFeed = args.DistFeeds * (idxFeed - 1);

        % impedEvalType = 'Total';
        % if distAtIdxFeed > args.MaxTotalEval && args.MaxTotalEval ~= 0
        %     impedEvalType = 'Residue';
        % end
        % 
        % impedArray(1, idxFeed, :) = evaluate_imped(freqPts, slotObj, distAtIdxFeed, impedEvalType);

        if distAtIdxFeed <= args.MaxTotalEval || args.MaxTotalEval == 0
            impedArray(1, idxFeed, :) = eval_Zm(freqPts, slotObj, distAtIdxFeed);
        else
            impedArray(1, idxFeed, :) = eval_Zmp(freqPts, slotObj, distAtIdxFeed);
        end
    end
    
    impedMatrix = NaN(args.NumFeeds, args.NumFeeds, numFreq);
    for idxFreq = 1 : numFreq
        impedArrayAtFreq = impedArray(:, :, idxFreq);
        impedMatrix(:, :, idxFreq) = toeplitz(real(impedArrayAtFreq)) + 1j * toeplitz(imag(impedArrayAtFreq));
    end
end

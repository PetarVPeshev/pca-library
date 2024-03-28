classdef TimeStepByImpressedCurrent < handle
    % TimeStepAlgorithmCoupling Add summary here
    %   Detailed explanation goes here

    properties (SetAccess = immutable)
        timeArray       (1, :)    double
        deltaTime       (1, 1)    double
        maxIdx          (1, 1)    double
    end

    properties
        numFeeds        (1, 1)    double
        constK          (1, 1)    double
        bias            (:, 1)    double
        weight          (1, :)    double
        impulResp       (:, :, :) double
        recTime         (1, 1)    double
        scatTime        (1, 1)    double
        feedDelay       (1, :)    double
        laserTimeStd    (1, 1)    double
    end

    properties (SetAccess = protected)
        weightAt0       (:, :)    double
        impulRespAt0    (:, :)    double
        pcmImpulResp    (:, :)    double
        imprCurrent     (:, :)    double
        % imprMatrix      (:, :, :) double
        decayForStep    (1, 1)    double
    end

    methods
        function obj = TimeStepByImpressedCurrent(timeArray, params)
            % Support name-value pair arguments when constructing object
            %   Detailed explanation goes here

            arguments
                timeArray           (1, :)    double
                params.numFeeds     (1, 1)    double = 0
                params.constK       (1, 1)    double = 0
                params.bias         (:, 1)    double = zeros(0, 1)
                params.weight       (1, :)    double = double.empty(1, 0)
                params.impulResp    (:, :, :) double = double.empty(1, 1, 0)
                params.recTime      (1, 1)    double = 0
                params.scatTime     (1, 1)    double = 0
                params.feedDelay    (1, :)    double = zeros(1, 0)
                params.laserTimeStd (1, 1)    double = 0
            end

            obj.timeArray    = timeArray;
            obj.deltaTime    = timeArray(2) - timeArray(1);
            obj.maxIdx       = length(timeArray);

            obj.numFeeds     = params.numFeeds;
            obj.constK       = params.constK;
            obj.bias         = params.bias;
            obj.impulResp    = params.impulResp;
            obj.weight       = params.weight;
            obj.recTime      = params.recTime;
            obj.scatTime     = params.scatTime;
            obj.feedDelay    = params.feedDelay;
            obj.laserTimeStd = params.laserTimeStd;
        end
    end

    methods(Access = protected)
        function setup(obj)
            obj.weightAt0    = obj.weight(1) * eye(obj.numFeeds);
            obj.impulRespAt0 = obj.impulResp(:, :, 1);
            obj.decayForStep = exp(- obj.deltaTime / obj.recTime) * exp(- obj.deltaTime / obj.scatTime);
            
            obj.pcmImpulResp = obj.compute_pcm_impul_resp(obj.timeArray, obj.constK, obj.laserTimeStd, ...
                obj.recTime, obj.feedDelay);
            obj.imprCurrent  = obj.compute_impr_current(obj.timeArray, obj.bias, obj.constK, ...
                obj.laserTimeStd, obj.recTime, obj.scatTime, obj.feedDelay);
            
            % TODO: Investigate the advantages / disadvantages of
            % pre-computing impressed matrix
            % obj.imprMatrix   = obj.create_impr_matrix(obj.imprCurrent, obj.impulResp);
        end

        function reset(obj)
            % TODO: Remove this function

            obj.weightAt0    = double.empty;
            obj.impulRespAt0 = double.empty;
            obj.decayForStep = double.empty;

            obj.pcmImpulResp = double.empty;
            obj.imprCurrent  = double.empty;
            % obj.imprMatrix   = double.empty;
        end

        function validate_properties(obj)
            mustBeValidTimeStepScalar(obj.numFeeds);
            mustBeValidTimeStepScalar(obj.constK);
            mustBeValidTimeStepScalar(obj.recTime);
            mustBeValidTimeStepScalar(obj.scatTime);
            mustBeValidTimeStepScalar(obj.laserTimeStd);

            mustBeValidTimeStepMatrix(obj.bias, [obj.numFeeds 1]);
            mustBeValidTimeStepMatrix(obj.weight, [1 obj.maxIdx]);
            mustBeValidTimeStepMatrix(obj.impulResp, [obj.numFeeds obj.numFeeds obj.maxIdx]);
            mustBeValidTimeStepMatrix(obj.feedDelay, [1 obj.numFeeds]);
        end
    end

    methods
        function [voltage, current] = simulate(obj)
            obj.validate_properties();
            % obj.reset();
            obj.setup();

            voltage    = NaN(obj.numFeeds, obj.maxIdx);
            intCurrent = voltage;

            % Initial conditions
            pcmImpulRespAtIdx = diag(obj.pcmImpulResp(:, 1));
            imprMatrix        = obj.create_impr_matrix(obj.imprCurrent(:, 1), obj.impulResp(:, :, 1));
            voltage(:, 1)     = (obj.weightAt0 + obj.impulRespAt0 * pcmImpulRespAtIdx) \ (imprMatrix * ones(obj.numFeeds, 1));
            intCurrent(:, 1)  = pcmImpulRespAtIdx * voltage(:, 1);

            for timeIdx = 2 : obj.maxIdx
                pcmImpulRespAtIdx = diag(obj.pcmImpulResp(:, timeIdx));
                imprMatrix        = obj.create_impr_matrix(obj.imprCurrent(:, 1 : timeIdx), obj.impulResp(:, :, 1 : timeIdx));
                intMatrix         = obj.create_int_matrix_at_idx(intCurrent(:, 1 : timeIdx - 1), obj.impulResp(:, :, 2 : timeIdx));
                voltageMatrix     = obj.create_v_matrix_at_idx(voltage(:, 1 : timeIdx - 1), obj.weight(1, 2 : timeIdx));

                % Transient voltage at time index
                voltage(:, timeIdx)    = obj.compute_voltage_at_idx(obj.decayForStep, obj.weightAt0, ...
                    obj.impulRespAt0, pcmImpulRespAtIdx, intCurrent(:, timeIdx - 1), ...
                    imprMatrix, intMatrix, voltageMatrix);

                % Internal current at time index
                intCurrent(:, timeIdx) = obj.compute_int_current_at_idx(obj.decayForStep, ...
                    intCurrent(:, timeIdx - 1), voltage(:, timeIdx), pcmImpulRespAtIdx);
            end

            current = obj.imprCurrent - intCurrent;
        end
    end

    methods (Static)
        function pcmImpulResp = compute_pcm_impul_resp(timeArray, constK, laserTimeStd, scatTime, feedDelay)
            %COMPUTE_FSM Summary of this method goes here
            %   Detailed explanation goes here
            maxIdx    = length(timeArray);
            numFeeds  = length(feedDelay);
            deltaTime = timeArray(2) - timeArray(1);

            [N, M]      = meshgrid(1 : maxIdx, 1 : maxIdx);
            timeMatrixN = timeArray(N);
            timeMatrixM = timeArray(M);

            pcmImpulResp = NaN(numFeeds, maxIdx);
            parfor feedIdx = 1 : numFeeds
                pcmRespMatrix = exp(- 0.5 * ((timeMatrixN - feedDelay(feedIdx)) / laserTimeStd) .^ 2) ...
                                .* exp(- (timeMatrixM - timeMatrixN) / scatTime);

                pcmImpulResp(feedIdx, :) = sum(tril(pcmRespMatrix), 2)';
            end

            pcmImpulResp = (deltaTime ^ 2) * constK * pcmImpulResp;
        end

        function imprCurrent = compute_impr_current(timeArray, bias, constK, laserTimeStd, recTime, ...
                scatTime, feedDelay)
            %COMPUTE_I_IMPR Summary of this method goes here
            %   Detailed explanation goes here
            maxIdx    = length(timeArray);
            numFeeds  = length(feedDelay);
            deltaTime = timeArray(2) - timeArray(1);

            [N, M]      = meshgrid(1 : maxIdx, 1 : maxIdx);
            timeMatrixN = timeArray(N);
            timeMatrixM = timeArray(M);

            imprCurrent = NaN(numFeeds, maxIdx);
            parfor feedIdx = 1 : numFeeds
                pcmRespMatrix = exp(- 0.5 * ((timeMatrixN - feedDelay(feedIdx)) / laserTimeStd) .^ 2) .* exp(- (timeMatrixM - timeMatrixN) / recTime) ...
                                .* (1 - exp(- (timeMatrixM - timeMatrixN) / scatTime));

                imprCurrent(feedIdx, :) = sum(tril(pcmRespMatrix), 2)';
            end

            imprCurrent = deltaTime * scatTime * constK * bias .* imprCurrent;
        end

        function imprMatrix = create_impr_matrix(imprCurrent, impulResp)
            %COMPUTE_C_IMPR Summary of this method goes here
            %   Detailed explanation goes here
            [numFeeds, ~, ~] = size(impulResp);

            impulResp   = permute(impulResp, [2 3 1]);
            imprCurrent = repmat(imprCurrent, [1 1 numFeeds]);

            imprMatrix = permute(sum(imprCurrent .* fliplr(impulResp), 2), [3 1 2]);
        end

        function intMatrixAtIdx = create_int_matrix_at_idx(intCurrent, impulResp)
            %COMPUTE_CM_INT Summary of this method goes here
            %   Detailed explanation goes here
            [numFeeds, ~, ~] = size(impulResp);

            impulResp  = permute(impulResp, [2 3 1]);
            intCurrent = repmat(intCurrent, [1 1 numFeeds]);

            intMatrixAtIdx = permute(sum(intCurrent .* fliplr(impulResp), 2), [3 1 2]);
        end

        function voltageMatrixAtIdx = create_v_matrix_at_idx(voltage, weight)
            %COMPUTE_CM_V Summary of this method goes here
            %   Detailed explanation goes here
            [numFeeds, ~] = size(voltage);

            weight = repmat(weight, [numFeeds 1]);

            voltageMatrixAtIdx = sum(voltage .* fliplr(weight), 2);
        end

        function voltageAtIdx = compute_voltage_at_idx(decayForStep, weightAt0, impulRespAt0, ...
                pcmImpulRespAtIdx, prevIntCurrent, imprMatrixAtIdx, intMatrixAtIdx, voltageMatrixAtIdx)
            %COMPUTE_VM Summary of this method goes here
            %   Detailed explanation goes here
            numFeeds  = size(prevIntCurrent, 1);
            oneVector = ones(numFeeds, 1);

            voltageAtIdx = (weightAt0 + impulRespAt0 * pcmImpulRespAtIdx) ...
                            \ ((imprMatrixAtIdx - intMatrixAtIdx) * oneVector ...
                               - decayForStep * impulRespAt0 * prevIntCurrent - voltageMatrixAtIdx);
        end

        function intCurrentAtIdx = compute_int_current_at_idx(decayForStep, prevIntCurrent, ...
                voltageAtIdx, pcmImpulRespAtIdx)
            %COMPUTE_IM Summary of this method goes here
            %   Detailed explanation goes here
            intCurrentAtIdx = decayForStep * prevIntCurrent + pcmImpulRespAtIdx * voltageAtIdx;
        end
    end
end
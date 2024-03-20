classdef ImpedanceMatrixTest < matlab.unittest.TestCase

    properties
        % TEST PARAMETERS
        numFFT        (1, 1) double = 4001       % NUMBER OF FFT POINTS : 4001
        limitsFFT     (1, 2) double = [1e9 3e12] % FREQUENCY LIMITS [minFrequency maxFrequency] : [0.001 3] THz
        feedDistance  (1, 1) double = 100e-6     % DISTANCE B/N FEEDS : 100 um
        feedLength    (1, 1) double = 10e-6      % FEED LENGTH : 10 um
        slotWidth     (1, 1) double = 5e-6       % SLOT WIDTH : 5 um
        relPermitUp   (1, 1) double = 11.7       % RELATIVE PERMITTIVITY UP : 11.7 (SILICON)
        relPermitDown (1, 1) double = 1          % RELATIVE PERMITTIVITY DOWN : 1 (FREE SPACE)
    end

    properties
        disperEqn = @(k0, K1, K2, width) ((K1 .^ 2) .* compute_J0H02(K1 * width / 4) ...
                                        + (K2 .^ 2) .* compute_J0H02(K2 * width / 4)) ...
                                        / (2 * k0 * get_phys_const('VacuumImpedance'));
    end

    properties (ClassSetupParameter)
        FFT           (1, :) double
        selfImpedance (1, :) double
        slot                 SlotInDielectrics
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function setEnvironment(testCase)
            originalPath = path;
            testCase.addTeardown(@path, originalPath);
            run('..\..\dispersive_pca\coupling\env.m');
        end

        function setDependentParameters(testCase)
            testCase.slot = SlotInDielectrics(testCase.feedLength, testCase.slotWidth, ...
                                              testCase.relPermitUp, testCase.relPermitDown);
        end

        function evaluateSelfImpedance(testCase)
            testCase.selfImpedance = NaN(1, testCase.numFFT);
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods

        function unimplementedTest(testCase)
            testCase.verifyFail("Unimplemented test");
        end
    end

end
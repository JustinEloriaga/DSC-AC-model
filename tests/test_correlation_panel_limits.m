function tests = test_correlation_panel_limits
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'scripts'));
testCase.TestData.root = root;
end

function testPositiveRange(testCase)
verifyEqual(testCase,correlation_panel_limits([0.2 0.4]),[0.19 0.41],'AbsTol',1e-12);
end

function testNegativeRange(testCase)
verifyEqual(testCase,correlation_panel_limits([-0.7 -0.3]),[-0.71 -0.29],'AbsTol',1e-12);
end

function testCrossZeroRange(testCase)
verifyEqual(testCase,correlation_panel_limits([-0.2 0.4]),[-0.21 0.41],'AbsTol',1e-12);
end

function testConstantRange(testCase)
verifyEqual(testCase,correlation_panel_limits(0.25,0.25),[0.24 0.26],'AbsTol',1e-12);
end

function testBandExtremaDefineRange(testCase)
medianPath = [0.1; 0.2];
lower = [-0.4; -0.1];
upper = [0.3; 0.8];
verifyEqual(testCase,correlation_panel_limits(medianPath,lower,upper), ...
    [-0.41 0.81],'AbsTol',1e-12);
end

function testNaNBandsAreIgnored(testCase)
verifyEqual(testCase,correlation_panel_limits([0.1; 0.2],[NaN; -0.3],[0.4; NaN]), ...
    [-0.31 0.41],'AbsTol',1e-12);
end

function testAllNonfiniteValuesReject(testCase)
verifyError(testCase,@()correlation_panel_limits([NaN Inf -Inf]),'dsc:PlotLimits');
verifyError(testCase,@()correlation_panel_limits(),'dsc:PlotLimits');
end

function testInvalidTypeRejects(testCase)
verifyError(testCase,@()correlation_panel_limits("0.2"),'dsc:PlotLimits');
verifyError(testCase,@()correlation_panel_limits(1+2i),'dsc:PlotLimits');
end

function tests=test_credible_bands
tests=functiontests(localfunctions);
end
function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'toolbox'));
end
function testPercentilesAndPerChainBurnin(testCase)
a=reshape([-.99 -.8 -.4 0 .4],1,1,5);
b=reshape([.99 .2 .6 .8],1,1,4);
out=dsc_credible_bands({a,b},{1:5,1:4},[2,1]);
expected=prctile([-.4 0 .4 .2 .6 .8],[16 50 84]);
verifyEqual(testCase,[out.lower out.median out.upper],expected,'AbsTol',1e-14);
verifyEqual(testCase,out.draws_per_chain,[3 3]);
verifyEqual(testCase,out.retained_iterations,{3:5,2:4});
verifyFalse(testCase,out.publication_ready);
end
function testPostburnChunksNotTrimmedTwice(testCase)
x=reshape([-.5 0 .5],1,1,3);
out=dsc_credible_bands({x},{[101 103 105]},100);
verifyEqual(testCase,out.draws_per_chain,3);
verifyEqual(testCase,out.median,0);
verifyEqual(testCase,out.probability,.68);
end
function testDatePairShapeAndBounds(testCase)
x=zeros(2,3,10);
for d=1:10, x(:,:,d)=d/20+[0 .05 .1;-.1 -.05 0]; end
out=dsc_credible_bands({x},{1:10},0);
verifySize(testCase,out.lower,[2,3]);
verifyEqual(testCase,out.lower,prctile(x,16,3),'AbsTol',1e-14);
verifyEqual(testCase,out.upper,prctile(x,84,3),'AbsTol',1e-14);
verifyTrue(testCase,all(out.lower<=out.median & out.median<=out.upper,'all'));
end
function testNoPilotBands(testCase)
verifyError(testCase,@()dsc_credible_bands({zeros(2,3,0)},{[]},20),'dsc:NoPosteriorDraws');
verifyError(testCase,@()dsc_credible_bands({zeros(2,3,13)},{1:13},20),'dsc:NoPosteriorDraws');
end
function testRejectBadInputs(testCase)
x=zeros(2,3,2);
verifyError(testCase,@()dsc_credible_bands({x},{[2 2]},0),'dsc:BandsIterations');
x(1)=1.5;
verifyError(testCase,@()dsc_credible_bands({x},{1:2},0),'dsc:BandsCorrelations');
end

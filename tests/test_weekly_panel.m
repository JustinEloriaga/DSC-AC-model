function tests=test_weekly_panel
%TEST_WEEKLY_PANEL Regression tests for the calendar and observation contract.
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'toolbox'));
testCase.TestData.root=root;
end

function testCanonicalCountsAndTransformation(testCase)
cfg=struct('source_file',fullfile(testCase.TestData.root,'data', ...
    'yad_tickers_no_usdcnh_from_20030804.csv'));
p=prepare_weekly_panel(cfg);
verifySize(testCase,p.returns,[1204,14]);
verifySize(testCase,p.levels,[1205,14]);
verifyEqual(testCase,p.dates([1,end]),[datetime(2003,8,15);datetime(2026,9,4)]);
verifyEqual(testCase,p.level_dates(1),datetime(2003,8,8));
verifyEqual(testCase,diff(p.dates),repmat(days(7),1203,1));
verifyEqual(testCase,sum(p.carried(:)),343);
verifyEqual(testCase,max(p.age_days(:)),7);
verifyEqual(testCase,sum(p.closure_level_mask(:)),1);
verifyEqual(testCase,sum(p.closure_return_mask(:)),2);
verifyTrue(testCase,all(p.observation_mask(:)));
verifyTrue(testCase,all(isfinite(p.returns(:))));
verifyTrue(testCase,all(p.source_dates<=p.level_dates,'all'));
verifyEqual(testCase,p.returns,100*diff(log(p.levels)), 'AbsTol',1e-12);
verifyEqual(testCase,numel(p.pair_i),91);
verifyEqual(testCase,[p.pair_i(1),p.pair_j(1)],[2,1]);
verifyFalse(testCase,any(contains(p.tickers,'USDCNH')));
verifyEqual(testCase,p.source_hash, ...
    'f11259e4c6347a54204cecb4f5ac3363c06ef2b05f787ab33ed9ab99e04c5bd7');
verifyEqual(testCase,p.excluded_final_rows,4);
end

function testFridayMissingUsesOwnThursdayOnly(testCase)
dates=weekdaysBetween(datetime(2020,1,6),datetime(2020,1,31));
levels=syntheticLevels(numel(dates));
friday=dates==datetime(2020,1,17);
levels(friday,2)=NaN;
[path,cleanup]=fixture(dates,levels); %#ok<ASGLU>
p=prepare_weekly_panel(struct('source_file',path));
k=find(p.level_dates==datetime(2020,1,17));
verifyEqual(testCase,p.source_dates(k,2),datetime(2020,1,16));
verifyEqual(testCase,p.age_days(k,2),1);
verifyEqual(testCase,p.source_dates(k,1),datetime(2020,1,17));
verifyEqual(testCase,p.levels(k,2),levels(dates==datetime(2020,1,16),2));
verifyEqual(testCase,p.returns(k-1,2), ...
    100*log(levels(dates==datetime(2020,1,16),2)/levels(dates==datetime(2020,1,10),2)), ...
    'AbsTol',1e-12);
end

function testCalendarYearBoundary(testCase)
dates=weekdaysBetween(datetime(2020,12,21),datetime(2021,1,8));
[path,cleanup]=fixture(dates,syntheticLevels(numel(dates))); %#ok<ASGLU>
p=prepare_weekly_panel(struct('source_file',path));
verifyEqual(testCase,p.level_dates, ...
    [datetime(2020,12,25);datetime(2021,1,1);datetime(2021,1,8)]);
verifyEqual(testCase,diff(p.dates),days(7));
end

function testUnexplainedWholeWeekGapIsRejected(testCase)
dates=weekdaysBetween(datetime(2020,1,6),datetime(2020,1,31));
levels=syntheticLevels(numel(dates));
levels(dates>=datetime(2020,1,13)&dates<=datetime(2020,1,17),2)=NaN;
[path,cleanup]=fixture(dates,levels); %#ok<ASGLU>
verifyError(testCase,@()prepare_weekly_panel(struct('source_file',path)), ...
    'weekly:UnexplainedGap');
end

function testDocumentedClosureMaskPreservesTimeAndValues(testCase)
dates=weekdaysBetween(datetime(2019,4,22),datetime(2019,5,17));
levels=syntheticLevels(numel(dates));
levels(dates>=datetime(2019,4,29)&dates<=datetime(2019,5,6),10)=NaN;
[path,cleanup]=fixture(dates,levels); %#ok<ASGLU>
cfg=struct('source_file',path,'closure_policy','mask');
p=prepare_weekly_panel(cfg);
verifyEqual(testCase,find(~p.observation_mask(:,10)),[1;2]);
verifyEqual(testCase,sum(~p.observation_mask(:)),2);
verifyEqual(testCase,p.returns(1,10),0,'AbsTol',0);
verifyTrue(testCase,all(p.returns(1,[1:9,11:14])>0));
verifyEqual(testCase,p.dates, ...
    [datetime(2019,5,3);datetime(2019,5,10);datetime(2019,5,17)]);
verifyTrue(testCase,all(isfinite(p.returns(:))));
end

function testIncompleteFinalWeekExcluded(testCase)
dates=weekdaysBetween(datetime(2020,1,6),datetime(2020,1,23));
[path,cleanup]=fixture(dates,syntheticLevels(numel(dates))); %#ok<ASGLU>
p=prepare_weekly_panel(struct('source_file',path));
verifyEqual(testCase,p.level_dates,[datetime(2020,1,10);datetime(2020,1,17)]);
verifyEqual(testCase,p.excluded_final_rows,4);
end

function testIncompleteInitialWeekExcluded(testCase)
dates=weekdaysBetween(datetime(2020,1,8),datetime(2020,1,31));
[path,cleanup]=fixture(dates,syntheticLevels(numel(dates))); %#ok<ASGLU>
p=prepare_weekly_panel(struct('source_file',path));
verifyEqual(testCase,p.level_dates(1),datetime(2020,1,17));
end

function testDuplicateDatesRejected(testCase)
dates=weekdaysBetween(datetime(2020,1,6),datetime(2020,1,31));
dates(4)=dates(3);
[path,cleanup]=fixture(dates,syntheticLevels(numel(dates))); %#ok<ASGLU>
verifyError(testCase,@()prepare_weekly_panel(struct('source_file',path)), 'weekly:Dates');
end

function testNonpositiveLevelsRejected(testCase)
dates=weekdaysBetween(datetime(2020,1,6),datetime(2020,1,31));
levels=syntheticLevels(numel(dates)); levels(2,3)=0;
[path,cleanup]=fixture(dates,levels); %#ok<ASGLU>
verifyError(testCase,@()prepare_weekly_panel(struct('source_file',path)), 'weekly:Levels');
end

function dates=weekdaysBetween(first,last)
dates=(first:caldays(1):last)';
dates=dates(~ismember(weekday(dates),[1,7]));
end

function levels=syntheticLevels(n)
levels=10+(1:n)'*(1:14)/100;
end

function [path,cleanup]=fixture(dates,levels)
% Generated test data live in one exact temporary CSV, never the input folder.
path=[tempname,'.csv'];
cleanup=onCleanup(@()delete(path));
tickers=["AUDUSD Index","BCOMCOT Index","BCOMGCTR Index","EURUSD Index", ...
    "GBPUSD Index","I02981JP Index","LEATTREU Index","LSG1TRGU Index", ...
    "LUATTRUU Index","NKYTR Index","SPXT Index","SX5T Index", ...
    "TUKXG Index","USDJPY Index"];
fid=fopen(path,'w');
closeFile=onCleanup(@()fclose(fid));
fprintf(fid,',%s\n',strjoin(repmat("px_last",1,14),','));
fprintf(fid,'security,%s\n',strjoin(tickers,','));
fprintf(fid,'date%s\n',repmat(',',1,14));
for t=1:numel(dates)
    fields=compose("%.15g",levels(t,:)); fields(isnan(levels(t,:)))="";
    fprintf(fid,'%s,%s\n',string(dates(t),'dd/MM/yyyy'),strjoin(fields,','));
end
clear closeFile
end

function tests=test_convergence_diagnostics
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'toolbox'));
testCase.TestData.root=tempname; mkdir(testCase.TestData.root);
end

function teardownOnce(testCase)
if isfolder(testCase.TestData.root), rmdir(testCase.TestData.root,'s'); end
end

function testIidHasGoodDiagnostics(testCase)
rng(810,'twister'); x=randn(4000,4); d=dsc_mcmc_diagnostics(x);
verifyTrue(testCase,d.available); verifyLessThan(testCase,d.rhat,1.01);
verifyGreaterThan(testCase,d.ess_bulk,10000); verifyGreaterThan(testCase,d.ess_tail,8000);
verifyLessThan(testCase,d.mcse_sd_ratio,.02);
verifyEqual(testCase,d.mcse_mean,std(x(:))/sqrt(1/d.mcse_sd_ratio^2),'RelTol',1e-12);
end

function testLocationAndScaleFailures(testCase)
rng(811,'twister'); x=randn(2500,4);
location=dsc_mcmc_diagnostics(x+[0 2 0 2]);
scale=dsc_mcmc_diagnostics(x.*[.1 1 4 12]);
verifyGreaterThan(testCase,location.rhat,1.1);
verifyGreaterThan(testCase,scale.rhat,1.1); % folded Rhat must find scale differences
end

function testAutocorrelationReducesEss(testCase)
rng(812,'twister'); innovation=randn(5000,4);
x=filter(1,[1 -.98],innovation); x=x(1001:end,:);
d=dsc_mcmc_diagnostics(x); iid=dsc_mcmc_diagnostics(innovation(1001:end,:));
verifyTrue(testCase,d.available);
verifyLessThan(testCase,d.ess_bulk,500); verifyLessThan(testCase,d.ess_tail,1500);
verifyGreaterThan(testCase,d.mcse_sd_ratio,4*iid.mcse_sd_ratio);
end

function testRawMeanMcseUsesRawEssAndDirectAutocovariance(testCase)
% Independent direct lag-dot-products oracle, applying the published
% positive/monotone sequence to raw split draws rather than rank ESS.
rng(813,'twister'); x=exp(filter(1,[1 -.8],randn(1200,4)));
d=dsc_mcmc_diagnostics(x); expected=direct_reference_ess(x);
verifyEqual(testCase,d.mcse_mean,std(x(:))/sqrt(expected),'RelTol',1e-10);
verifyEqual(testCase,d.mcse_sd_ratio,1/sqrt(expected),'RelTol',1e-10);
verifyGreaterThan(testCase,abs(expected-d.ess_bulk),10);
end

function testOddDrawCountAndTiedRanks(testCase)
rng(814,'twister'); x=round(3*randn(1201,4));
d=dsc_mcmc_diagnostics(x);
verifyTrue(testCase,d.available);
verifyEqual(testCase,d.mcse_sd_ratio,1/sqrt(direct_reference_ess(x)),'RelTol',1e-10);
% Strictly increasing transforms preserve both pooled ranks and tail ESS.
transformed=dsc_mcmc_diagnostics(exp(x/10));
verifyEqual(testCase,d.ess_bulk,transformed.ess_bulk,'RelTol',1e-12);
verifyEqual(testCase,d.ess_tail,transformed.ess_tail,'RelTol',1e-12);
end

function testUnavailableAndSingleChain(testCase)
rng(815,'twister'); samples={ones(100,4),[ones(100,1) randn(100,3)], ...
    randn(5,4),[NaN;randn(99,1)],[Inf;randn(99,1)]};
for k=1:numel(samples)
    d=dsc_mcmc_diagnostics(samples{k}); verifyFalse(testCase,d.available);
    verifyTrue(testCase,isnan(d.rhat)); verifyTrue(testCase,isnan(d.ess_bulk));
    verifyTrue(testCase,isnan(d.ess_tail)); verifyTrue(testCase,isnan(d.mcse_mean));
end
d=dsc_mcmc_diagnostics(randn(1500,1));
verifyTrue(testCase,d.available); verifyEqual(testCase,d.original_chains,1);
end

function testShortSplitBoundaryMatchesPosteriorReference(testCase)
% R's 1:0 includes element 1; posterior's max_t==0 branch therefore has
% tau=2 for split chains of length 3..5, not the antithetic upper ESS cap.
rng(816,'twister');
for n=[6 7 10]
    x=randn(n,4); d=dsc_mcmc_diagnostics(x);
    expected=4*floor(n/2);
    verifyTrue(testCase,d.available);
    verifyEqual(testCase,d.ess_bulk,expected,'AbsTol',1e-12);
    verifyEqual(testCase,d.ess_tail,expected,'AbsTol',1e-12);
    verifyEqual(testCase,d.mcse_mean,std(x(:))/sqrt(expected),'RelTol',1e-12);
end
end

function testSavedIidChunksPassAndWriteArtifacts(testCase)
[dirs,cfg,out]=fixture(testCase,4,1600,true,'-v7.3');
% Sampler chunks are already retained. An absurd burn-in must not trim them.
cfg.burnin=1e9; cfg.thin=7;
[report,details]=dsc_diagnose_chains(dirs,cfg,out);
verifyTrue(testCase,report.passed); verifyEqual(testCase,report.status,'passed');
verifyEqual(testCase,report.actual_draws_per_chain,[1600 1600 1600 1600]);
verifyEqual(testCase,report.draws_per_chain,1600);
verifyEqual(testCase,report.chain_ids,[11 12 13 14]);
verifyEqual(testCase,report.retained_iteration_ranges,repmat([102 3300],4,1));
verifyEqual(testCase,report.quantities_checked,16); verifyEqual(testCase,height(details),16);
verifyEqual(testCase,nnz(details.family=="V"),3); % includes the off-diagonal
verifyEqual(testCase,nnz(details.family=="P_pairs"),2);
verifyEqual(testCase,nnz(details.date=="2020-01-03"),5);
verifyEqual(testCase,nnz(details.date=="2020-01-10"),5);
verifyTrue(testCase,all(details.passed));
verifyTrue(testCase,isfile(fullfile(out,'convergence_diagnostics.json')));
verifyTrue(testCase,isfile(fullfile(out,'convergence_diagnostics.csv')));
saved=load(fullfile(out,'convergence_diagnostics.mat'),'report','details');
verifyEqual(testCase,saved.report,report); verifyEqual(testCase,saved.details,details);
decoded=jsondecode(fileread(fullfile(out,'convergence_diagnostics.json')));
verifyTrue(testCase,decoded.passed);
% Temporary staging must be cleaned; only the requested reports remain.
entries=dir(out); verifyEqual(testCase,nnz(~[entries.isdir]),3);
end

function testChunkNumericOrderAndEarliestAlignment(testCase)
[dirs,cfg,out]=fixture(testCase,4,[800 900 1000 1100],true,'-v7');
[report,details]=dsc_diagnose_chains(dirs,cfg,out);
verifyEqual(testCase,report.draws_per_chain,800);
verifyEqual(testCase,report.actual_draws_per_chain,[800 900 1000 1100]);
verifyFalse(testCase,report.passed); verifyEqual(testCase,report.status,'failed');
verifyTrue(testCase,contains(report.reason,'Unequal retained counts'));
% Fixture files have numeric suffixes 2 and 10, whose lexical order is wrong.
x=zeros(800,4);
for c=1:4
    a=load(fullfile(dirs{c},'posterior_chunk_2.mat'),'chunk');
    b=load(fullfile(dirs{c},'posterior_chunk_10.mat'),'chunk');
    values=[a.chunk.sig2h;b.chunk.sig2h]; x(:,c)=values(1:800,1);
end
expected=dsc_mcmc_diagnostics(x); row=find(details.quantity=="sig2h(1)");
verifyEqual(testCase,details.rhat(row),expected.rhat,'AbsTol',1e-12);
verifyEqual(testCase,details.ess_bulk(row),expected.ess_bulk,'RelTol',1e-12);
end

function testFixedConstantsAreExcludedButAllStateDatesChecked(testCase)
[dirs,cfg,out]=fixture(testCase,4,1400,false,'-v7');
[report,details]=dsc_diagnose_chains(dirs,cfg,out);
verifyTrue(testCase,report.passed); verifyEqual(testCase,height(details),10);
verifyFalse(testCase,any(ismember(details.family,["V","sig2h","sig2r"])));
% Corrupt only the final date of P to ensure every state date is checked.
for c=1:4
    files=dir(fullfile(dirs{c},'posterior_chunk_*.mat'));
    for k=1:numel(files)
        file=fullfile(files(k).folder,files(k).name); s=load(file,'chunk'); chunk=s.chunk;
        chunk.P_pairs(end,1,:)=.5; save(file,'chunk','-v7');
    end
end
[report,details]=dsc_diagnose_chains(dirs,cfg,out);
verifyFalse(testCase,report.passed);
verifyTrue(testCase,ismember('P(2,1)[2020-01-10]',report.unavailable_quantities));
verifyFalse(testCase,details.available(details.quantity=="P(2,1)[2020-01-10]"));
end

function testOffDiagonalCovarianceMustPass(testCase)
[dirs,cfg,out]=fixture(testCase,4,900,true,'-v7');
files=dir(fullfile(dirs{4},'posterior_chunk_*.mat'));
for k=1:numel(files)
    file=fullfile(files(k).folder,files(k).name); s=load(file,'chunk'); chunk=s.chunk;
    chunk.V(2,1,:)=chunk.V(2,1,:)+.5; chunk.V(1,2,:)=chunk.V(2,1,:);
    save(file,'chunk','-v7');
end
[report,details]=dsc_diagnose_chains(dirs,cfg,out);
verifyFalse(testCase,report.passed); verifyTrue(testCase,ismember('V(2,1)',report.failing_quantities));
verifyGreaterThan(testCase,details.rhat(details.quantity=="V(2,1)"),1.1);
end

function testOneChainAndShortDrawsCannotPass(testCase)
[dirs,cfg,out]=fixture(testCase,1,1400,true,'-v7');
report=dsc_diagnose_chains(dirs,cfg,out);
verifyFalse(testCase,report.passed); verifyEqual(testCase,report.status,'insufficient_chains');
[dirs,cfg,out]=fixture(testCase,4,20,true,'-v7');
report=dsc_diagnose_chains(dirs,cfg,out);
verifyFalse(testCase,report.passed); verifyEqual(testCase,report.status,'insufficient_draws');
dirs{4}=fullfile(out,'not_launched');
report=dsc_diagnose_chains(dirs,cfg,out);
verifyEqual(testCase,report.draws_per_chain,0); verifyEqual(testCase,report.status,'insufficient_draws');
end

function testOffAndCyclingDrawsAreNotChecked(testCase)
out=tempname(testCase.TestData.root); cfg=struct('convergence_mode','off');
[report,details]=dsc_diagnose_chains({'nonexistent_one','nonexistent_two'},cfg,out);
verifyFalse(testCase,report.passed); verifyEqual(testCase,report.status,'not_checked');
verifyEqual(testCase,height(details),0); verifyFalse(testCase,isfield(report,'actual_draws_per_chain'));
cfg=struct('convergence_mode','report','estimate_parameters',false,'parameter_smoothing','draws');
report=dsc_diagnose_chains({'nonexistent'},cfg,out);
verifyFalse(testCase,report.passed); verifyEqual(testCase,report.status,'not_checked');
verifyTrue(testCase,contains(report.reason,'conditional transition kernel'));
end

function testRemovedConvergenceModeIsRejectedBeforeDiagnosticWrites(testCase)
out=tempname(testCase.TestData.root); cfg=struct('convergence_mode','required');
verifyError(testCase,@()dsc_diagnose_chains({'not_launched'},cfg,out), ...
    'MATLAB:unrecognizedStringChoice');
verifyFalse(testCase,isfolder(out));
end

function testRejectsDuplicateChainIdentities(testCase)
[dirs,cfg,out]=fixture(testCase,2,120,true,'-v7');
verifyError(testCase,@()dsc_diagnose_chains(dirs([1 1]),cfg,out),'dsc:DiagnosticIdentity');
end

function testExpiredDiagnosticBudgetDoesNotReadChunks(testCase)
out=tempname(testCase.TestData.root);
cfg=struct('diagnostic_max_seconds',0,'convergence_mode','report');
[report,details]=dsc_diagnose_chains({'not_launched'},cfg,out);
verifyEqual(testCase,report.status,'not_checked'); verifyFalse(testCase,report.passed);
verifyTrue(testCase,contains(report.reason,'time budget exhausted'));
verifyFalse(testCase,isfield(report,'actual_draws_per_chain'));
verifyEqual(testCase,height(details),0);
end

function testRejectsChangedDatesAndIterationOverlap(testCase)
[dirs,cfg,out]=fixture(testCase,2,120,true,'-v7');
file=fullfile(dirs{2},'posterior_chunk_10.mat'); s=load(file,'chunk'); original=s.chunk;
chunk=original; chunk.dates(end)=chunk.dates(end)+calweeks(1); save(file,'chunk','-v7');
verifyError(testCase,@()dsc_diagnose_chains(dirs,cfg,out),'dsc:DiagnosticIdentity');
chunk=original; chunk.iterations=chunk.iterations-120; save(file,'chunk','-v7');
verifyError(testCase,@()dsc_diagnose_chains(dirs,cfg,out),'dsc:DiagnosticIterations');
end

function testRejectsDifferentTargetsRepeatedSeedsAndMissingDraws(testCase)
[dirs,cfg,out]=fixture(testCase,2,120,true,'-v7');
path=fullfile(dirs{2},'checkpoint.mat'); loaded=load(path,'checkpoint'); original=loaded.checkpoint;
checkpoint=original; checkpoint.identity.priors.target_marker=2; save(path,'checkpoint','-v7');
verifyError(testCase,@()dsc_diagnose_chains(dirs,cfg,out),'dsc:DiagnosticIdentity');
checkpoint=original; checkpoint.identity.cfg.seed=checkpoint.identity.cfg.seed-1;
save(path,'checkpoint','-v7'); % distinct IDs now share the same effective seed
verifyError(testCase,@()dsc_diagnose_chains(dirs,cfg,out),'dsc:DiagnosticIdentity');
checkpoint=original; checkpoint.saved=checkpoint.saved+1; save(path,'checkpoint','-v7');
verifyError(testCase,@()dsc_diagnose_chains(dirs,cfg,out),'dsc:DiagnosticIterations');
delete(path);
verifyError(testCase,@()dsc_diagnose_chains(dirs,cfg,out),'dsc:DiagnosticIdentity');
end

function [dirs,cfg,out]=fixture(testCase,chains,counts,estimate,version)
rng(841,'twister'); root=tempname(testCase.TestData.root); mkdir(root);
if isscalar(counts), counts=repmat(counts,1,chains); end
dirs=cell(1,chains); out=fullfile(root,'diagnostics');
cfg=struct('convergence_mode','report','estimate_parameters',estimate, ...
    'rhat_threshold',1.01,'min_ess',400,'max_mcse_ratio',.05,'min_diagnostic_draws',100);
for c=1:chains
    dirs{c}=fullfile(root,sprintf('chain_%d',c)); mkdir(dirs{c});
    n=counts(c); B=randn(2,2,n); h=randn(2,2,n); P=tanh(.1*randn(2,1,n));
    V=zeros(2,2,n); V(1,1,:)=2+exp(.1*randn(1,1,n)); V(2,2,:)=2+exp(.1*randn(1,1,n));
    V(2,1,:)=.05*randn(1,1,n); V(1,2,:)=V(2,1,:);
    sig2h=exp(.1*randn(n,2)); sig2r=exp(.1*randn(n,1));
    if ~estimate, V=repmat(eye(2),1,1,n); sig2h(:)=1; sig2r(:)=1; end
    cuts={1:floor(n/2),floor(n/2)+1:n}; suffix=[2 10];
    for k=1:2
        ix=cuts{k}; chunk=struct('B',B(:,:,ix),'h',h(:,:,ix),'P_pairs',P(:,:,ix), ...
            'V',V(:,:,ix),'sig2h',sig2h(ix,:),'sig2r',sig2r(ix,:), ...
            'iterations',100+2*ix,'chain_id',10+c, ...
            'dates',[datetime(2020,1,3);datetime(2020,1,10)], ...
            'tickers',["A","B"],'pair_i',2,'pair_j',1);
        save(fullfile(dirs{c},sprintf('posterior_chunk_%d.mat',suffix(k))),'chunk',version);
    end
    retention=struct('seed',1000,'chain_id',10+c,'burnin',100,'thin',2, ...
        'estimate_parameters',estimate);
    if ~estimate
        retention.fixed_parameters=struct('V',eye(2),'sig2h',[1 1],'sig2r',1);
    end
    identity=struct('dates',chunk.dates,'tickers',chunk.tickers,'returns',zeros(2,2), ...
        'mask',true(2,2),'source_hash','synthetic-diagnostic-target','cfg',retention, ...
        'priors',struct('target_marker',1));
    checkpoint=struct('identity',identity,'saved',n,'completed',100+2*n);
    save(fullfile(dirs{c},'checkpoint.mat'),'checkpoint',version);
end
end

function ess=direct_reference_ess(original)
% Independent O(N^2) autocovariances check the FFT implementation against
% stan-dev/posterior's biased 1/N autocovariance and exact lag truncation.
n=floor(size(original,1)/2); x=[original(1:n,:) original(end-n+1:end,:)];
m=size(x,2); centered=x-mean(x,1); ac=zeros(n,1);
for lag=0:n-1, ac(lag+1)=sum(sum(centered(1:n-lag,:).*centered(1+lag:n,:)))/(n*m); end
W=ac(1)*n/(n-1); vp=ac(1)+var(mean(x,1)); rho=1-(W-ac)/vp; rho(1)=1;
% Keep the initial nonnegative pair sums; the final even lag participates
% alone in the improved antithetic estimator used by posterior.
retained=zeros(n,1); retained(1:2)=rho(1:2); lag=0;
while lag<n-5&&sum(rho(lag+(1:2)))>0
    lag=lag+2;
    if sum(rho(lag+(1:2)))>=0, retained(lag+(1:2))=rho(lag+(1:2)); end
end
if rho(lag+1)>0, retained(lag+1)=rho(lag+1); end
for at=2:2:lag-2
    previous=sum(retained(at-1:at)); current=sum(retained(at+1:at+2));
    if current>previous, retained(at+1:at+2)=previous/2; end
end
if lag==0, truncated=retained(1); else, truncated=sum(retained(1:lag)); end
tau=-1+2*truncated+retained(lag+1);
ess=n*m/max(tau,1/log10(n*m));
end

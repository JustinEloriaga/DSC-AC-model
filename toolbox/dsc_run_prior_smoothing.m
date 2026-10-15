function result=dsc_run_prior_smoothing(panel,priors,cfg,run_dir)
%DSC_RUN_PRIOR_SMOOTHING Equal-weight prior mixture of conditional smoothers.
% Every component keeps one independent prior theta fixed for all sweeps.
% Component chains have distinct IDs/RNG seeds, independent warm-up and the
% same retained quota. max_seconds is ONE wall budget for all components.
started=tic;
[cfg,K,C,R]=validate_config(cfg);
if cfg.resume
    error('dsc:PriorResumeUnsupported', ...
        ['Prior-mixture resume is not implemented. Existing priors and checkpoints remain ' ...
         'unchanged; use a new run for a new authorized budget.']);
end
if ~isnumeric(panel.returns)||size(panel.returns,2)<2||isempty(panel.dates)
    error('dsc:PriorSmoothingConfig','The modeled panel must contain at least two variables and one date.');
end
parameter_file=fullfile(run_dir,'prior_parameter_draws.mat');
if isfile(parameter_file)||isfile(fullfile(run_dir,'result.mat'))|| ...
        ~isempty(dir(fullfile(run_dir,'prior_component_*')))
    error('dsc:ExistingRun','The prior smoothing directory already contains a run; it is preserved.');
end
% All theta are committed before a conditional chain can consume RNG/work.
prior_seed=mod(double(cfg.seed)+104729,2^32);
draws=dsc_draw_prior_parameters(priors,size(panel.returns,2),K,prior_seed);
if ~isfolder(run_dir), mkdir(run_dir); end
prior_parameters=struct('schema_version',1,'parameter_source','prior', ...
    'parameter_draw_seed',prior_seed,'prior_draws',draws,'priors',priors, ...
    'created_at',char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")), ...
    'source_hash',panel.source_hash,'tickers',panel.tickers,'dates',panel.dates, ...
    'max_seconds',cfg.max_seconds,'conditional_chains_per_parameter',C, ...
    'requested_retained_draws_per_chain',R,'burnin',cfg.burnin,'thin',cfg.thin, ...
    'chain_ids',cfg.chain_id+(0:K*C-1));
atomic_save(parameter_file,'prior_parameters',prior_parameters);
chains=cell(1,K*C); dirs=cell(1,K*C); components=cell(1,K);
conditional=cell(1,K); component_elapsed=zeros(1,K);
for d=1:K
    component_dir=fullfile(run_dir,sprintf('prior_component_%03d',d));
    conditional{d}=struct('status','not_started','passed',false, ...
        'reason','This fixed-parameter conditional component has not started.');
    for c=1:C
        k=(d-1)*C+c;
        if C==1, dirs{k}=component_dir;
        else, dirs{k}=fullfile(component_dir,'chains',sprintf('chain_%03d',c)); end
        chains{k}=empty_chain(panel,priors,cfg,dirs{k},cfg.chain_id+k-1);
    end
end
status='running'; reason='Running the fixed-parameter conditional components.'; failure='';
result=assemble(); persist();
for d=1:K
    if remaining_seconds()<=0
        status='time_limit'; reason='The shared prior-mixture wall-clock budget was reached.'; break
    end
    local=cfg; local.estimate_parameters=false; local.parameter_smoothing='mean';
    local.fixed_parameters=struct('V',draws.V(:,:,d), ...
        'sig2h',draws.sig2h(d,:),'sig2r',draws.sig2r(d,:));
    if isfield(local,'fixed_parameter_draws'), local=rmfield(local,'fixed_parameter_draws'); end
    local.chain_id=cfg.chain_id+(d-1)*C;
    local.max_seconds=max(realmin,remaining_seconds());
    component_dir=fullfile(run_dir,sprintf('prior_component_%03d',d));
    component_clock=tic;
    try
        current=dsc_run_chains(panel,priors,local,component_dir);
        components{d}=current;
        component_elapsed(d)=toc(component_clock);
        conditional{d}=current.convergence;
        for c=1:C, chains{(d-1)*C+c}=current.chain_results{c}; end
    catch problem
        component_elapsed(d)=toc(component_clock);
        % Recover already committed work without advancing or redrawing it.
        for c=1:C
            k=(d-1)*C+c; chains{k}=recover_chain(chains{k},dirs{k});
        end
        status='failed'; reason=problem.message; failure=problem.identifier;
        conditional{d}=struct('status','failed','passed',false,'reason',problem.message);
        result=assemble(); persist(); break
    end
    if strcmp(current.status,'failed')
        status='failed'; reason=current.reason; failure=current.failure_identifier;
    elseif any(current.saved_draws_per_chain~=R)|| ...
            any(current.completed_iterations_per_chain~=cfg.max_iterations)
        status='time_limit';
        reason='A conditional component did not meet its equal retained quota within the shared budget.';
    end
    result=assemble(); persist();
    if ~strcmp(status,'running'), break; end
end
if strcmp(status,'running')
    status='iteration_limit'; reason='Every conditional chain completed the requested equal retained quota.';
end
result=assemble(); persist();

    function value=remaining_seconds()
        value=cfg.max_seconds-toc(started);
    end
    function combined=assemble()
        combined=aggregate(chains,panel,priors,cfg,run_dir,dirs);
        for component=components
            if ~isempty(component{1})
                combined.diagnostic_seconds=combined.diagnostic_seconds+component{1}.diagnostic_seconds;
            end
        end
        counts=combined.saved_draws_per_chain;
        completed=reshape(counts==R & combined.completed_iterations_per_chain==cfg.max_iterations,C,K);
        completed_components=all(completed,1);
        all_passed=all(cellfun(@(report)isfield(report,'passed')&&logical(report.passed),conditional));
        mixture_complete=all(completed_components)&&~strcmp(status,'failed');
        combined.status=status; combined.reason=reason; combined.failure_identifier=failure;
        combined.inference_mode='prior_parameter_smoothing';
        combined.estimate_parameters=false; combined.parameter_smoothing='prior_draws';
        combined.parameter_source='prior'; combined.parameter_estimate=[];
        combined.prior_sampling=struct('requested_parameter_draws',K, ...
            'completed_parameter_draws',sum(completed_components), ...
            'conditional_chains_per_parameter',C, ...
            'parameter_draw_ids_per_chain',repelem(draws.draw_ids,C), ...
            'mixture_weights',draws.weights,'mixture_complete',mixture_complete, ...
            'conditional_convergence',{conditional}, ...
            'all_conditional_convergence_passed',all_passed, ...
            'parameter_draw_seed',prior_seed,'parameter_file',parameter_file, ...
            'requested_retained_draws_per_chain',R, ...
            'conditional_component_elapsed_seconds',component_elapsed, ...
            'elapsed_component_wall_seconds',sum(component_elapsed), ...
            'actual_retained_draws_per_chain',counts, ...
            'target','Integral p(states | observed returns, theta) p(theta) dtheta.', ...
            'weighting','Each saved prior theta has weight 1/K; no likelihood reweighting or redraw.', ...
            'outer_prior_integration_convergence','not_established');
        combined.convergence=struct('status','not_checked','passed',false, ...
            'reason',['Different theta are different targets; see conditional diagnostics; ' ...
                      'outer prior integration is not established.']);
        combined.convergence_established=false; combined.convergence_mode=cfg.convergence_mode;
        combined.convergence_details=struct();
        combined.auto_extend=false; combined.requested_retained_draws=R;
        combined.requested_retained_draws_per_chain=R;
        combined.max_retained_draws=R;
        combined.session_seconds=toc(started); combined.total_seconds=toc(started);
        combined.stage_elapsed_seconds=combined.total_seconds;
        combined.stage_seconds.elapsed=combined.total_seconds;
        combined.stage_seconds.diagnostics=combined.diagnostic_seconds;
        combined.paths.parameters=parameter_file;
        combined.paths.prior_parameters=parameter_file;
        combined.paths.prior_sampling_metadata=fullfile(run_dir,'prior_sampling_metadata.json');
        combined.paths.summary=fullfile(run_dir,'pilot_summary.json');
        combined.paths.result=fullfile(run_dir,'result.mat');
    end
    function persist()
        atomic_json(result.paths.prior_sampling_metadata,result.prior_sampling);
        exclude={'chain_results','priors','dates','tickers','pair_i','pair_j','convergence_details'};
        atomic_json(result.paths.summary,rmfield(result,exclude));
        atomic_save(result.paths.result,'result',result);
    end
end

function [cfg,K,C,R]=validate_config(cfg)
defaults=struct('prior_parameter_draws',10,'num_chains',1,'chain_id',1, ...
    'resume',false,'auto_extend',false,'parallel_chains',false,'convergence_mode','report');
for name=fieldnames(defaults)'
    if ~isfield(cfg,name{1}), cfg.(name{1})=defaults.(name{1}); end
end
for name={'prior_parameter_draws','num_chains','chain_id','max_iterations','thin','chunk_size'}
    if ~isfield(cfg,name{1})||~valid_integer(cfg.(name{1}),1)
        error('dsc:PriorSmoothingConfig','%s must be a positive integer.',name{1});
    end
end
for name={'burnin','seed'}
    if ~isfield(cfg,name{1})||~valid_integer(cfg.(name{1}),0)
        error('dsc:PriorSmoothingConfig','%s must be a nonnegative integer.',name{1});
    end
end
for name={'resume','auto_extend','parallel_chains'}
    value=cfg.(name{1});
    if ~(isnumeric(value)||islogical(value))||~isreal(value)||~isscalar(value)||~ismember(value,[0 1])
        error('dsc:PriorSmoothingConfig','%s must be true or false.',name{1});
    end
    cfg.(name{1})=logical(value);
end
if cfg.auto_extend
    error('dsc:PriorSmoothingConfig', ...
        'AutoExtend must be false: every prior theta must receive the same conditional retained quota.');
end
if ~isfield(cfg,'max_seconds')||~isnumeric(cfg.max_seconds)||~isreal(cfg.max_seconds)|| ...
        ~isscalar(cfg.max_seconds)||~isfinite(cfg.max_seconds)||cfg.max_seconds<=0
    error('dsc:PriorSmoothingConfig','max_seconds must be a positive finite shared wall budget.');
end
cfg.convergence_mode=validatestring(cfg.convergence_mode,{'off','report'});
K=cfg.prior_parameter_draws; C=cfg.num_chains;
if K*C>flintmax||cfg.seed+cfg.chain_id+K*C-2>2^32-1
    error('dsc:PriorSmoothingConfig','The full collection needs distinct valid chain IDs and RNG seeds.');
end
if cfg.max_iterations<=cfg.burnin||mod(cfg.max_iterations-cfg.burnin,cfg.thin)~=0
    error('dsc:PriorSmoothingConfig', ...
        'max_iterations must equal burnin + retained_draws*thin with at least one retained draw.');
end
R=(cfg.max_iterations-cfg.burnin)/cfg.thin;
end

function tf=valid_integer(value,minimum)
tf=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&& ...
    value>=minimum&&value==floor(value);
end

function result=empty_chain(panel,priors,cfg,folder,id)
[T,m]=size(panel.returns); q=m*(m-1)/2;
timing=struct('mean',0,'mean_variance',0,'correlation',0,'correlation_variance',0, ...
    'volatility',0,'volatility_variance',0,'checkpoint',0);
result=struct('status','not_started','reason','Conditional component has not started.', ...
    'failure_identifier','','chain_id',id,'seed',cfg.seed+id-1, ...
    'estimate_parameters',false,'inference_mode','fixed_parameter_smoothing', ...
    'parameter_smoothing','mean','parameter_draws',0,'parameter_estimate',[], ...
    'completed_iterations',0,'burnin',cfg.burnin,'saved_draws',0, ...
    'session_completed_iterations',0,'unfinished_stage','','unfinished_sweep_seconds',0, ...
    'T',T,'m',m,'pairs',q,'first_date',char(string(panel.dates(1),'yyyy-MM-dd')), ...
    'last_date',char(string(panel.dates(end),'yyyy-MM-dd')), ...
    'session_seconds',0,'total_seconds',0,'stage_seconds',timing, ...
    'r_proposals',0,'h_proposals',0,'r_invalid_proposals',0,'h_invalid_proposals',0, ...
    'sweep_seconds',zeros(0,1),'r_proposals_by_sweep',zeros(0,1), ...
    'h_proposals_by_sweep',zeros(0,1),'run_dir',folder,'priors',priors, ...
    'dates',panel.dates,'tickers',panel.tickers,'pair_i',panel.pair_i,'pair_j',panel.pair_j, ...
    'paths',struct());
end

function chain=recover_chain(chain,folder)
result_file=fullfile(folder,'result.mat');
if isfile(result_file)
    saved=load(result_file,'result');
    if isfield(saved.result,'chain_results'), chain=saved.result.chain_results{1};
    else, chain=saved.result; end
end
checkpoint_file=fullfile(folder,'checkpoint.mat');
if ~isfile(checkpoint_file), return; end
saved=load(checkpoint_file,'checkpoint'); cp=saved.checkpoint;
chain.completed_iterations=cp.completed; chain.saved_draws=cp.saved;
chain.parameter_draws=cp.parameter_sums.count;
chain.total_seconds=cp.total_seconds; chain.stage_seconds=cp.timing;
chain.sweep_seconds=cp.diagnostics.sweep_seconds;
chain.r_proposals_by_sweep=cp.diagnostics.r_proposals;
chain.h_proposals_by_sweep=cp.diagnostics.h_proposals;
chain.r_proposals=sum(cp.diagnostics.r_proposals); chain.h_proposals=sum(cp.diagnostics.h_proposals);
chain.r_invalid_proposals=sum(cp.diagnostics.r_invalid_proposals);
chain.h_invalid_proposals=sum(cp.diagnostics.h_invalid_proposals);
chain.session_completed_iterations=cp.completed;
end

function result=aggregate(chains,panel,priors,cfg,folder,dirs)
result=chains{1};
result=rmfield(result,{'chain_id','seed'});
result.chain_results=chains; result.chain_dirs=dirs;
result.chain_ids=cellfun(@(r)r.chain_id,chains); result.seeds=cellfun(@(r)r.seed,chains);
result.saved_draws_per_chain=cellfun(@(r)r.saved_draws,chains);
result.completed_iterations_per_chain=cellfun(@(r)r.completed_iterations,chains);
for name={'saved_draws','parameter_draws','completed_iterations','r_proposals', ...
        'h_proposals','r_invalid_proposals','h_invalid_proposals','session_completed_iterations'}
    result.(name{1})=sum(cellfun(@(r)r.(name{1}),chains));
end
names=fieldnames(chains{1}.stage_seconds);
for k=1:numel(names)
    name=names{k}; result.stage_seconds.(name)=sum(cellfun(@(r)r.stage_seconds.(name),chains));
end
for name={'sweep_seconds','r_proposals_by_sweep','h_proposals_by_sweep'}
    values=cellfun(@(r)r.(name{1}),chains,'UniformOutput',false);
    result.(name{1})=vertcat(values{:});
end
result.mean_sweep_seconds=mean(result.sweep_seconds);
result.median_sweep_seconds=median(result.sweep_seconds);
result.num_chains=numel(chains); result.parallel_chains=cfg.parallel_chains;
result.run_dir=folder; result.priors=priors; result.dates=panel.dates;
result.tickers=panel.tickers; result.pair_i=panel.pair_i; result.pair_j=panel.pair_j;
result.paths=struct('posterior_dirs',{dirs}); result.diagnostic_seconds=0;
% Component walls are measured separately by the caller. Sampler timers above
% are work totals across chains and must never be labeled elapsed wall time.
end

function atomic_save(path,name,value)
temporary=[tempname(fileparts(path)),'.mat'];
container=struct(); container.(name)=value;
save(temporary,'-struct','container','-v7.3'); movefile(temporary,path,'f');
end

function atomic_json(path,value)
temporary=[tempname(fileparts(path)),'.json'];
fid=fopen(temporary,'w'); assert(fid>=0,'dsc:OutputFile','Cannot write %s.',path);
guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true)); clear guard
movefile(temporary,path,'f');
end

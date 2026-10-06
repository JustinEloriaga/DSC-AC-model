function result=dsc_run_chains(panel,priors,cfg,run_dir)
%DSC_RUN_CHAINS Independent, resumable chains with bounded convergence checks.
% Retained targets and extension limits are PER CHAIN. max_seconds covers the
% entire stage, including sequential chains, batches and diagnostic work.
started=tic;
cfg=validate_config(cfg);
% Do this before directories, checkpoints, estimation or pool creation.
if cfg.parallel_chains
    [available,reason]=dsc_parallel_available();
    if ~available, error('dsc:ParallelUnavailable','%s',reason); end
end
chain_ids=cfg.chain_id+(0:cfg.num_chains-1);
chain_dirs=cell(1,cfg.num_chains);
for c=1:cfg.num_chains
    if cfg.num_chains==1, chain_dirs{c}=run_dir;
    else, chain_dirs{c}=fullfile(run_dir,'chains',sprintf('chain_%03d',c)); end
end
previous_seconds=0; previous_diagnostic_seconds=0;
root_result=fullfile(run_dir,'result.mat');
if cfg.resume&&isfile(root_result)
    loaded=load(root_result,'result'); prior_result=loaded.result;
    if isfield(prior_result,'chain_ids')&&~isequal(prior_result.chain_ids,chain_ids)
        error('dsc:ResumeMismatch','The resumed stage must preserve chain identities.');
    end
    if isfield(prior_result,'stage_elapsed_seconds')
        previous_seconds=prior_result.stage_elapsed_seconds;
    end
    if isfield(prior_result,'diagnostic_seconds')
        previous_diagnostic_seconds=prior_result.diagnostic_seconds;
    end
elseif ~cfg.resume&&isfile(root_result)
    error('dsc:ExistingRun','Run directory contains a result; explicitly resume it.');
end
sample_cfg=sampler_config(cfg);
chain_results=cell(1,cfg.num_chains);
checkpoint_seconds=zeros(1,cfg.num_chains);
for c=1:cfg.num_chains
    local_cfg=sample_cfg; local_cfg.chain_id=chain_ids(c);
    chain_results{c}=empty_result(panel,priors,local_cfg,chain_dirs{c});
    file=fullfile(chain_dirs{c},'result.mat');
    checkpoint=fullfile(chain_dirs{c},'checkpoint.mat');
    if ~cfg.resume&&(isfile(checkpoint)||~isempty(dir(fullfile(chain_dirs{c},'posterior_chunk_*.mat'))))
        error('dsc:ExistingRun','Chain %d already contains sampler output; explicitly resume it.',chain_ids(c));
    end
    if cfg.resume&&isfile(checkpoint)
        local_cfg.resume=true; local_cfg.resume_validate_only=true;
        validated=dsc_sample(panel,priors,local_cfg,chain_dirs{c});
        checkpoint_seconds(c)=validated.total_seconds;
    end
    if cfg.resume&&isfile(file)
        loaded=load(file,'result'); chain_results{c}=loaded.result;
        % The single-chain root also holds the coordinator's result.
        if isfield(chain_results{c},'chain_results')
            chain_results{c}=chain_results{c}.chain_results{1};
        end
    end
    if cfg.resume&&isfile(checkpoint)&&(~isfile(file)|| ...
            chain_results{c}.completed_iterations~=validated.completed_iterations|| ...
            chain_results{c}.saved_draws~=validated.saved_draws)
        % result.mat is a batch summary; a committed checkpoint can be newer
        % after interruption. Finalize its retained buffer without advancing
        % the chain, then use that authoritative state for targets and pooling.
        local_cfg.resume_validate_only=false; local_cfg.resume_finalize_only=true;
        chain_results{c}=dsc_sample(panel,priors,local_cfg,chain_dirs{c});
    end
end
if cfg.resume
    if cfg.parallel_chains, sampler_elapsed=max(checkpoint_seconds);
    else, sampler_elapsed=sum(checkpoint_seconds); end
    % A crash can precede the stage summary; never replenish the time already
    % recorded in chain checkpoints merely because that summary is missing.
    previous_seconds=max(previous_seconds,sampler_elapsed+previous_diagnostic_seconds);
end
if ~exist(run_dir,'dir'), mkdir(run_dir); end
pool=[];
if cfg.parallel_chains
    pool=gcp('nocreate');
    if ~isempty(pool)&&isa(pool,'parallel.ThreadPool')
        error('dsc:ProcessPoolRequired', ...
            'Parallel chains require a process pool; close the existing thread pool first.');
    end
    if isempty(pool), pool=parpool('Processes',cfg.num_chains); end
    if pool.NumWorkers<cfg.num_chains
        error('dsc:ParallelWorkers','The process pool needs at least NumChains=%d workers.',cfg.num_chains);
    end
end
requested=max(0,floor((cfg.max_iterations-cfg.burnin)/cfg.thin));
initial_iterations=cfg.max_iterations;
counts=cellfun(@(r) r.saved_draws,chain_results);
initial_completed=cellfun(@(r) r.completed_iterations,chain_results);
target=min(requested,max(cfg.check_every,min(counts)+cfg.check_every));
diagnostic_seconds=previous_diagnostic_seconds; report=struct('passed',false,'status','not_checked');
details=struct(); stopped='iteration_limit'; failure=''; reason='';
last_checked_counts=[];
while true
    iteration_target=cfg.burnin+target*cfg.thin;
    if target==requested, iteration_target=initial_iterations; end
    if requested==0, iteration_target=initial_iterations; end
    iteration_target=max(1,iteration_target);
    active=find(cellfun(@(r) r.completed_iterations<iteration_target,chain_results));
    if remaining_seconds()<=0
        stopped='time_limit'; reason='The cumulative stage wall-clock limit was reached.'; break
    end
    if cfg.parallel_chains&&~isempty(active)
        futures=cell(1,numel(active));
        for k=1:numel(active)
            c=active(k); local_cfg=batch_config(c,iteration_target);
            futures{k}=parfeval(pool,@dsc_sample,1,panel,priors,local_cfg,chain_dirs{c});
        end
        try
            for k=1:numel(active)
                chain_results{active(k)}=fetchOutputs(futures{k});
            end
        catch problem
            for k=1:numel(futures), cancel(futures{k}); end
            rethrow(problem)
        end
    else
        for c=active
            if remaining_seconds()<=0, break; end
            local_cfg=batch_config(c,iteration_target);
            chain_results{c}=dsc_sample(panel,priors,local_cfg,chain_dirs{c});
            if any(strcmp(chain_results{c}.status,{'failed','time_limit'})), break; end
        end
    end
    counts=cellfun(@(r) r.saved_draws,chain_results);
    completed=cellfun(@(r) r.completed_iterations,chain_results);
    failed=find(cellfun(@(r) strcmp(r.status,'failed'),chain_results),1);
    if ~isempty(failed)
        stopped='failed'; reason=chain_results{failed}.reason;
        failure=chain_results{failed}.failure_identifier; break
    end
    if all(counts>=cfg.min_diagnostic_draws)&&remaining_seconds()>0
        evaluate_diagnostics();
    end
    if any(completed<iteration_target)||remaining_seconds()<=0
        stopped='time_limit'; reason='The cumulative stage wall-clock limit was reached.'; break
    end
    % A favorable early check never shortens the user's initial draw target.
    if target<requested
        target=min(requested,target+cfg.check_every); continue
    end
    if report.passed
        stopped='converged'; reason='All configured convergence criteria passed.'; break
    end
    if ~cfg.auto_extend
        reason='Requested iterations completed; convergence criteria have not been established.'; break
    end
    if target>=cfg.max_retained_draws
        stopped='retained_draw_limit'; reason='The per-chain retained-draw cap was reached before convergence.'; break
    end
    target=min(cfg.max_retained_draws,target+cfg.check_every);
end
% Persist a diagnostic report even when zero draws or a wall cap prevented a
% regular check. Diagnostic work cannot launch another batch after the cap.
counts=cellfun(@(r) r.saved_draws,chain_results);
if ~isequal(last_checked_counts,counts), evaluate_diagnostics(); end
if isempty(reason), reason='Requested iterations completed.'; end
result=aggregate_results(chain_results,panel,priors,run_dir,chain_dirs,chain_ids);
result.session_completed_iterations=result.completed_iterations-sum(initial_completed);
result.status=stopped; result.reason=reason; result.failure_identifier=failure;
result.convergence=report; result.convergence_established=logical(report.passed);
result.convergence_details=details;
result.num_chains=cfg.num_chains; result.parallel_chains=cfg.parallel_chains;
result.convergence_mode=cfg.convergence_mode; result.auto_extend=cfg.auto_extend;
result.requested_retained_draws=requested;
result.requested_retained_draws_per_chain=requested;
result.max_retained_draws=cfg.max_retained_draws;
result.session_seconds=toc(started);
result.stage_elapsed_seconds=previous_seconds+toc(started);
result.total_seconds=result.stage_elapsed_seconds;
result.diagnostic_seconds=diagnostic_seconds;
result.stage_seconds.elapsed=result.stage_elapsed_seconds;
result.stage_seconds.diagnostics=diagnostic_seconds;
result.paths.summary=fullfile(run_dir,'pilot_summary.json');
result.paths.result=root_result;
result.paths.convergence=fullfile(run_dir,'convergence','convergence_diagnostics.json');
summary=rmfield(result,{'chain_results','priors','dates','tickers','pair_i','pair_j','convergence_details'});
atomic_json(result.paths.summary,summary);
atomic_save(root_result,'result',result);

    function seconds=remaining_seconds()
        seconds=cfg.max_seconds-previous_seconds-toc(started);
    end
    function local_cfg=batch_config(c,limit)
        local_cfg=sample_cfg; local_cfg.chain_id=chain_ids(c);
        local_cfg.max_iterations=limit;
        local_cfg.max_seconds=max(realmin,remaining_seconds());
        local_cfg.resume=isfile(fullfile(chain_dirs{c},'checkpoint.mat'));
    end
    function evaluate_diagnostics()
        diagnostic_clock=tic;
        diagnostic_cfg=cfg;
        diagnostic_cfg.diagnostic_max_seconds=max(0,remaining_seconds());
        [report,details]=dsc_diagnose_chains(chain_dirs,diagnostic_cfg,run_dir);
        diagnostic_seconds=diagnostic_seconds+toc(diagnostic_clock);
        last_checked_counts=cellfun(@(r) r.saved_draws,chain_results);
    end
end

function cfg=validate_config(cfg)
defaults=struct('num_chains',1,'parallel_chains',false,'convergence_mode','report', ...
    'auto_extend',false,'check_every',250,'max_retained_draws',2000, ...
    'rhat_threshold',1.01,'min_ess',400,'max_mcse_ratio',.05, ...
    'min_diagnostic_draws',100,'chain_id',1,'resume',false);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(cfg,names{k}), cfg.(names{k})=defaults.(names{k}); end
end
if isfield(cfg,'chunk_size') && isfinite(cfg.chunk_size) && cfg.chunk_size>0
    cfg.check_every=cfg.chunk_size;
    cfg.min_diagnostic_draws=cfg.chunk_size;
end
for name={'num_chains','check_every','max_retained_draws','min_diagnostic_draws','chain_id'}
    value=cfg.(name{1});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<1||value~=floor(value)
        error('dsc:ChainConfig','%s must be a positive integer.',name{1});
    end
end
for name={'parallel_chains','auto_extend','resume'}
    value=cfg.(name{1});
    if ~(islogical(value)||isnumeric(value))||~isreal(value)||~isscalar(value)||~ismember(value,[0 1])
        error('dsc:ChainConfig','%s must be true or false.',name{1});
    end
    cfg.(name{1})=logical(value);
end
cfg.convergence_mode=validatestring(cfg.convergence_mode,{'off','report'});
for name={'rhat_threshold','min_ess','max_mcse_ratio','max_seconds'}
    value=cfg.(name{1});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<=0
        error('dsc:ChainConfig','%s must be positive and finite.',name{1});
    end
end
for name={'max_iterations','burnin','thin','seed'}
    value=cfg.(name{1});
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<0||value~=floor(value)
        error('dsc:ChainConfig','%s must be a nonnegative integer.',name{1});
    end
end
if cfg.max_iterations<1||cfg.thin<1||cfg.rhat_threshold<=1|| ...
        cfg.min_diagnostic_draws<6||cfg.seed+cfg.chain_id+cfg.num_chains-2>2^32-1
    error('dsc:ChainConfig','Invalid iteration, thinning, R-hat threshold or chain seed.');
end
requested=max(0,floor((cfg.max_iterations-cfg.burnin)/cfg.thin));
if cfg.auto_extend&&(cfg.max_retained_draws<requested||strcmp(cfg.convergence_mode,'off'))
    error('dsc:ChainConfig','AutoExtend requires diagnostics and a cap at least as large as the requested per-chain target.');
end
end

function cfg=sampler_config(cfg)
% Scheduling and diagnostic choices may change on resume without changing
% a chain's target distribution, warm-up or random sequence.
names={'num_chains','parallel_chains','convergence_mode','auto_extend', ...
    'check_every','max_retained_draws','rhat_threshold','min_ess', ...
    'max_mcse_ratio','min_diagnostic_draws'};
dispersed=cfg.num_chains>1;
cfg=rmfield(cfg,names);
cfg.randomize_initial_state=dispersed;
end

function result=empty_result(panel,priors,cfg,folder)
[T,m]=size(panel.returns); nr=m*(m-1)/2;
estimate=true; if isfield(cfg,'estimate_parameters'), estimate=cfg.estimate_parameters; end
mode='parameter_estimation'; if ~estimate, mode='fixed_parameter_smoothing'; end
smoothing='mean'; if isfield(cfg,'parameter_smoothing')&&~isempty(cfg.parameter_smoothing), smoothing=cfg.parameter_smoothing; end
timing=struct('mean',0,'mean_variance',0,'correlation',0,'correlation_variance',0, ...
    'volatility',0,'volatility_variance',0,'checkpoint',0);
result=struct('status','not_started','reason','The stage ended before this chain started.', ...
    'failure_identifier','','convergence_established',false,'chain_id',cfg.chain_id, ...
    'seed',cfg.seed+cfg.chain_id-1,'randomize_initial_state',cfg.randomize_initial_state, ...
    'estimate_parameters',logical(estimate),'inference_mode',mode,'parameter_smoothing',smoothing, ...
    'parameter_draws',0,'parameter_estimate',[],'completed_iterations',0,'burnin',cfg.burnin, ...
    'saved_draws',0,'session_completed_iterations',0,'unfinished_stage','', ...
    'unfinished_sweep_seconds',0,'T',T,'m',m,'pairs',nr, ...
    'first_date',char(string(panel.dates(1),'yyyy-MM-dd')), ...
    'last_date',char(string(panel.dates(end),'yyyy-MM-dd')), ...
    'session_seconds',0,'total_seconds',0,'stage_seconds',timing, ...
    'r_proposals',0,'h_proposals',0,'r_invalid_proposals',0,'h_invalid_proposals',0, ...
    'sweep_seconds',zeros(0,1),'r_proposals_by_sweep',zeros(0,1), ...
    'h_proposals_by_sweep',zeros(0,1),'run_dir',folder,'priors',priors, ...
    'dates',panel.dates,'tickers',panel.tickers,'pair_i',panel.pair_i,'pair_j',panel.pair_j, ...
    'paths',struct());
end

function result=aggregate_results(chains,panel,priors,folder,dirs,ids)
result=chains{1};
result.chain_results=chains; result.chain_dirs=dirs; result.chain_ids=ids;
result.saved_draws_per_chain=cellfun(@(r) r.saved_draws,chains);
result.completed_iterations_per_chain=cellfun(@(r) r.completed_iterations,chains);
result.seeds=cellfun(@(r) r.seed,chains);
for name={'saved_draws','parameter_draws','completed_iterations','r_proposals', ...
        'h_proposals','r_invalid_proposals','h_invalid_proposals','session_completed_iterations'}
    result.(name{1})=sum(cellfun(@(r) r.(name{1}),chains));
end
names=fieldnames(chains{1}.stage_seconds);
for k=1:numel(names)
    result.stage_seconds.(names{k})=sum(cellfun(@(r) r.stage_seconds.(names{k}),chains));
end
for name={'sweep_seconds','r_proposals_by_sweep','h_proposals_by_sweep'}
    values=cellfun(@(r) r.(name{1}),chains,'UniformOutput',false);
    result.(name{1})=vertcat(values{:});
end
result.mean_sweep_seconds=mean(result.sweep_seconds);
result.median_sweep_seconds=median(result.sweep_seconds);
result.parameter_estimate=[];
if result.parameter_draws>0
    m=size(panel.returns,2); nr=m*(m-1)/2;
    estimate=struct('V',zeros(m),'sig2h',zeros(1,m),'sig2r',zeros(1,nr));
    for c=1:numel(chains)
        n=chains{c}.parameter_draws;
        if n==0, continue; end
        for name={'V','sig2h','sig2r'}
            estimate.(name{1})=estimate.(name{1})+n*chains{c}.parameter_estimate.(name{1})/result.parameter_draws;
        end
    end
    result.parameter_estimate=estimate;
end
result.run_dir=folder; result.priors=priors;
result.dates=panel.dates; result.tickers=panel.tickers;
result.pair_i=panel.pair_i; result.pair_j=panel.pair_j;
if numel(chains)>1
    result=rmfield(result,{'chain_id','seed'});
    result.paths=struct();
end
result.paths.posterior_dirs=dirs;
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

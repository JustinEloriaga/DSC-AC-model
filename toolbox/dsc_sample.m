function result = dsc_sample(panel,priors,cfg,run_dir)
% Constant-mean DSC smoother. All files belong to one chain.
% Observation innovations are conditionally independent through time given
% stochastic volatility/correlations; no time-varying mean or mean covariance.
% A time limit aborts an unfinished sweep; only completed sweeps commit.
started=tic;
Y=panel.returns; [T,m]=size(Y); nr=m*(m-1)/2;
if ~isequal(size(panel.observation_mask),size(Y))||m<2, error('dsc:PanelShape','Invalid panel/mask dimensions.'); end
mask=logical(panel.observation_mask)&isfinite(Y);
if numel(panel.dates)~=T||any(diff(panel.dates)<=seconds(0))
    error('dsc:PanelDates','Dates must increase strictly and match returns.');
end
if isfield(priors,'likelihood_start_date')&& ...
        (isempty(panel.dates)||panel.dates(1)~=priors.likelihood_start_date)
    error('dsc:CalibrationOverlap','The model sample must start after the reserved initial-prior calibration weeks.');
end
if any(~isfinite(Y(logical(panel.observation_mask))))
    error('dsc:ObservedNaN','Observed returns must be finite.');
end
if ~isequal(priors.tickers,panel.tickers)||~isequal(priors.source_hash,panel.source_hash)
    error('dsc:PriorIdentity','Priors must match the panel tickers and source hash.');
end
if ~isfield(cfg,'chain_id'), cfg.chain_id=1; end
if ~isfield(cfg,'resume'), cfg.resume=false; end
if ~isfield(cfg,'mean_model'), cfg.mean_model='constant'; end
if ~strcmp(string(cfg.mean_model),'constant')
    error('dsc:MeanModel','This sampler estimates a constant mean; historical dynamic-mean runs require the historical sampler.');
end
cfg.mean_model='constant';
schema='constant_mu_v1';
[mu0,mu_cov0]=mean_prior(priors,m);
[prior_mean_h,prior_scale_h]=innovation_prior(priors,'h',m);
[prior_mean_r,prior_scale_r]=innovation_prior(priors,'r',nr);
if ~isfield(cfg,'mean_backend'), cfg.mean_backend='auto'; end
cfg.mean_backend=validatestring(cfg.mean_backend,{'auto','native','matlab'});
if strcmp(cfg.mean_backend,'auto')
    if exist('dsc_constant_mean_mex','file')==3, cfg.mean_backend='native';
    else, cfg.mean_backend='matlab'; end
end
if strcmp(cfg.mean_backend,'native')&&exist('dsc_constant_mean_mex','file')~=3
    error('dsc:MissingKernel','Build dsc_constant_mean_mex before selecting the native mean backend.');
end
if ~isfield(cfg,'file_prefix'), cfg.file_prefix=''; end
assert((ischar(cfg.file_prefix)&&isrow(cfg.file_prefix))|| ...
    (isstring(cfg.file_prefix)&&isscalar(cfg.file_prefix))||isempty(cfg.file_prefix), ...
    'dsc:OutputPrefix','file_prefix must be a scalar filename prefix.');
cfg.file_prefix=char(cfg.file_prefix);
assert(isempty(cfg.file_prefix)||~isempty(regexp(cfg.file_prefix,'^[A-Za-z0-9_-]+$','once')), ...
    'dsc:OutputPrefix','file_prefix may contain only letters, numbers, underscores and hyphens.');
if ~isfield(cfg,'randomize_initial_state'), cfg.randomize_initial_state=false; end
if ~isscalar(cfg.randomize_initial_state)|| ...
        ~(islogical(cfg.randomize_initial_state)||isnumeric(cfg.randomize_initial_state))|| ...
        ~isreal(cfg.randomize_initial_state)||~ismember(cfg.randomize_initial_state,[0 1])
    error('dsc:SamplerConfig','randomize_initial_state must be true or false.');
end
cfg.randomize_initial_state=logical(cfg.randomize_initial_state);
if ~isfield(cfg,'estimate_parameters'), cfg.estimate_parameters=true; end
if ~isfield(cfg,'parameter_smoothing')||strlength(string(cfg.parameter_smoothing))==0
    cfg.parameter_smoothing='mean';
end
cfg.parameter_smoothing=validatestring(cfg.parameter_smoothing,{'draws','mean'});
if ~isscalar(cfg.estimate_parameters)|| ...
        ~(islogical(cfg.estimate_parameters)||isnumeric(cfg.estimate_parameters))|| ...
        ~isreal(cfg.estimate_parameters)||~ismember(cfg.estimate_parameters,[0 1])
    error('dsc:SamplerConfig','estimate_parameters must be true or false.');
end
cfg.estimate_parameters=logical(cfg.estimate_parameters);
if cfg.estimate_parameters
    inference_mode='parameter_estimation';
else
    inference_mode='fixed_parameter_smoothing';
    if strcmp(cfg.parameter_smoothing,'draws')
        if ~isfield(cfg,'fixed_parameter_draws')
            error('dsc:FixedParameters','Draw-based smoothing requires fixed_parameter_draws.sig2h and sig2r.');
        end
        cfg.fixed_parameter_draws=validate_parameter_draws(cfg.fixed_parameter_draws,m,nr);
        cfg.fixed_parameters=parameter_draw_at(cfg.fixed_parameter_draws,1);
    else
        if ~isfield(cfg,'fixed_parameters')
            error('dsc:FixedParameters','Mean-based smoothing requires fixed_parameters.sig2h and sig2r.');
        end
        cfg.fixed_parameters=validate_fixed_parameters(cfg.fixed_parameters,m,nr);
    end
end
if ~isfield(cfg,'correlation_backend'), cfg.correlation_backend='mex_new'; end
if ~isfield(cfg,'correlation_threads'), cfg.correlation_threads='auto'; end
% A direct sampler/pilot call runs one chain. Coordinated parallel chains
% arrive with their numeric allocation already resolved on the client.
cfg.correlation_threads=dsc_resolve_correlation_threads(cfg.correlation_threads,1,false);
if ~isfield(cfg,'correlation_warm_start'), cfg.correlation_warm_start=true; end
validateattributes(cfg.correlation_warm_start,{'logical','numeric'},{'scalar','binary'});
cfg.correlation_warm_start=logical(cfg.correlation_warm_start);
if ~isscalar(cfg.correlation_threads)||~isfinite(cfg.correlation_threads)|| ...
        cfg.correlation_threads<1||cfg.correlation_threads>64||cfg.correlation_threads~=floor(cfg.correlation_threads)
    error('dsc:SamplerConfig','Correlation thread count must be an integer from 1 to 64.');
end
cfg.correlation_backend=validatestring(cfg.correlation_backend,{'auto','matlab','mex','matlab_new','mex_new'});
if strcmp(cfg.correlation_backend,'auto')
    if exist('dsc_corr_batch_mex_new','file')==3
        cfg.correlation_backend='mex_new';
    elseif exist('dsc_corr_batch_mex','file')==3
        cfg.correlation_backend='mex';
    else
        cfg.correlation_backend='matlab';
    end
end
if strcmp(cfg.correlation_backend,'mex')&&exist('dsc_corr_batch_mex','file')~=3
    error('dsc:MissingKernel','Run build_dsc_mex before selecting the MEX backend.');
end
if strcmp(cfg.correlation_backend,'mex_new')&&exist('dsc_corr_batch_mex_new','file')~=3
    error('dsc:MissingKernel','Run build_dsc_mex before selecting the new MEX backend.');
end
implementation=dsc_sampler_signature(cfg.correlation_backend,cfg.mean_backend);
implementation.mean_model=cfg.mean_model;
if cfg.burnin<0||cfg.thin<1||cfg.chunk_size<1||cfg.max_iterations<1||cfg.max_seconds<=0
    error('dsc:SamplerConfig','Invalid iteration, thinning, chunk or time limit.');
end
if any([cfg.burnin cfg.thin cfg.chunk_size cfg.max_iterations]~=floor([cfg.burnin cfg.thin cfg.chunk_size cfg.max_iterations]))
    error('dsc:SamplerConfig','Iteration counts must be integers.');
end
if ~exist(run_dir,'dir')
    error('dsc:OutputDirectory','Sampler requires an existing output directory: %s',run_dir);
end
output_path=@(name) fullfile(run_dir,[cfg.file_prefix name]);
checkpoint_path=output_path('checkpoint.mat');
fixed_cfg=cfg;
for field={'max_iterations','max_seconds','resume','output_root','resume_validate_only','resume_finalize_only'}
    if isfield(fixed_cfg,field{1}), fixed_cfg=rmfield(fixed_cfg,field{1}); end
end
identity=struct('dates',panel.dates,'tickers',panel.tickers,'returns',Y, ...
    'mask',mask,'source_hash',panel.source_hash,'cfg',fixed_cfg,'priors',priors, ...
    'implementation',implementation,'schema',schema);
obs=cell(T,1); for t=1:T, obs{t}=find(mask(t,:)); end
sel=tril(true(m),-1);
[expected_i,expected_j]=find(sel);
if ~isequal(panel.pair_i(:),expected_i)||~isequal(panel.pair_j(:),expected_j)
    error('dsc:PairOrder','Pair indices must use column-major strict lower triangle order.');
end
timing=struct('mean',0,'correlation',0,'correlation_variance',0, ...
    'volatility',0,'volatility_variance',0,'checkpoint',0);
empty_buffer=struct('P_pairs',zeros(T,nr,0),'h',zeros(T,m,0), ...
    'mu',zeros(0,m),'sig2h',zeros(0,m), ...
    'sig2r',zeros(0,nr),'iterations',zeros(1,0));
if cfg.resume
    if ~exist(checkpoint_path,'file'), error('dsc:ResumeMissing','Checkpoint does not exist.'); end
    loaded=load(checkpoint_path,'checkpoint'); cp=loaded.checkpoint;
    if ~isfield(cp,'schema')||~strcmp(cp.schema,schema)|| ...
            isfield(cp.state,'B')||isfield(cp.state,'V')
        error('dsc:LegacyCheckpoint','This constant-mu sampler cannot resume historical B/V checkpoints.');
    end
    if ~isequaln(cp.identity,identity), error('dsc:ResumeMismatch','Panel, priors or fixed configuration differs.'); end
    % The coordinator must validate every resumed chain even if its existing
    % draws already meet this call's target or its stage budget has expired.
    % Return before restoring RNG or writing anything for this read-only path.
    if isfield(cfg,'resume_validate_only')&&cfg.resume_validate_only
        result=struct('completed_iterations',cp.completed,'saved_draws',cp.saved, ...
            'total_seconds',cp.total_seconds);
        return
    end
    state=cp.state; rng(cp.rng); completed=cp.completed; saved=cp.saved;
    next_chunk=cp.next_chunk; buffer=cp.buffer; performance=cp.performance;
    parameter_sums=cp.parameter_sums;
    timing=cp.timing; previous_seconds=cp.total_seconds;
else
    if exist(checkpoint_path,'file')||~isempty(dir(output_path('posterior_chunk_*.mat')))
        error('dsc:ExistingRun','Run directory already contains sampler output; explicitly resume it.');
    end
    rng(cfg.seed+cfg.chain_id-1,'twister');
    state.mu=mu0';
    state.h=repmat(priors.mh0',T,1);
    state.r=repmat(priors.mr0',T,1);
    if cfg.estimate_parameters
        state.sig2h=prior_mean_h;
        state.sig2r=prior_mean_r;
    else
        state.sig2h=cfg.fixed_parameters.sig2h;
        state.sig2r=cfg.fixed_parameters.sig2r;
        if isfield(cfg.fixed_parameters,'mu'), state.mu=cfg.fixed_parameters.mu; end
    end
    if cfg.randomize_initial_state
        state=dispersed_initial_state(state,priors,cfg.estimate_parameters,T,m,nr,mu0,mu_cov0, ...
            prior_scale_h,prior_scale_r);
    end
    completed=0; saved=0; next_chunk=1; buffer=empty_buffer;
    parameter_sums=struct('mu',zeros(1,m),'sig2h',zeros(1,m), ...
        'sig2r',zeros(1,nr),'count',0);
    performance=struct( ...
        'r_proposals',zeros(0,1),'h_proposals',zeros(0,1), ...
        'r_invalid_proposals',zeros(0,1),'h_invalid_proposals',zeros(0,1), ...
        'sweep_seconds',zeros(0,1));
    previous_seconds=0;
end
committed_rng=rng;
status='iteration_limit'; stop_reason='Configured iteration limit reached.';
failure_identifier='';
session_start_iterations=completed;
unfinished_stage=''; unfinished_sweep_seconds=0;
% Initial cache contains no random draws. Even a tiny budget receives a valid
% zero-sweep checkpoint and explicit time_limit status from the loop below.
cache=build_cache(state.r,obs,m,cfg.correlation_backend,mask,cfg.correlation_threads,@() []);
if cfg.randomize_initial_state&&~isfinite(cached_likelihood(Y-state.mu,state.h,cache,mask))
    error('dsc:InitialLikelihood','Dispersed initial state has a nonfinite log likelihood.');
end
if ~cfg.resume, commit_checkpoint(); end
while completed<cfg.max_iterations&&~(isfield(cfg,'resume_finalize_only')&&cfg.resume_finalize_only)
    before=state; rng_before=rng; timing_before=timing;
    iteration_clock=tic;
    stage_name='mean';
    try
        check_time(); stage=tic;
        if ~cfg.estimate_parameters && strcmp(cfg.parameter_smoothing,'draws')
            state=apply_parameter_draw(state,cfg.fixed_parameter_draws,completed+1);
        end
        zmean=randn(m,1);
        if strcmp(cfg.mean_backend,'native')
            mu=dsc_constant_mean_mex(Y,state.h,cache.P,mask,mu0,mu_cov0,zmean);
            check_time();
        else
            mu=draw_constant_mean_cached(Y,state.h,cache,mu0,mu_cov0,zmean,@check_time);
        end
        state.mu=mu';
        timing.mean=timing.mean+toc(stage);
        residual=Y-state.mu;
        stage_name='correlation'; stage=tic; ll=cached_likelihood(residual,state.h,cache,mask);
        % E and h stay fixed throughout all correlation-coordinate updates.
        % Standardize once, rather than once per date and proposal.
        Z=residual.*exp(-state.h/2); Z(~mask)=0;
        correlation_constant=-.5*(nnz(mask)*log(2*pi)+sum(state.h(mask)));
        r_proposals=0; r_invalid=0;
        for j=1:nr
            check_time();
            diagonal_guess=[];
            if cfg.correlation_warm_start, diagonal_guess=cache.diagonal; end
            f=@(x) correlation_likelihood(x,j,state.r,Z,correlation_constant, ...
                obs,m,cfg.correlation_backend,mask,cfg.correlation_threads,diagonal_guess,@check_time);
            [state.r(:,j),ll,ntry,ninvalid]=ellipse(state.r(:,j),priors.mr0(j), ...
                sqrt(state.sig2r(j)),priors.r0_scale,ll,f,@check_time);
            r_proposals=r_proposals+ntry; r_invalid=r_invalid+ninvalid;
        end
        cache_new=build_cache(state.r,obs,m,cfg.correlation_backend,mask,cfg.correlation_threads,@check_time);
        timing.correlation=timing.correlation+toc(stage);
        if cfg.estimate_parameters
            stage_name='correlation_variance'; stage=tic;
            state.sig2r=draw_variances(state.r,priors.mr0, ...
                priors.r0_scale,priors.ig_shape,prior_scale_r);
            timing.correlation_variance=timing.correlation_variance+toc(stage);
        end
        stage_name='volatility'; stage=tic; ll=cached_likelihood(residual,state.h,cache_new,mask); h_proposals=0; h_invalid=0;
        [directions,z,u]=volatility_cache(residual,state.h,cache_new,mask);
        for j=1:m
            check_time();
            direction=directions(:,j,:);
            cross=reshape(sum(direction.*u,1),T,1);
            precision=reshape(sum(direction.^2,1),T,1);
            old_h=state.h(:,j); old_z=z(:,j);
            f=@(x) dsc_volatility_likelihood(x,old_h,residual(:,j), ...
                old_z,cross,precision,mask(:,j),ll);
            [state.h(:,j),ll,ntry,ninvalid]=ellipse(state.h(:,j),priors.mh0(j), ...
                sqrt(state.sig2h(j)),priors.h0_scale,ll,f,@check_time);
            z(:,j)=residual(:,j).*exp(-state.h(:,j)/2); z(~mask(:,j),j)=0;
            u=u+direction.*reshape(z(:,j)-old_z,1,1,T);
            h_proposals=h_proposals+ntry; h_invalid=h_invalid+ninvalid;
        end
        timing.volatility=timing.volatility+toc(stage);
        if cfg.estimate_parameters
            stage_name='volatility_variance'; stage=tic;
            state.sig2h=draw_variances(state.h,priors.mh0, ...
                priors.h0_scale,priors.ig_shape,prior_scale_h);
            timing.volatility_variance=timing.volatility_variance+toc(stage);
        end
        check_time();
    catch problem
        unfinished_stage=stage_name; unfinished_sweep_seconds=toc(iteration_clock);
        state=before; rng(rng_before); timing=timing_before;
        failure_identifier=problem.identifier;
        if strcmp(problem.identifier,'dsc:TimeLimit')
            status='time_limit'; stop_reason='Time limit reached; unfinished sweep discarded.';
        else
            status='failed'; stop_reason=problem.message;
        end
        break
    end
    cache=cache_new; completed=completed+1; committed_rng=rng;
    performance.r_proposals(end+1,1)=r_proposals;
    performance.h_proposals(end+1,1)=h_proposals;
    performance.r_invalid_proposals(end+1,1)=r_invalid;
    performance.h_invalid_proposals(end+1,1)=h_invalid;
    performance.sweep_seconds(end+1,1)=toc(iteration_clock);
    if completed>cfg.burnin&&mod(completed-cfg.burnin,cfg.thin)==0
        packed=reshape(cache.P,m*m,T);
        packed=packed(sel(:),:)';
        buffer.P_pairs(:,:,end+1)=packed;
        buffer.h(:,:,end+1)=state.h;
        buffer.mu(end+1,:)=state.mu;
        buffer.sig2h(end+1,:)=state.sig2h;
        buffer.sig2r(end+1,:)=state.sig2r;
        parameter_sums.mu=parameter_sums.mu+state.mu;
        parameter_sums.sig2h=parameter_sums.sig2h+state.sig2h;
        parameter_sums.sig2r=parameter_sums.sig2r+state.sig2r;
        parameter_sums.count=parameter_sums.count+1;
        buffer.iterations(end+1)=completed; saved=saved+1;
        if numel(buffer.iterations)>=cfg.chunk_size
            flush_buffer();
            commit_checkpoint();
        end
    end
    % The live state and RNG already delimit a completed sweep in memory.
    % Persist at chunk/call boundaries instead of serializing the expanding
    % retained buffer after every sweep. The coordinator resumes only after
    % this call returns; time limits and caught failures use the final save.
    fprintf('DSC chain %d: sweep %d/%d, %.2fs; R/H proposals %d/%d\n', ...
        cfg.chain_id,completed,cfg.max_iterations,performance.sweep_seconds(end),r_proposals,h_proposals);
end
% Flush even a short final chunk. Checkpoint remembers its new chunk index.
if ~isempty(buffer.iterations), flush_buffer(); end
commit_checkpoint();
parameter_estimate=[];
if parameter_sums.count>0
    parameter_estimate=struct('mu',parameter_sums.mu/parameter_sums.count, ...
        'sig2h',parameter_sums.sig2h/parameter_sums.count, ...
        'sig2r',parameter_sums.sig2r/parameter_sums.count);
    parameter_estimate.delta=parameter_estimate.sig2h;
    parameter_estimate.gamma=parameter_estimate.sig2r;
end
summary=struct('status',status,'reason',stop_reason,'failure_identifier',failure_identifier, ...
    'chain_id',cfg.chain_id,'seed',cfg.seed+cfg.chain_id-1, ...
    'randomize_initial_state',cfg.randomize_initial_state, ...
    'estimate_parameters',cfg.estimate_parameters,'inference_mode',inference_mode, ...
    'parameter_smoothing',cfg.parameter_smoothing, ...
    'mean_model',cfg.mean_model,'schema',schema,'file_prefix',cfg.file_prefix, ...
    'parameter_draws',parameter_sums.count,'parameter_estimate',parameter_estimate, ...
    'completed_iterations',completed,'burnin',cfg.burnin,'saved_draws',saved, ...
    'session_completed_iterations',completed-session_start_iterations, ...
    'unfinished_stage',unfinished_stage,'unfinished_sweep_seconds',unfinished_sweep_seconds, ...
    'T',T,'m',m,'pairs',nr,'first_date',char(string(panel.dates(1),'yyyy-MM-dd')), ...
    'last_date',char(string(panel.dates(end),'yyyy-MM-dd')), ...
    'session_seconds',toc(started),'total_seconds',previous_seconds+toc(started), ...
    'stage_seconds',timing,'r_proposals',sum(performance.r_proposals), ...
    'h_proposals',sum(performance.h_proposals), ...
    'r_invalid_proposals',sum(performance.r_invalid_proposals), ...
    'h_invalid_proposals',sum(performance.h_invalid_proposals), ...
    'implementation',implementation, ...
    'correlation_threads',cfg.correlation_threads, ...
    'mean_sweep_seconds',mean(performance.sweep_seconds), ...
    'median_sweep_seconds',median(performance.sweep_seconds), ...
    'sweep_seconds',performance.sweep_seconds, ...
    'r_proposals_by_sweep',performance.r_proposals, ...
    'h_proposals_by_sweep',performance.h_proposals, ...
    'old_dense_brownian_GiB',8*T*T*(m+nr)/2^30, ...
    'old_4000_draw_h_r_P_GiB',8*T*4000*(m+nr+m*m)/2^30, ...
    'packed_P_2000_draw_GiB',8*T*2000*nr/2^30, ...
    'storage','Postburn-in packed P_pairs and h paths, static mu and sig2h/sig2r draws; no B/V or r history; one checkpoint; chain identity preserved.');
if isfield(priors,'excluded_initial_weeks')
    summary.excluded_initial_weeks=priors.excluded_initial_weeks;
    summary.calibration_start=char(string(priors.calibration_dates(1),'yyyy-MM-dd'));
    summary.calibration_end=char(string(priors.calibration_dates(end),'yyyy-MM-dd'));
end
summary_path=output_path('pilot_summary.json'); atomic_json(summary_path,summary);
result=summary; result.run_dir=run_dir; result.priors=priors;
result.dates=panel.dates; result.tickers=panel.tickers;
result.pair_i=panel.pair_i; result.pair_j=panel.pair_j;
result.paths=struct('checkpoint',checkpoint_path, ...
    'summary',summary_path,'posterior_pattern',output_path('posterior_chunk_*.mat'));
atomic_save(output_path('result.mat'),'result',result);

    function check_time()
        if toc(started)>=cfg.max_seconds, error('dsc:TimeLimit','Sampler wall-clock limit reached.'); end
    end
    function flush_buffer()
        chunk=buffer; chunk.dates=panel.dates; chunk.tickers=panel.tickers;
        chunk.pair_i=panel.pair_i; chunk.pair_j=panel.pair_j; chunk.chain_id=cfg.chain_id;
        chunk.mean_model=cfg.mean_model; chunk.schema=schema;
        atomic_save(output_path(sprintf('posterior_chunk_%06d.mat',next_chunk)),'chunk',chunk);
        next_chunk=next_chunk+1; buffer=empty_buffer;
    end
    function commit_checkpoint()
        saveclock=tic;
        checkpoint=struct('schema',schema,'identity',identity,'state',state,'rng',committed_rng, ...
            'completed',completed,'saved',saved,'next_chunk',next_chunk,'buffer',buffer, ...
            'parameter_sums',parameter_sums, ...
            'performance',performance,'timing',timing,'total_seconds',previous_seconds+toc(started));
        atomic_save(checkpoint_path,'checkpoint',checkpoint);
        timing.checkpoint=timing.checkpoint+toc(saveclock);
    end
end

function [mu0,S0]=mean_prior(priors,m)
% Preserve the historical initial-mean prior while allowing semantic names.
if isfield(priors,'mu0'), mu0=priors.mu0; else, mu0=priors.Bbar; end
if isfield(priors,'mu_cov0'), S0=priors.mu_cov0; else, S0=4*priors.VBbar; end
mu0=mu0(:);
if ~isnumeric(mu0)||~isreal(mu0)||numel(mu0)~=m||any(~isfinite(mu0))|| ...
        ~isnumeric(S0)||~isreal(S0)||~isequal(size(S0),[m m])|| ...
        any(~isfinite(S0(:)))||norm(S0-S0','fro')>1e-12*max(1,norm(S0,'fro'))
    error('dsc:ConstantMeanPrior','A finite m-vector prior mean and symmetric m-by-m covariance are required.');
end
S0=(double(S0)+double(S0)')/2;
[~,bad]=chol(S0);
if bad, error('dsc:ConstantMeanPrior','The constant-mean prior covariance must be positive definite.'); end
end

function [prior_mean,scale]=innovation_prior(priors,kind,width)
% Normalize working vectors only: keep legacy scalar priors intact in saved
% identities, while new calibration supplies one scale for every coordinate.
mean_name=['sig2' kind '_mean']; scale_name=['ig_scale_' kind];
if ~all(isfield(priors,{mean_name,scale_name,'ig_shape'}))|| ...
        ~isnumeric(priors.ig_shape)||~isreal(priors.ig_shape)|| ...
        ~isscalar(priors.ig_shape)||~isfinite(priors.ig_shape)||priors.ig_shape<=1
    error('dsc:InnovationPrior','A finite inverse-gamma shape greater than one is required.');
end
prior_mean=prior_vector(priors.(mean_name),width);
scale=prior_vector(priors.(scale_name),width);
if any(abs(prior_mean-scale/(priors.ig_shape-1))>1e-12*max(1,abs(prior_mean)))
    error('dsc:InnovationPrior','Innovation prior means and inverse-gamma scales disagree.');
end
end

function value=prior_vector(value,width)
if ~isnumeric(value)||~isreal(value)||~isvector(value)|| ...
        ~any(numel(value)==[1 width])||any(~isfinite(value(:)))||any(value(:)<=0)
    error('dsc:InnovationPrior','Innovation prior means/scales must be positive finite scalars or coordinate vectors.');
end
value=double(value(:)');
if isscalar(value), value=repmat(value,1,width); end
end

function state=dispersed_initial_state(state,priors,estimate_parameters,T,m,nr,mu0,S0,scale_h,scale_r)
% Independently seeded parameter draws and bounded state perturbations.
if estimate_parameters
    state.sig2h=1./gamrnd(priors.ig_shape*ones(m,1),1./scale_h(:))';
    state.sig2r=1./gamrnd(priors.ig_shape*ones(nr,1),1./scale_r(:))';
end
parameters=validate_fixed_parameters(struct('sig2h',state.sig2h,'sig2r',state.sig2r),m,nr);
state.sig2h=parameters.sig2h; state.sig2r=parameters.sig2r;
state.mu=(mu0+chol(S0,'lower')*(2*tanh(randn(m,1)/2)))';
hscale=min(sqrt((priors.h0_scale+1)*state.sig2h),1);
rscale=min(sqrt((priors.r0_scale+1)*state.sig2r),.15);
h0=priors.mh0(:)'+hscale.*tanh(randn(1,m));
r0=priors.mr0(:)'+rscale.*tanh(randn(1,nr));
state.h=repmat(h0,T,1); state.r=repmat(r0,T,1);
if any(~isfinite([state.mu(:);state.h(:);state.r(:)]))
    error('dsc:InitialState','Dispersed initial states must be finite.');
end
end

function parameters=validate_fixed_parameters(parameters,m,nr)
% Fixed-parameter smoothing fixes evolution variances, never the mean.
if ~isstruct(parameters)||~isscalar(parameters)|| ...
        ~all(isfield(parameters,{'sig2h','sig2r'}))
    error('dsc:FixedParameters','fixed_parameters must contain sig2h and sig2r.');
end
if isfield(parameters,'V')
    error('dsc:FixedParameters','Historical V is unsupported; use compact constant-mean parameters.');
end
for entry={{'sig2h',m},{'sig2r',nr}}
    name=entry{1}{1}; width=entry{1}{2}; values=parameters.(name);
    if ~isnumeric(values)||~isreal(values)||~isequal(size(values),[1 width])|| ...
            any(~isfinite(values))||any(values<=0)
        error('dsc:FixedParameters','Fixed %s must be a finite positive 1-by-%d vector.',name,width);
    end
end
if isfield(parameters,'mu')
    mu=parameters.mu;
    if ~isnumeric(mu)||~isreal(mu)||numel(mu)~=m||any(~isfinite(mu(:)))
        error('dsc:FixedParameters','Saved mu must be a finite m-vector; it initializes but does not fix the mean.');
    end
    mu=double(mu(:)');
else, mu=[]; end
parameters=struct('sig2h',double(parameters.sig2h),'sig2r',double(parameters.sig2r));
if ~isempty(mu), parameters.mu=mu; end
end

function draws=validate_parameter_draws(draws,m,nr)
if ~isstruct(draws)||~isscalar(draws)||~all(isfield(draws,{'sig2h','sig2r'}))
    error('dsc:FixedParameters','fixed_parameter_draws must contain sig2h and sig2r.');
end
if isfield(draws,'V')
    error('dsc:FixedParameters','Historical V draws are unsupported; use compact constant-mean parameters.');
end
sig2h=draws.sig2h; sig2r=draws.sig2r; n=size(sig2h,1);
if n<1||~isnumeric(sig2h)||~isreal(sig2h)||~isequal(size(sig2h),[n m])|| ...
        ~isnumeric(sig2r)||~isreal(sig2r)||~isequal(size(sig2r),[n nr])|| ...
        any(~isfinite(sig2h(:)))||any(sig2h(:)<=0)|| ...
        any(~isfinite(sig2r(:)))||any(sig2r(:)<=0)
    error('dsc:FixedParameters','Fixed variance draws must be positive finite draw-by-coordinate arrays.');
end
if isfield(draws,'mu')
    mu=draws.mu;
    if ~isnumeric(mu)||~isreal(mu)||~isequal(size(mu),[n m])||any(~isfinite(mu(:)))
        error('dsc:FixedParameters','Saved mu draws must be a finite draw-by-series array.');
    end
    mu=double(mu);
else, mu=[]; end
draws=struct('sig2h',double(sig2h),'sig2r',double(sig2r));
if ~isempty(mu), draws.mu=mu; end
end

function state=apply_parameter_draw(state,draws,iteration)
params=parameter_draw_at(draws,iteration);
state.sig2h=params.sig2h; state.sig2r=params.sig2r;
end

function params=parameter_draw_at(draws,iteration)
n=size(draws.sig2h,1); d=mod(iteration-1,n)+1;
params=struct('sig2h',draws.sig2h(d,:),'sig2r',draws.sig2r(d,:));
if isfield(draws,'mu'), params.mu=draws.mu(d,:); end
end

function mu=draw_constant_mean_cached(Y,h,cache,mu0,S0,z,check)
% Exact static Gaussian update, reusing observed-subspace correlation factors.
% Q=L_P\diag(exp(-h/2)) whitens each observed covariance without refactoring.
check(); m=size(Y,2); prior_factor=chol(S0);
precision=prior_factor\(prior_factor'\eye(m)); rhs=precision*mu0;
for k=1:numel(cache.pattern_indices)
    ix=cache.pattern_indices{k}; dates=cache.pattern_dates{k}; n=numel(ix);
    if n==0, continue; end
    check();
    Q=pagemldivide(cache.L(1:n,1:n,dates),eye(n));
    Q=Q.*reshape(exp(-h(dates,ix)'/2),1,n,numel(dates));
    Qt=pagetranspose(Q);
    precision(ix,ix)=precision(ix,ix)+sum(pagemtimes(Qt,Q),3);
    whitened=pagemtimes(Q,reshape(Y(dates,ix)',n,1,numel(dates)));
    rhs(ix)=rhs(ix)+sum(pagemtimes(Qt,whitened),3);
end
factor=chol((precision+precision')/2);
mu=factor\(factor'\rhs)+factor\z;
if any(~isfinite(mu)), error('dsc:ConstantMeanDraw','The sampled mean must be finite.'); end
check();
end

function cache=build_cache(r,obs,m,backend,mask,threads,check)
% Cholesky factors occupy the observed subspace at each cube's upper left.
T=size(r,1); cache.P=zeros(m,m,T); cache.L=zeros(m,m,T); cache.logdet=zeros(T,1); cache.diagonal=[]; x=[];
cache.all_observed=all(mask(:));
if cache.all_observed
    cache.pattern_indices={1:m}; cache.pattern_dates={1:T};
else
    [patterns,~,groups]=unique(mask,'rows');
    cache.pattern_indices=cell(size(patterns,1),1); cache.pattern_dates=cache.pattern_indices;
    for k=1:size(patterns,1)
        cache.pattern_indices{k}=find(patterns(k,:)); cache.pattern_dates{k}=find(groups==k);
    end
end
if strcmp(backend,'mex')||strcmp(backend,'mex_new')
    check();
    if strcmp(backend,'mex_new')
        % Build the accepted-state cache from a deterministic cold solve.
        % This makes exact checkpoint/RNG replay independent of prior caches.
        [~,cache.P,cache.L,cache.logdet,~,cache.diagonal]=dsc_corr_batch_mex_new(r,zeros(T,m),mask,threads);
    else
        [~,cache.P,cache.L,cache.logdet]=dsc_corr_batch_mex(r,zeros(T,m),mask,threads);
    end
    check(); return
end
for t=1:T
    if mod(t,16)==1, check(); end
    if strcmp(backend,'matlab_new'), [C,~,x]=dsc_corr_matrix_new(r(t,:),x);
    else, [C,~,x]=dsc_corr_matrix(r(t,:),x); end
    cache.P(:,:,t)=C; ix=obs{t}; n=numel(ix);
    if n==0, continue; end
    L=chol(C(ix,ix),'lower'); cache.L(1:n,1:n,t)=L;
    cache.logdet(t)=2*sum(log(diag(L)));
end
end

function ll=cached_likelihood(E,h,cache,mask)
% Pattern batches keep all dates, including missing dates, in the state walk.
[T,m]=size(E); z=E.*exp(-h/2); z(~mask)=0; quadratic=0;
if cache.all_observed
    u=pagemldivide(cache.L,reshape(z',m,1,T)); quadratic=sum(u(:).^2);
else
    for k=1:numel(cache.pattern_indices)
        ix=cache.pattern_indices{k}; dates=cache.pattern_dates{k}; n=numel(ix);
        if n==0, continue; end
        u=pagemldivide(cache.L(1:n,1:n,dates),reshape(z(dates,ix)',n,1,numel(dates)));
        quadratic=quadratic+sum(u(:).^2);
    end
end
ll=-.5*(nnz(mask)*log(2*pi)+sum(h(mask))+sum(cache.logdet)+quadratic);
if ~isfinite(ll), ll=-Inf; end
end

function [directions,z,u]=volatility_cache(E,h,cache,mask)
[T,m]=size(E);
if cache.all_observed
    directions=pagemldivide(cache.L,eye(m));
else
    directions=zeros(m,m,T);
    for k=1:numel(cache.pattern_indices)
        ix=cache.pattern_indices{k}; dates=cache.pattern_dates{k}; n=numel(ix);
        if n==0, continue; end
        directions(1:n,ix,dates)=pagemldivide(cache.L(1:n,1:n,dates),eye(n));
    end
end
z=E.*exp(-h/2); z(~mask)=0;
u=pagemtimes(directions,reshape(z',m,1,T));
end

function ll=correlation_likelihood(candidate,j,r,Z,constant,obs,m,backend,mask,threads,diagonal_guess,check)
r(:,j)=candidate; ll=0; x=[];
try
    if strcmp(backend,'mex')||strcmp(backend,'mex_new')
        check();
        if strcmp(backend,'mex_new')
            if isempty(diagonal_guess)
                ll=dsc_corr_batch_mex_new(r,Z,mask,threads)+constant;
            else
                ll=dsc_corr_batch_mex_new(r,Z,mask,threads,diagonal_guess)+constant;
            end
        else
            ll=dsc_corr_batch_mex(r,Z,mask,threads)+constant;
        end
        check();
    else
        for t=1:size(Z,1)
            if mod(t,16)==1, check(); end
            if strcmp(backend,'matlab_new')
                [C,~,x]=dsc_corr_matrix_new(r(t,:),x);
            else
                [C,~,x]=dsc_corr_matrix(r(t,:),x);
            end
            ix=obs{t};
            if isempty(ix), continue; end
            L=chol(C(ix,ix),'lower'); u=L\Z(t,ix)';
            ll=ll-.5*(2*sum(log(diag(L)))+u'*u);
        end
        ll=ll+constant;
    end
catch problem
    if strcmp(problem.identifier,'dsc:TimeLimit'), rethrow(problem); end
    if startsWith(problem.identifier,'dsc:Correlation'), ll=-Inf; else, rethrow(problem); end
end
if ~isfinite(ll), ll=-Inf; end
end

function [x,ll,ntry,ninvalid]=ellipse(x,mu,sd,c,ll,f,check)
if ~isfinite(ll), error('dsc:InitialLikelihood','Current state has a nonfinite log likelihood.'); end
v=dsc_random_walk_draw(randn(numel(x),1),c,sd); centered=x-mu;
threshold=ll+log(rand); theta=2*pi*rand; low=theta-2*pi; high=theta;
ninvalid=0;
for ntry=1:1000
    check(); proposal=centered*cos(theta)+v*sin(theta)+mu;
    candidate_ll=f(proposal);
    if ~isfinite(candidate_ll), ninvalid=ninvalid+1; end
    if candidate_ll>threshold, x=proposal; ll=candidate_ll; return; end
    if theta<0, low=theta; else, high=theta; end
    theta=low+(high-low)*rand;
end
error('dsc:SliceLimit','Elliptical slice update exceeded 1000 proposals.');
end

function s=draw_variances(path,mu,c,shape,scale)
% Exact conjugate draw: the initial state's variance (c+1)*s is also
% inverse-gamma in s, so all T path terms enter the posterior directly.
T=size(path,1); changes=diff(path,1,1); posterior_shape=T/2+shape;
s=zeros(1,size(path,2));
for j=1:size(path,2)
    d=path(1,j)-mu(j);
    posterior_scale=sum(changes(:,j).^2)/2+d*d/(2*(c+1))+scale(j);
    s(j)=1/gamrnd(posterior_shape,1/posterior_scale);
end
end

function atomic_save(path,name,value)
folder=fileparts(path); temporary=[tempname(folder),'.mat'];
container=struct(); container.(name)=value;
save(temporary,'-struct','container','-v7.3');
movefile(temporary,path,'f');
end

function atomic_json(path,value)
folder=fileparts(path); temporary=[tempname(folder),'.json'];
fid=fopen(temporary,'w');
if fid<0, error('dsc:OutputFile','Cannot write summary.'); end
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true)); clear cleanup
movefile(temporary,path,'f');
end

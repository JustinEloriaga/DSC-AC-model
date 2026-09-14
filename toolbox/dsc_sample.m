function result = dsc_sample(panel,priors,cfg,run_dir)
% Joint p=0 full-history DSC smoother. All files belong to one chain.
% A time limit aborts an unfinished sweep; only completed sweeps commit.
started=tic;
if cfg.p~=0, error('dsc:UnsupportedVAR','Only p=0 is implemented.'); end
Y=panel.returns; [T,m]=size(Y); nr=m*(m-1)/2;
if ~isequal(size(panel.observation_mask),size(Y))||m<2, error('dsc:PanelShape','Invalid panel/mask dimensions.'); end
mask=logical(panel.observation_mask)&isfinite(Y);
if numel(panel.dates)~=T||any(diff(panel.dates)<=seconds(0))
    error('dsc:PanelDates','Dates must increase strictly and match returns.');
end
if any(~isfinite(Y(logical(panel.observation_mask))))
    error('dsc:ObservedNaN','Observed returns must be finite.');
end
if ~isequal(priors.tickers,panel.tickers)||~isequal(priors.source_hash,panel.source_hash)
    error('dsc:PriorIdentity','Priors must match the panel tickers and source hash.');
end
if ~isfield(cfg,'chain_id'), cfg.chain_id=1; end
if ~isfield(cfg,'resume'), cfg.resume=false; end
if ~isfield(cfg,'correlation_backend'), cfg.correlation_backend='matlab'; end
if ~isfield(cfg,'correlation_threads'), cfg.correlation_threads=1; end
if ~isscalar(cfg.correlation_threads)||~isfinite(cfg.correlation_threads)|| ...
        cfg.correlation_threads<1||cfg.correlation_threads>64||cfg.correlation_threads~=floor(cfg.correlation_threads)
    error('dsc:SamplerConfig','Correlation thread count must be an integer from 1 to 64.');
end
cfg.correlation_backend=validatestring(cfg.correlation_backend,{'auto','matlab','mex'});
if strcmp(cfg.correlation_backend,'auto')
    if exist('dsc_corr_batch_mex','file')==3
        cfg.correlation_backend='mex';
    else
        cfg.correlation_backend='matlab';
    end
end
if strcmp(cfg.correlation_backend,'mex')&&exist('dsc_corr_batch_mex','file')~=3
    error('dsc:MissingKernel','Run build_dsc_mex before selecting the MEX backend.');
end
implementation=dsc_sampler_signature(cfg.correlation_backend);
if cfg.burnin<0||cfg.thin<1||cfg.chunk_size<1||cfg.max_iterations<1||cfg.max_seconds<=0
    error('dsc:SamplerConfig','Invalid iteration, thinning, chunk or time limit.');
end
if any([cfg.burnin cfg.thin cfg.chunk_size cfg.max_iterations]~=floor([cfg.burnin cfg.thin cfg.chunk_size cfg.max_iterations]))
    error('dsc:SamplerConfig','Iteration counts must be integers.');
end
if ~exist(run_dir,'dir'), mkdir(run_dir); end
checkpoint_path=fullfile(run_dir,'checkpoint.mat');
fixed_cfg=cfg;
for field={'max_iterations','max_seconds','resume','output_root'}
    if isfield(fixed_cfg,field{1}), fixed_cfg=rmfield(fixed_cfg,field{1}); end
end
identity=struct('dates',panel.dates,'tickers',panel.tickers,'returns',Y, ...
    'mask',mask,'source_hash',panel.source_hash,'cfg',fixed_cfg,'priors',priors, ...
    'implementation',implementation);
obs=cell(T,1); for t=1:T, obs{t}=find(mask(t,:)); end
sel=tril(true(m),-1);
[expected_i,expected_j]=find(sel);
if ~isequal(panel.pair_i(:),expected_i)||~isequal(panel.pair_j(:),expected_j)
    error('dsc:PairOrder','Pair indices must use column-major strict lower triangle order.');
end
timing=struct('mean',0,'mean_variance',0,'correlation',0,'correlation_variance',0, ...
    'volatility',0,'volatility_variance',0,'checkpoint',0);
empty_buffer=struct('P_pairs',zeros(T,nr,0),'h',zeros(T,m,0), ...
    'B',zeros(T,m,0),'iterations',zeros(1,0));
if cfg.resume
    if ~exist(checkpoint_path,'file'), error('dsc:ResumeMissing','Checkpoint does not exist.'); end
    loaded=load(checkpoint_path,'checkpoint'); cp=loaded.checkpoint;
    if ~isequaln(cp.identity,identity), error('dsc:ResumeMismatch','Panel, priors or fixed configuration differs.'); end
    state=cp.state; rng(cp.rng); completed=cp.completed; saved=cp.saved;
    next_chunk=cp.next_chunk; buffer=cp.buffer; diagnostics=cp.diagnostics;
    timing=cp.timing; previous_seconds=cp.total_seconds;
else
    if exist(checkpoint_path,'file')||~isempty(dir(fullfile(run_dir,'posterior_chunk_*.mat')))
        error('dsc:ExistingRun','Run directory already contains sampler output; explicitly resume it.');
    end
    rng(cfg.seed+cfg.chain_id-1,'twister');
    state.B=repmat(priors.Bbar',T,1);
    state.V=priors.V0B/(priors.nuB-m-1);
    state.h=repmat(priors.mh0',T,1);
    state.r=repmat(priors.mr0',T,1);
    state.sig2h=priors.sig2h_mean*ones(1,m);
    state.sig2r=priors.sig2r_mean*ones(1,nr);
    completed=0; saved=0; next_chunk=1; buffer=empty_buffer;
    diagnostics=struct('iteration',zeros(0,1),'sig2h',zeros(0,m), ...
        'sig2r',zeros(0,nr),'V_diag',zeros(0,m),'loglik',zeros(0,1), ...
        'r_proposals',zeros(0,1),'h_proposals',zeros(0,1), ...
        'r_invalid_proposals',zeros(0,1),'h_invalid_proposals',zeros(0,1), ...
        'sweep_seconds',zeros(0,1));
    previous_seconds=0;
end
committed_rng=rng;
status='iteration_limit'; stop_reason='Configured iteration limit reached; convergence has not been established.';
failure_identifier='';
session_start_iterations=completed;
unfinished_stage=''; unfinished_sweep_seconds=0;
% Initial cache contains no random draws. Even a tiny budget receives a valid
% zero-sweep checkpoint and explicit time_limit status from the loop below.
cache=build_cache(state.r,obs,m,cfg.correlation_backend,mask,cfg.correlation_threads,@() []);
if ~cfg.resume, commit_checkpoint(); end
while completed<cfg.max_iterations
    before=state; rng_before=rng; timing_before=timing;
    iteration_clock=tic;
    stage_name='mean';
    try
        check_time(); stage=tic;
        state.B=draw_mean(Y,obs,state.h,cache,state.V,priors,@check_time);
        timing.mean=timing.mean+toc(stage);
        stage_name='mean_variance'; stage=tic; err=diff(state.B,1,1); scale=err'*err+priors.V0B;
        state.V=iwishrnd((scale+scale')/2,T-1+priors.nuB);
        timing.mean_variance=timing.mean_variance+toc(stage);
        residual=Y-state.B;
        stage_name='correlation'; stage=tic; ll=cached_likelihood(residual,state.h,cache,obs);
        % E and h stay fixed throughout all correlation-coordinate updates.
        % Standardize once, rather than once per date and proposal.
        Z=residual.*exp(-state.h/2); Z(~mask)=0;
        correlation_constant=-.5*(nnz(mask)*log(2*pi)+sum(state.h(mask)));
        r_proposals=0; r_invalid=0;
        for j=1:nr
            check_time();
            f=@(x) correlation_likelihood(x,j,state.r,Z,correlation_constant, ...
                obs,m,cfg.correlation_backend,mask,cfg.correlation_threads,@check_time);
            [state.r(:,j),ll,ntry,ninvalid]=ellipse(state.r(:,j),priors.mr0(j), ...
                sqrt(state.sig2r(j)),priors.r0_scale,ll,f,@check_time);
            r_proposals=r_proposals+ntry; r_invalid=r_invalid+ninvalid;
        end
        cache_new=build_cache(state.r,obs,m,cfg.correlation_backend,mask,cfg.correlation_threads,@check_time);
        timing.correlation=timing.correlation+toc(stage);
        stage_name='correlation_variance'; stage=tic;
        state.sig2r=draw_variances(state.r,state.sig2r,priors.mr0, ...
            priors.r0_scale,priors.ig_shape,priors.ig_scale_r);
        timing.correlation_variance=timing.correlation_variance+toc(stage);
        stage_name='volatility'; stage=tic; ll=cached_likelihood(residual,state.h,cache_new,obs); h_proposals=0; h_invalid=0;
        for j=1:m
            check_time();
            f=@(x) volatility_likelihood(x,j,state.h,residual,cache_new,obs);
            [state.h(:,j),ll,ntry,ninvalid]=ellipse(state.h(:,j),priors.mh0(j), ...
                sqrt(state.sig2h(j)),priors.h0_scale,ll,f,@check_time);
            h_proposals=h_proposals+ntry; h_invalid=h_invalid+ninvalid;
        end
        timing.volatility=timing.volatility+toc(stage);
        stage_name='volatility_variance'; stage=tic;
        state.sig2h=draw_variances(state.h,state.sig2h,priors.mh0, ...
            priors.h0_scale,priors.ig_shape,priors.ig_scale_h);
        timing.volatility_variance=timing.volatility_variance+toc(stage);
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
    diagnostics.iteration(end+1,1)=completed;
    diagnostics.sig2h(end+1,:)=state.sig2h;
    diagnostics.sig2r(end+1,:)=state.sig2r;
    diagnostics.V_diag(end+1,:)=diag(state.V)';
    diagnostics.loglik(end+1,1)=ll;
    diagnostics.r_proposals(end+1,1)=r_proposals;
    diagnostics.h_proposals(end+1,1)=h_proposals;
    diagnostics.r_invalid_proposals(end+1,1)=r_invalid;
    diagnostics.h_invalid_proposals(end+1,1)=h_invalid;
    diagnostics.sweep_seconds(end+1,1)=toc(iteration_clock);
    if completed>cfg.burnin&&mod(completed-cfg.burnin,cfg.thin)==0
        packed=zeros(T,nr);
        for t=1:T, C=cache.P(:,:,t); packed(t,:)=C(sel)'; end
        buffer.P_pairs(:,:,end+1)=packed;
        buffer.h(:,:,end+1)=state.h;
        buffer.B(:,:,end+1)=state.B;
        buffer.iterations(end+1)=completed; saved=saved+1;
        if numel(buffer.iterations)>=cfg.chunk_size, flush_buffer(); end
    end
    commit_checkpoint();
    fprintf('DSC chain %d: sweep %d/%d, %.2fs; R/H proposals %d/%d\n', ...
        cfg.chain_id,completed,cfg.max_iterations,diagnostics.sweep_seconds(end),r_proposals,h_proposals);
end
% Flush even a short final chunk. Checkpoint remembers its new chunk index.
if ~isempty(buffer.iterations), flush_buffer(); end
commit_checkpoint();
atomic_save(fullfile(run_dir,'diagnostics.mat'),'diagnostics',diagnostics);
summary=struct('status',status,'reason',stop_reason,'failure_identifier',failure_identifier, ...
    'convergence_established',false,'chain_id',cfg.chain_id,'seed',cfg.seed+cfg.chain_id-1, ...
    'completed_iterations',completed,'burnin',cfg.burnin,'saved_draws',saved, ...
    'session_completed_iterations',completed-session_start_iterations, ...
    'unfinished_stage',unfinished_stage,'unfinished_sweep_seconds',unfinished_sweep_seconds, ...
    'T',T,'m',m,'pairs',nr,'first_date',char(string(panel.dates(1),'yyyy-MM-dd')), ...
    'last_date',char(string(panel.dates(end),'yyyy-MM-dd')), ...
    'session_seconds',toc(started),'total_seconds',previous_seconds+toc(started), ...
    'stage_seconds',timing,'r_proposals',sum(diagnostics.r_proposals), ...
    'h_proposals',sum(diagnostics.h_proposals), ...
    'r_invalid_proposals',sum(diagnostics.r_invalid_proposals), ...
    'h_invalid_proposals',sum(diagnostics.h_invalid_proposals), ...
    'implementation',implementation, ...
    'correlation_threads',cfg.correlation_threads, ...
    'mean_sweep_seconds',mean(diagnostics.sweep_seconds), ...
    'median_sweep_seconds',median(diagnostics.sweep_seconds), ...
    'sweep_seconds',diagnostics.sweep_seconds, ...
    'r_proposals_by_sweep',diagnostics.r_proposals, ...
    'h_proposals_by_sweep',diagnostics.h_proposals, ...
    'old_dense_brownian_GiB',8*T*T*(m+nr)/2^30, ...
    'old_4000_draw_h_r_P_GiB',8*T*4000*(m+nr+m*m)/2^30, ...
    'packed_P_2000_draw_GiB',8*T*2000*nr/2^30, ...
    'storage','Postburn-in packed P_pairs plus h/B; no r history; one checkpoint; chain identity preserved.');
summary_path=fullfile(run_dir,'pilot_summary.json'); atomic_json(summary_path,summary);
result=summary; result.run_dir=run_dir; result.priors=priors;
result.dates=panel.dates; result.tickers=panel.tickers;
result.pair_i=panel.pair_i; result.pair_j=panel.pair_j;
result.paths=struct('checkpoint',checkpoint_path,'diagnostics',fullfile(run_dir,'diagnostics.mat'), ...
    'summary',summary_path,'posterior_pattern',fullfile(run_dir,'posterior_chunk_*.mat'));
atomic_save(fullfile(run_dir,'result.mat'),'result',result);

    function check_time()
        if toc(started)>=cfg.max_seconds, error('dsc:TimeLimit','Sampler wall-clock limit reached.'); end
    end
    function flush_buffer()
        chunk=buffer; chunk.dates=panel.dates; chunk.tickers=panel.tickers;
        chunk.pair_i=panel.pair_i; chunk.pair_j=panel.pair_j; chunk.chain_id=cfg.chain_id;
        atomic_save(fullfile(run_dir,sprintf('posterior_chunk_%06d.mat',next_chunk)),'chunk',chunk);
        next_chunk=next_chunk+1; buffer=empty_buffer;
    end
    function commit_checkpoint()
        saveclock=tic;
        checkpoint=struct('identity',identity,'state',state,'rng',committed_rng, ...
            'completed',completed,'saved',saved,'next_chunk',next_chunk,'buffer',buffer, ...
            'diagnostics',diagnostics,'timing',timing,'total_seconds',previous_seconds+toc(started));
        atomic_save(checkpoint_path,'checkpoint',checkpoint);
        timing.checkpoint=timing.checkpoint+toc(saveclock);
    end
end

function B=draw_mean(Y,obs,h,cache,V,priors,check)
[T,m]=size(Y); means=zeros(T,m); variances=zeros(m,m,T);
mu=priors.Bbar; S=4*priors.VBbar;
for t=1:T
    if mod(t,16)==1, check(); end
    % B_1 ~ N(Bbar,4*VBbar), independently of V. Only transitions
    % B_t-B_(t-1), t>=2, enter the inverse-Wishart full conditional.
    if t==1, predicted=S; else, predicted=S+V; end
    ix=obs{t};
    if ~isempty(ix)
        sd=exp(h(t,ix)/2); R=(sd'*sd).*cache.P(ix,ix,t);
        F=predicted(ix,ix)+R; gain=predicted(:,ix)/F;
        mu=mu+gain*(Y(t,ix)'-mu(ix));
        S=predicted-gain*predicted(ix,:); S=(S+S')/2;
    else
        S=predicted;
    end
    means(t,:)=mu'; variances(:,:,t)=S;
end
B=zeros(T,m); B(T,:)=normal_draw(mu,S)';
for t=T-1:-1:1
    if mod(t,16)==1, check(); end
    St=variances(:,:,t); gain=St/(St+V);
    mu=means(t,:)'+gain*(B(t+1,:)'-means(t,:)');
    C=St-gain*St; B(t,:)=normal_draw(mu,(C+C')/2)';
end
end

function x=normal_draw(mu,S)
S=(S+S')/2; [L,bad]=chol(S,'lower');
if bad
    [U,D]=eig(S,'vector'); tol=1e-10*max(1,max(abs(D)));
    if min(D)<-tol, error('dsc:StateCovariance','State covariance has a material negative eigenvalue.'); end
    L=U*diag(sqrt(max(D,0)));
end
x=mu+L*randn(numel(mu),1);
end

function cache=build_cache(r,obs,m,backend,mask,threads,check)
T=size(r,1); cache.P=zeros(m,m,T); cache.L=cell(T,1); cache.logdet=zeros(T,1); x=[];
if strcmp(backend,'mex')
    check();
    [~,cache.P,L,cache.logdet]=dsc_corr_batch_mex(r,zeros(T,m),mask,threads);
    check();
    for t=1:T
        n=numel(obs{t}); cache.L{t}=L(1:n,1:n,t);
    end
    return
end
for t=1:T
    if mod(t,16)==1, check(); end
    [C,~,x]=dsc_corr_matrix(r(t,:),x); cache.P(:,:,t)=C;
    ix=obs{t}; cache.L{t}=chol(C(ix,ix),'lower');
    cache.logdet(t)=2*sum(log(diag(cache.L{t})));
end
end

function ll=cached_likelihood(E,h,cache,obs)
ll=0;
for t=1:size(E,1)
    ix=obs{t};
    if isempty(ix), continue; end
    z=E(t,ix)'.*exp(-h(t,ix)'/2); u=cache.L{t}\z;
    ll=ll-.5*(numel(ix)*log(2*pi)+sum(h(t,ix))+cache.logdet(t)+u'*u);
end
if ~isfinite(ll), ll=-Inf; end
end

function ll=volatility_likelihood(x,j,h,E,cache,obs)
h(:,j)=x; ll=cached_likelihood(E,h,cache,obs);
end

function ll=correlation_likelihood(candidate,j,r,Z,constant,obs,m,backend,mask,threads,check)
r(:,j)=candidate; ll=0; x=[];
try
    if strcmp(backend,'mex')
        check(); ll=dsc_corr_batch_mex(r,Z,mask,threads)+constant; check();
    else
        for t=1:size(Z,1)
            if mod(t,16)==1, check(); end
            [C,~,x]=dsc_corr_matrix(r(t,:),x); ix=obs{t};
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

function s=draw_variances(path,s,mu,c,shape,scale)
T=size(path,1); changes=diff(path,1,1); posterior_shape=(T-1)/2+shape;
for j=1:size(path,2)
    posterior_scale=sum(changes(:,j).^2)/2+scale;
    proposal=1/gamrnd(posterior_shape,1/posterior_scale);
    d=path(1,j)-mu(j);
    logratio=-.5*log(proposal/s(j))-.5*d*d/(c+1)*(1/proposal-1/s(j));
    if log(rand)<logratio, s(j)=proposal; end
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

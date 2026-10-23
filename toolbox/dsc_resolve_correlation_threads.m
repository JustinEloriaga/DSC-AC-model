function [threads,info]=dsc_resolve_correlation_threads(requested,num_chains,parallel_chains,panel,priors,backend)
% Resolve native threads once, before dispatching chains to workers.
% Sequential chains share the full CPU budget; parallel chains divide it.
% The original three-argument API remains an allocation policy. Optional
% panel/priors/backend context tunes auto mex_new counts inside that budget;
% explicit numeric counts never benchmark. The calibration workload is deterministic,
% consumes no RNG draws and writes no files. Resolve once and pass the
% resulting numeric count to every sampler batch and checkpoint.
if nargin<1, requested='auto'; end
if nargin<2, num_chains=1; end
if nargin<3, parallel_chains=false; end
validateattributes(num_chains,{'numeric'},{'real','scalar','finite','integer','positive'});
validateattributes(parallel_chains,{'logical','numeric'},{'real','scalar','binary'});
info=struct('status','explicit','budget',[],'candidates',[], ...
    'selected',[],'cold_seconds',[],'warm_seconds',[],'score_seconds',[], ...
    'calibration_seconds',0,'calibration_dates',0,'reason','');
if (ischar(requested)&&isrow(requested))|| ...
        (isstring(requested)&&isscalar(requested)&&~ismissing(requested))
    if strcmp(requested,'auto')
        try
            cores=feature('numcores');
        catch
            cores=1;
        end
        if ~isnumeric(cores)||~isreal(cores)||~isscalar(cores)|| ...
                ~isfinite(cores)||cores<1
            cores=1;
        end
        simultaneous=1;
        if parallel_chains
            simultaneous=num_chains;
            % Client cores are a shared budget only for a local process pool.
            % Never start a pool or infer remote worker resources from them.
            if exist('gcp','file')~=0
                try
                    pool=gcp('nocreate');
                    if ~isempty(pool)&&~isa(pool.Cluster,'parallel.cluster.Local')
                        cores=1;
                    end
                catch
                    cores=1;
                end
            end
        end
        threads=max(1,min(64,floor(cores/simultaneous)));
        info.status='allocation'; info.budget=threads; info.selected=threads;
        if nargin>=6&&threads>1&& ...
                ((ischar(backend)&&isrow(backend))|| ...
                (isstring(backend)&&isscalar(backend)&&~ismissing(backend)))&& ...
                strcmp(backend,'mex_new')&&exist('dsc_corr_batch_mex_new','file')==3
            [threads,info]=tune_native_threads(threads,panel,priors,info);
        end
        return
    end
    requested=str2double(requested);
end
validateattributes(requested,{'numeric'},{'real','scalar','finite','integer','>=',1,'<=',64});
threads=double(requested);
info.selected=threads;
end

function [selected,info]=tune_native_threads(budget,panel,priors,info)
selected=budget;
started=tic;
% An incomplete optional context leaves the three-argument policy intact.
if ~isstruct(panel)||~isscalar(panel)||~isfield(panel,'returns')|| ...
        ~isfield(panel,'observation_mask')||~isstruct(priors)||~isscalar(priors)|| ...
        ~isfield(priors,'mh0')||~isfield(priors,'mr0')|| ...
        ~isfield(priors,'sig2r_mean')
    info.reason='Optional calibration context is incomplete.'; return
end
Y=panel.returns; mask=panel.observation_mask;
if ~isnumeric(Y)||~isreal(Y)||issparse(Y)||ndims(Y)~=2|| ...
        ~isequal(size(mask),size(Y))||isempty(Y)||size(Y,2)<2|| ...
        ~(islogical(mask)||isnumeric(mask))||~isreal(mask)||issparse(mask)
    info.reason='Optional panel context is invalid.'; return
end
[T,m]=size(Y); q=m*(m-1)/2;
if isfield(priors,'mu0'), mu=priors.mu0;
elseif isfield(priors,'Bbar'), mu=priors.Bbar;
else, info.reason='Optional constant-mean prior is absent.'; return
end
h=priors.mh0; r0=priors.mr0; variance=priors.sig2r_mean;
if ~isnumeric(mu)||~isreal(mu)||numel(mu)~=m||any(~isfinite(mu(:)))|| ...
        ~isnumeric(h)||~isreal(h)||numel(h)~=m||any(~isfinite(h(:)))|| ...
        ~isnumeric(r0)||~isreal(r0)||numel(r0)~=q||any(~isfinite(r0(:)))|| ...
        ~isnumeric(variance)||~isreal(variance)||~isscalar(variance)|| ...
        ~isfinite(variance)||variance<=0
    info.reason='Optional prior context is invalid.'; return
end
if T<32
    selected=1; info.status='small_panel'; info.selected=1; info.candidates=1;
    info.calibration_dates=T;
    info.reason='Fewer than 32 dates: avoid calibration and native worker overhead.';
    return
end
% Preserve the date workload for supported panels; cap only long histories.
count=min(T,2048);
ix=unique(round(linspace(1,T,count))); count=numel(ix);
mask=logical(mask(ix,:))&isfinite(Y(ix,:));
Z=(double(Y(ix,:))-double(mu(:)')).*exp(-double(h(:)')/2); Z(~mask)=0;
if any(~isfinite(Z(mask)))
    info.reason='Standardized calibration observations are nonfinite.'; return
end
% Smooth prior-scale paths plus small date-to-date changes exercise the
% inverse solver beyond the unusually cheap repeated constant initial path.
% These waves are deterministic arithmetic, not posterior draws.
t=(1:count)'; u=t/max(1,count); sd=sqrt(double(variance));
r=repmat(double(r0(:)'),count,1);
for j=1:q
    phase=j*sqrt(2);
    r(:,j)=r(:,j)+sd*sqrt(count)*(.6*sin(2*pi*u+phase)+ ...
        .3*sin(11*pi*u+phase)+.1*sin(47*pi*u+phase))+ ...
        .5*sd*sin(t*(sqrt(5)+j/7)+phase);
end
proposal=r; proposal(:,1)=proposal(:,1)+.1*sd*sin(t/11);
candidates=unique([1 2.^(1:floor(log2(budget))) budget]);
cold=zeros(3,numel(candidates)); warm=cold;
try
    [~,~,~,~,~,guess]=dsc_corr_batch_mex_new(r,Z,mask,1);
    % Warm every code path before alternating the timed candidate order.
    for k=candidates
        cold_call(r,Z,mask,k); value=dsc_corr_batch_mex_new(proposal,Z,mask,k,guess); %#ok<NASGU>
    end
    repetitions=max(1,min(3,floor(12/m)));
    for trial=1:3
        order=1:numel(candidates); if mod(trial,2)==0, order=fliplr(order); end
        for index=order
            k=candidates(index); clock=tic;
            for rep=1:repetitions, cold_call(r,Z,mask,k); end
            cold(trial,index)=toc(clock)/repetitions;
            clock=tic;
            for rep=1:repetitions, value=dsc_corr_batch_mex_new(proposal,Z,mask,k,guess); end %#ok<NASGU>
            warm(trial,index)=toc(clock)/repetitions;
        end
    end
catch problem
    % A synthetic calibration path cannot invalidate the actual model.
    % Preserve the allocation policy if this auxiliary workload fails.
    if any(strcmp(problem.identifier,{'dsc:KernelMemory','dsc:KernelInternal'}))
        selected=1; info.selected=1; info.status='resource_fallback';
        info.reason=['Auxiliary calibration resource failure: ' problem.identifier];
        info.calibration_seconds=toc(started); return
    end
    if ~startsWith(problem.identifier,'dsc:Correlation'), rethrow(problem); end
    info.reason=['Auxiliary calibration failed: ' problem.identifier];
    info.calibration_seconds=toc(started); return
end
cold=median(cold,1); warm=median(warm,1);
score=.2*cold+.8*warm;
% A five-percent timing tie favors fewer native workers and less contention.
best=min(score); index=find(score<=1.05*best,1,'first');
selected=candidates(index);
info.status='measured'; info.candidates=candidates; info.selected=selected;
info.cold_seconds=cold; info.warm_seconds=warm; info.score_seconds=score;
info.calibration_seconds=toc(started); info.calibration_dates=count;
info.reason='Median of three trials; 20% cache build / 80% nearby proposal; five-percent ties prefer fewer threads.';
end

function value=cold_call(r,Z,mask,threads)
% The production correlation-cache build requests all six native outputs.
[value,~,~,~,~,~]=dsc_corr_batch_mex_new(r,Z,mask,threads);
end

function path = dsc_save_parameters(result,panel,priors,cfg,parameter_dir,calibration_panel)
%DSC_SAVE_PARAMETERS Save retained parameter estimates for later smoothing.
%   Retained parameter draws are saved for reuse. Their posterior averages
%   are also stored as a plug-in summary. A checkpoint's last state
%   (including a warm-up state) is never a parameter estimate.
%   Writes a run-linked, schema-versioned MAT file atomically, preserving
%   the exact likelihood sample, separate excluded calibration weeks and
%   calibrated priors for future checks. Both dated blocks are mandatory.
if ~isstruct(result)||~isscalar(result)||~isfield(result,'status')|| ...
        ~isfield(result,'parameter_estimate')||~isfield(result,'parameter_draws')
    error('dsc:ParametersResult','Sampler result has no retained parameter estimate.');
end
if strcmp(result.status,'failed')
    error('dsc:ParametersFailedRun','Cannot save parameter estimates from a failed sampler run.');
end
if (isstruct(cfg)&&isscalar(cfg)&&isfield(cfg,'estimate_parameters')&& ...
        isequal(cfg.estimate_parameters,false))|| ...
        (isfield(result,'estimate_parameters')&&isequal(result.estimate_parameters,false))|| ...
        (isfield(result,'inference_mode')&&strcmp(result.inference_mode,'fixed_parameter_smoothing'))
    error('dsc:ParametersFixedMode','Fixed-parameter smoothing cannot publish a new parameter estimate.');
end
count=result.parameter_draws;
if ~(isnumeric(count)&&isreal(count)&&isscalar(count)&&isfinite(count)&&count>=2&&count==floor(count))
    error('dsc:ParametersDraws','At least two retained parameter draws are required; warm-up states cannot be saved as estimates.');
end
if isfield(result,'saved_draws')&&~isequal(result.saved_draws,count)
    error('dsc:ParametersDraws','Retained parameter count differs from the saved posterior draw count.');
end
if ~isstruct(cfg)||~isscalar(cfg)||~isfield(cfg,'p')||~isequal(cfg.p,0)
    error('dsc:ParametersConfig','Parameter reuse currently supports the joint p=0 model.');
end
fields={'dates','returns','observation_mask','tickers','pair_i','pair_j','source_hash'};
if nargin<6||~isstruct(calibration_panel)||~isscalar(calibration_panel)|| ...
        ~all(isfield(calibration_panel,fields))||~isstruct(panel)||~isscalar(panel)|| ...
        ~all(isfield(panel,fields))||isempty(panel.dates)||isempty(calibration_panel.dates)
    error('dsc:ParametersCalibration','Provide a separate calibration_panel and a nonempty likelihood training panel.');
end
if ~isfield(cfg,'prior_weeks')||~isnumeric(cfg.prior_weeks)||~isreal(cfg.prior_weeks)||~isscalar(cfg.prior_weeks)|| ...
        ~isfinite(cfg.prior_weeks)||cfg.prior_weeks<1||cfg.prior_weeks~=floor(cfg.prior_weeks)|| ...
        numel(calibration_panel.dates)~=cfg.prior_weeks|| ...
        (isfield(priors,'initial_weeks')&&~isequal(priors.initial_weeks,cfg.prior_weeks))
    error('dsc:ParametersCalibration','Calibration week count must equal cfg.prior_weeks and the saved initial-prior window.');
end
for block={calibration_panel,panel}
    values=block{1}.returns; mask=block{1}.observation_mask;
    rows=numel(block{1}.dates); columns=numel(block{1}.tickers);
    if ~isnumeric(values)||~isreal(values)||~isequal(size(values),[rows columns])|| ...
            ~isequal(size(mask),size(values))|| ...
            ~(islogical(mask)||(isnumeric(mask)&&isreal(mask)&&all(ismember(mask(:),[0 1]))))|| ...
            any(isinf(values(:)))||any(~isfinite(values(logical(mask))))
        error('dsc:ParametersCalibration','Calibration and training blocks must have valid return dimensions, masks and finite observed values.');
    end
end
for field={'tickers','pair_i','pair_j','source_hash'}
    if ~isequal(calibration_panel.(field{1}),panel.(field{1}))
        error('dsc:ParametersCalibration','Calibration and training blocks must have identical series, pair order and source hash.');
    end
end
if ~isfield(result,'run_dir')||~is_text_scalar(result.run_dir)
    error('dsc:ParametersResult','Sampler result must identify its source run directory.');
end
if ~is_text_scalar(parameter_dir)
    error('dsc:ParametersDirectory','parameter_dir must be a nonempty directory path.');
end
parameter_dir=char(parameter_dir);
if ~isfolder(parameter_dir)
    error('dsc:ParametersDirectory','Parameter directory must already exist: %s',parameter_dir);
end
parameter_dir=char(java.io.File(parameter_dir).getCanonicalPath());
% Round once to milliseconds so ISO and epoch timestamps describe one instant.
epoch=round(posixtime(datetime('now','TimeZone','UTC'))*1000)/1000;
estimated=datetime(epoch,'ConvertFrom','posixtime','TimeZone','UTC');
[draws,chain_ids,chain_counts,chain_dirs,chain_file_prefixes,seeds]=collect_parameter_draws(result,count,panel,priors,cfg);
estimate=parameter_draw_mean(draws);
if ~isempty(result.parameter_estimate)&&~estimates_close(canonical_estimate(result.parameter_estimate),estimate)
    error('dsc:ParametersValues','Saved posterior chunks disagree with the sampler parameter estimate.');
end
parameters=struct();
parameters.schema_version=5;
parameters.model='joint_dsc_constant_mean';
parameters.mean_model='constant';
parameters.innovation_variance_aliases=struct('delta','sig2h','gamma','sig2r');
parameters.estimator='posterior_draws';
parameters.estimated_at=char(string(estimated,"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
parameters.estimated_at_posix=epoch;
parameters.estimation_start=panel.dates(1);
parameters.estimation_end=panel.dates(end);
parameters.tickers=panel.tickers;
parameters.pair_i=panel.pair_i;
parameters.pair_j=panel.pair_j;
parameters.source_hash=panel.source_hash;
parameters.calibration_dates=calibration_panel.dates;
parameters.calibration_returns=calibration_panel.returns;
parameters.calibration_mask=logical(calibration_panel.observation_mask);
parameters.calibration_weeks=cfg.prior_weeks;
parameters.calibration_start=calibration_panel.dates(1);
parameters.calibration_end=calibration_panel.dates(end);
parameters.training_dates=panel.dates;
parameters.training_returns=panel.returns;
parameters.training_mask=logical(panel.observation_mask);
parameters.priors=priors;
parameters.parameter_draws=draws;
parameters.parameter_estimate=estimate;
parameters.retained_draws=count;
parameters.chain_ids=chain_ids;
parameters.chain_dirs=chain_dirs;
parameters.chain_file_prefixes=chain_file_prefixes;
parameters.chain_seeds=seeds;
parameters.retained_draws_per_chain=chain_counts;
parameters.estimation_status=result.status;
if isfield(result,'completed_iterations'), parameters.completed_iterations=result.completed_iterations; end
if isfield(result,'implementation'), parameters.implementation=result.implementation; end
if all(isfield(cfg,{'max_iterations','burnin','thin'}))
    parameters.requested_retained_draws=numel(chain_ids)*floor((cfg.max_iterations-cfg.burnin)/cfg.thin);
end
parameters.prior_posterior_summary=dsc_parameter_summary(draws,priors, ...
    panel.tickers,panel.pair_i,panel.pair_j);
parameters.parameter_summary_metadata=struct( ...
    'sample_count',count, ...
    'posterior_std_convention','Sample standard deviation with n-1 denominator.', ...
    'posterior_source','Actual retained draws pooled across chains after per-chain warm-up and thinning.', ...
    'prior_source','Analytic marginal moments of the saved calibrated empirical-Bayes priors; calibration uses data.', ...
    'inverse_gamma_parameterization','p(x) proportional to x^(-shape-1) exp(-scale/x), x > 0.');
source_run_dir=char(java.io.File(char(result.run_dir)).getCanonicalPath());
[source_parent,source_name]=fileparts(source_run_dir);
if strcmp(source_name,'estimation')
    source_run_dir=source_parent;
    [~,source_name]=fileparts(source_run_dir);
end
parameters.run_dir=char(result.run_dir);
parameters.file_prefix='';
if isfield(result,'file_prefix'), parameters.file_prefix=char(result.file_prefix); end
if ~isempty(parameters.file_prefix), source_name=regexprep(parameters.file_prefix,'_+$',''); end
if isempty(source_name)||isempty(regexp(source_name,'^[A-Za-z0-9_-]+$','once'))
    error('dsc:ParametersIdentity','A nonempty safe run identifier is required for the parameter filename.');
end
% The filename prefix is the sole run identity. Keep the legacy field out of
% new bundles so every run-linked artifact is identified consistently.
parameters.estimation_config=struct();
for field={'p','seed','chain_id','burnin','thin','max_iterations','prior_weeks', ...
        'shrinkage','calibration_windows','calibration_window','ig_shape','mean_model','mean_backend', ...
        'num_chains','parallel_chains', ...
        'correlation_backend','correlation_threads','file_prefix'}
    if isfield(cfg,field{1})
        parameters.estimation_config.(field{1})=cfg.(field{1});
    end
end
unique=char(java.util.UUID.randomUUID());
path=fullfile(parameter_dir,[parameters.file_prefix 'parameters.mat']);
temporary=fullfile(parameter_dir,['parameters_pending_' unique '.mat']);
cleanup=onCleanup(@() remove_temporary(temporary));
save(temporary,'parameters','-v7');
% The loader always receives the full data history, including calibration.
% Use exactly its validation before atomically publishing the new artifact.
validation_panel=struct('dates',[calibration_panel.dates(:);panel.dates(:)], ...
    'returns',[calibration_panel.returns;panel.returns], ...
    'observation_mask',[logical(calibration_panel.observation_mask);logical(panel.observation_mask)], ...
    'tickers',panel.tickers,'pair_i',panel.pair_i,'pair_j',panel.pair_j,'source_hash',panel.source_hash);
dsc_load_parameters(parameter_dir,validation_panel,temporary);
[ok,message]=movefile(temporary,path);
if ~ok, error('dsc:ParametersWrite','Cannot publish parameter file: %s',message); end
clear cleanup
fprintf('Saved model parameters estimated at %s UTC through %s: %s\n', ...
    parameters.estimated_at(1:end-1),string(parameters.estimation_end,'yyyy-MM-dd'),path);
end

function tf=is_text_scalar(value)
tf=(ischar(value)&&isrow(value)&&~isempty(value))|| ...
    (isstring(value)&&isscalar(value)&&~ismissing(value)&&strlength(value)>0);
end

function remove_temporary(path)
if isfile(path), delete(path); end
end

function [draws,ids,counts,dirs,prefixes,seeds]=collect_parameter_draws(result,count,panel,priors,cfg)
dirs={result.run_dir};
if isfield(result,'chain_dirs'), dirs=cellstr(string(result.chain_dirs)); end
dirs=cellfun(@(p)char(java.io.File(p).getCanonicalPath()),dirs,'UniformOutput',false);
dirs=dirs(:)';
prefixes={''};
if isfield(result,'chain_file_prefixes'), prefixes=cellstr(string(result.chain_file_prefixes));
elseif isfield(result,'file_prefix'), prefixes={char(result.file_prefix)}; end
prefixes=prefixes(:)';
if isempty(dirs)||numel(prefixes)~=numel(dirs)|| ...
        numel(unique(strcat(string(dirs),filesep,string(prefixes))))~=numel(dirs)
    error('dsc:ParametersIdentity','Chain directory/prefix identities must be nonempty and distinct.');
end
for c=1:numel(prefixes)
    if ~isempty(prefixes{c})&&isempty(regexp(prefixes{c},'^[A-Za-z0-9_-]+$','once'))
        error('dsc:ParametersIdentity','Invalid source chain filename prefix.');
    end
end
m=numel(panel.tickers); nr=numel(panel.pair_i);
mu=zeros(count,m); sig2h=zeros(count,m); sig2r=zeros(count,nr);
draw_ids=zeros(count,1); iterations=zeros(count,1); offset=0;
ids=zeros(1,numel(dirs)); counts=ids; seeds=ids; target=[];
for c=1:numel(dirs)
    cp_path=fullfile(dirs{c},[prefixes{c} 'checkpoint.mat']);
    if ~isfile(cp_path), error('dsc:ParametersIdentity','Missing source chain checkpoint: %s.',cp_path); end
    loaded=load(cp_path,'checkpoint'); cp=loaded.checkpoint; identity=cp.identity;
    if ~isfield(cp,'schema')||~strcmp(cp.schema,'constant_mu_v1')|| ...
            ~isfield(identity.cfg,'mean_model')||~strcmp(identity.cfg.mean_model,'constant')|| ...
            isfield(cp.state,'B')||isfield(cp.state,'V')
        error('dsc:ParametersMeanModel','Only compact constant-mean checkpoints can publish new parameter bundles.');
    end
    if ~isequal(identity.dates,panel.dates)||~isequaln(identity.returns,panel.returns)|| ...
            ~isequal(identity.mask,logical(panel.observation_mask)&isfinite(panel.returns))|| ...
            ~isequal(identity.tickers,panel.tickers)||~isequaln(identity.priors,priors)|| ...
            ~strcmp(identity.source_hash,panel.source_hash)|| ...
            ~identity.cfg.estimate_parameters||identity.cfg.burnin~=cfg.burnin||identity.cfg.thin~=cfg.thin|| ...
            ~isfield(identity.cfg,'file_prefix')||~strcmp(identity.cfg.file_prefix,prefixes{c})
        error('dsc:ParametersIdentity','Source chains must share the saved estimation sample, priors and retention settings.');
    end
    ids(c)=identity.cfg.chain_id; seeds(c)=identity.cfg.seed+ids(c)-1;
    current=identity.cfg;
    for field={'seed','chain_id','max_iterations','max_seconds','resume','output_root','file_prefix'}
        if isfield(current,field{1}), current=rmfield(current,field{1}); end
    end
    if c==1, target=current;
    elseif ~isequaln(current,target), error('dsc:ParametersIdentity','Chains have different fixed configurations.'); end
    files=dir(fullfile(dirs{c},[prefixes{c} 'posterior_chunk_*.mat']));
    if isempty(files), error('dsc:ParametersDraws','No posterior chunks found in %s.',dirs{c}); end
    [~,order]=sort({files.name}); files=files(order); previous=cfg.burnin; start=offset;
    for k=1:numel(files)
        loaded=load(fullfile(files(k).folder,files(k).name),'chunk'); chunk=loaded.chunk;
        n=numel(chunk.iterations); ix=chunk.iterations(:);
        if n<1||any(~isfinite(ix))||any(ix~=floor(ix))|| ...
                any(diff([previous;ix])~=cfg.thin)||ix(end)>cp.completed|| ...
                ~isequal(chunk.chain_id,ids(c))||~isequal(chunk.dates,panel.dates)|| ...
                ~isequal(chunk.tickers,panel.tickers)||~isequal(chunk.pair_i,panel.pair_i)|| ...
                ~isequal(chunk.pair_j,panel.pair_j)
            error('dsc:ParametersIdentity','Invalid chain identity or retained iteration sequence in %s.',files(k).name);
        end
        if ~isfield(chunk,'schema')||~strcmp(chunk.schema,'constant_mu_v1')|| ...
                isfield(chunk,'V')||isfield(chunk,'B')||~isfield(chunk,'mu')|| ...
                ~isequal(size(chunk.mu),[n m])|| ...
                ~isequal(size(chunk.sig2h),[n m])||~isequal(size(chunk.sig2r),[n nr])||offset+n>count
            error('dsc:ParametersDraws','Posterior chunk %s has invalid draw dimensions/count.',files(k).name);
        end
        rows=offset+(1:n); mu(rows,:)=chunk.mu; sig2h(rows,:)=chunk.sig2h; sig2r(rows,:)=chunk.sig2r;
        draw_ids(rows)=ids(c); iterations(rows)=ix; offset=offset+n; previous=ix(end);
    end
    counts(c)=offset-start;
    if counts(c)~=cp.saved, error('dsc:ParametersDraws','Chain draw count disagrees with its checkpoint.'); end
end
if numel(unique(ids))~=numel(ids)||numel(unique(seeds))~=numel(seeds)
    error('dsc:ParametersIdentity','Original chain IDs and effective seeds must be distinct.');
end
if offset~=count
    error('dsc:ParametersDraws','Posterior chunks contain %d parameter draws, expected %d.',offset,count);
end
draws=struct('mu',mu,'sig2h',sig2h,'sig2r',sig2r,'chain_id',draw_ids,'iteration',iterations);
validate_parameter_draws_local(draws,m,nr,count);
end

function estimate=parameter_draw_mean(draws)
estimate=struct('mu',mean(draws.mu,1),'sig2h',mean(draws.sig2h,1), ...
    'sig2r',mean(draws.sig2r,1));
estimate.delta=estimate.sig2h; estimate.gamma=estimate.sig2r;
end

function estimate=canonical_estimate(estimate)
if ~isstruct(estimate)||~isscalar(estimate)||~all(isfield(estimate,{'mu','sig2h','sig2r'}))|| ...
        isfield(estimate,'B')||isfield(estimate,'V')
    error('dsc:ParametersValues','A compact constant-mean parameter estimate is required.');
end
estimate=struct('mu',double(estimate.mu(:))', ...
    'sig2h',double(estimate.sig2h(:))','sig2r',double(estimate.sig2r(:))');
end

function tf=estimates_close(a,b)
tf=isequal(size(a.mu),size(b.mu))&&isequal(size(a.sig2h),size(b.sig2h))&& ...
    isequal(size(a.sig2r),size(b.sig2r))&& ...
    max(abs(a.mu(:)-b.mu(:)))<=1e-12*max(1,max(abs(b.mu(:))))&& ...
    max(abs(a.sig2h(:)-b.sig2h(:)))<=1e-12*max(1,max(abs(b.sig2h(:))))&& ...
    max(abs(a.sig2r(:)-b.sig2r(:)))<=1e-12*max(1,max(abs(b.sig2r(:))));
end

function validate_parameter_draws_local(draws,m,nr,count)
if ~isstruct(draws)||~isscalar(draws)||~all(isfield(draws,{'mu','sig2h','sig2r'}))|| ...
        isfield(draws,'B')||isfield(draws,'V')|| ...
        ~isequal(size(draws.mu),[count m])||~isequal(size(draws.sig2h),[count m])|| ...
        ~isequal(size(draws.sig2r),[count nr])||any(~isfinite(draws.mu(:)))|| ...
        any(~isfinite(draws.sig2h(:)))||any(draws.sig2h(:)<=0)|| ...
        any(~isfinite(draws.sig2r(:)))||any(draws.sig2r(:)<=0)
    error('dsc:ParametersValues','Mean draws must be finite and innovation-variance draws positive with matching dimensions.');
end
end

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
    [ok,message]=mkdir(parameter_dir);
    if ~ok, error('dsc:ParametersDirectory','Cannot create parameter directory: %s',message); end
end
parameter_dir=char(java.io.File(parameter_dir).getCanonicalPath());
% Round once to milliseconds so ISO and epoch timestamps describe one instant.
epoch=round(posixtime(datetime('now','TimeZone','UTC'))*1000)/1000;
estimated=datetime(epoch,'ConvertFrom','posixtime','TimeZone','UTC');
[draws,chain_ids,chain_counts,chain_dirs,seeds]=collect_parameter_draws(result,count,panel,priors,cfg);
estimate=parameter_draw_mean(draws);
if ~isempty(result.parameter_estimate)&&~estimates_close(canonical_estimate(result.parameter_estimate),estimate)
    error('dsc:ParametersValues','Saved posterior chunks disagree with the sampler parameter estimate.');
end
parameters=struct();
parameters.schema_version=4;
parameters.model='joint_dsc_p0';
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
parameters.chain_seeds=seeds;
parameters.retained_draws_per_chain=chain_counts;
parameters.estimation_status=result.status;
if isfield(result,'completed_iterations'), parameters.completed_iterations=result.completed_iterations; end
if isfield(result,'implementation'), parameters.implementation=result.implementation; end
if all(isfield(cfg,{'max_iterations','burnin','thin'}))
    parameters.requested_retained_draws=numel(chain_ids)*floor((cfg.max_iterations-cfg.burnin)/cfg.thin);
end
parameters.convergence=struct('status','not_checked','passed',false);
if isfield(result,'convergence'), parameters.convergence=result.convergence; end
m=numel(panel.tickers); nr=numel(panel.pair_i);
expected_quantities=numel(panel.dates)*(m+nr);
dsc_validate_convergence(parameters.convergence,chain_ids,chain_counts,expected_quantities);
parameters.convergence_established=parameters.convergence.passed;
source_run_dir=char(java.io.File(char(result.run_dir)).getCanonicalPath());
[source_parent,source_name]=fileparts(source_run_dir);
if strcmp(source_name,'estimation')
    source_run_dir=source_parent;
    [~,source_name]=fileparts(source_run_dir);
end
parameters.run_dir=char(result.run_dir);
parameters.run_id=source_name;
parameters.estimation_config=struct();
for field={'p','seed','chain_id','burnin','thin','max_iterations','prior_weeks', ...
        'shrinkage','calibration_windows','calibration_window','ig_shape','kB', ...
        'num_chains','parallel_chains','convergence_mode','auto_extend','check_every', ...
        'max_retained_draws','rhat_threshold','min_ess','max_mcse_ratio','min_diagnostic_draws'}
    if isfield(cfg,field{1})
        parameters.estimation_config.(field{1})=cfg.(field{1});
    end
end
unique=char(java.util.UUID.randomUUID());
path=fullfile(parameter_dir,['parameters_' source_name '.mat']);
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

function [draws,ids,counts,dirs,seeds]=collect_parameter_draws(result,count,panel,priors,cfg)
dirs={result.run_dir};
if isfield(result,'chain_dirs'), dirs=cellstr(string(result.chain_dirs)); end
dirs=cellfun(@(p)char(java.io.File(p).getCanonicalPath()),dirs,'UniformOutput',false);
if isempty(dirs)||numel(unique(dirs))~=numel(dirs)
    error('dsc:ParametersIdentity','Chain directories must be nonempty and distinct.');
end
m=numel(panel.tickers); nr=numel(panel.pair_i);
V=zeros(m,m,count); sig2h=zeros(count,m); sig2r=zeros(count,nr);
draw_ids=zeros(count,1); iterations=zeros(count,1); offset=0;
ids=zeros(1,numel(dirs)); counts=ids; seeds=ids; target=[];
for c=1:numel(dirs)
    cp_path=fullfile(dirs{c},'checkpoint.mat');
    if ~isfile(cp_path), error('dsc:ParametersIdentity','Missing source chain checkpoint: %s.',cp_path); end
    loaded=load(cp_path,'checkpoint'); cp=loaded.checkpoint; identity=cp.identity;
    if ~isequal(identity.dates,panel.dates)||~isequaln(identity.returns,panel.returns)|| ...
            ~isequal(identity.mask,logical(panel.observation_mask)&isfinite(panel.returns))|| ...
            ~isequal(identity.tickers,panel.tickers)||~isequaln(identity.priors,priors)|| ...
            ~strcmp(identity.source_hash,panel.source_hash)|| ...
            ~identity.cfg.estimate_parameters||identity.cfg.burnin~=cfg.burnin||identity.cfg.thin~=cfg.thin
        error('dsc:ParametersIdentity','Source chains must share the saved estimation sample, priors and retention settings.');
    end
    ids(c)=identity.cfg.chain_id; seeds(c)=identity.cfg.seed+ids(c)-1;
    current=identity.cfg;
    for field={'seed','chain_id','max_iterations','max_seconds','resume','output_root'}
        if isfield(current,field{1}), current=rmfield(current,field{1}); end
    end
    if c==1, target=current;
    elseif ~isequaln(current,target), error('dsc:ParametersIdentity','Chains have different fixed configurations.'); end
    files=dir(fullfile(dirs{c},'posterior_chunk_*.mat'));
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
        if size(chunk.V,1)~=m||size(chunk.V,2)~=m||size(chunk.V,3)~=n|| ...
                ~isequal(size(chunk.sig2h),[n m])||~isequal(size(chunk.sig2r),[n nr])||offset+n>count
            error('dsc:ParametersDraws','Posterior chunk %s has invalid draw dimensions/count.',files(k).name);
        end
        rows=offset+(1:n); V(:,:,rows)=chunk.V; sig2h(rows,:)=chunk.sig2h; sig2r(rows,:)=chunk.sig2r;
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
draws=struct('V',V,'sig2h',sig2h,'sig2r',sig2r,'chain_id',draw_ids,'iteration',iterations);
validate_parameter_draws_local(draws,m,nr,count);
end

function estimate=parameter_draw_mean(draws)
estimate=struct('V',mean(draws.V,3),'sig2h',mean(draws.sig2h,1), ...
    'sig2r',mean(draws.sig2r,1));
end

function estimate=canonical_estimate(estimate)
estimate=struct('V',(double(estimate.V)+double(estimate.V)')/2, ...
    'sig2h',double(estimate.sig2h(:))','sig2r',double(estimate.sig2r(:))');
end

function tf=estimates_close(a,b)
tf=max(abs(a.V(:)-b.V(:)))<=1e-12*max(1,max(abs(b.V(:))))&& ...
    max(abs(a.sig2h(:)-b.sig2h(:)))<=1e-12*max(1,max(abs(b.sig2h(:))))&& ...
    max(abs(a.sig2r(:)-b.sig2r(:)))<=1e-12*max(1,max(abs(b.sig2r(:))));
end

function validate_parameter_draws_local(draws,m,nr,count)
if ~isstruct(draws)||~isscalar(draws)||~all(isfield(draws,{'V','sig2h','sig2r'}))|| ...
        ~isequal(size(draws.V),[m m count])||~isequal(size(draws.sig2h),[count m])|| ...
        ~isequal(size(draws.sig2r),[count nr])||any(~isfinite(draws.V(:)))|| ...
        any(~isfinite(draws.sig2h(:)))||any(draws.sig2h(:)<=0)|| ...
        any(~isfinite(draws.sig2r(:)))||any(draws.sig2r(:)<=0)
    error('dsc:ParametersValues','Parameter draws must be finite and positive with matching dimensions.');
end
for d=1:count
    V=(draws.V(:,:,d)+draws.V(:,:,d)')/2;
    if norm(draws.V(:,:,d)-draws.V(:,:,d)','fro')>1e-10*max(1,norm(draws.V(:,:,d),'fro'))
        error('dsc:ParametersValues','Parameter V draw %d is not symmetric.',d);
    end
    [~,bad]=chol(V);
    if bad, error('dsc:ParametersValues','Parameter V draw %d is not positive definite.',d); end
end
end

function path = dsc_save_parameters(result,panel,priors,cfg,parameter_dir,calibration_panel)
%DSC_SAVE_PARAMETERS Save retained parameter estimates for later smoothing.
%   Retained parameter draws are saved for reuse. Their posterior averages
%   are also stored as a plug-in summary. A checkpoint's last state
%   (including a warm-up state) is never a parameter estimate.
%   Writes a timestamped schema-versioned MAT file atomically, preserving
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
draws=collect_parameter_draws(result.run_dir,count,numel(panel.tickers),numel(panel.pair_i));
estimate=parameter_draw_mean(draws);
if ~isempty(result.parameter_estimate)&&~estimates_close(canonical_estimate(result.parameter_estimate),estimate)
    error('dsc:ParametersValues','Saved posterior chunks disagree with the sampler parameter estimate.');
end
parameters=struct();
parameters.schema_version=3;
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
parameters.estimation_status=result.status;
if isfield(result,'completed_iterations'), parameters.completed_iterations=result.completed_iterations; end
if isfield(result,'implementation'), parameters.implementation=result.implementation; end
if all(isfield(cfg,{'max_iterations','burnin','thin'}))
    parameters.requested_retained_draws=floor((cfg.max_iterations-cfg.burnin)/cfg.thin);
end
parameters.convergence_established=false;
parameters.run_dir=char(result.run_dir);
parameters.estimation_config=struct();
for field={'p','seed','chain_id','burnin','thin','max_iterations','prior_weeks', ...
        'shrinkage','calibration_windows','calibration_window','ig_shape','kB'}
    if isfield(cfg,field{1})
        parameters.estimation_config.(field{1})=cfg.(field{1});
    end
end
stamp=char(string(estimated,"yyyyMMdd'T'HHmmssSSS'Z'"));
unique=char(java.util.UUID.randomUUID());
path=fullfile(parameter_dir,['dsc_parameters_' stamp '_' unique '.mat']);
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

function draws=collect_parameter_draws(run_dir,count,m,nr)
files=dir(fullfile(run_dir,'posterior_chunk_*.mat'));
if isempty(files)
    error('dsc:ParametersDraws','No posterior chunks found in %s; cannot save parameter draws.',run_dir);
end
[~,order]=sort({files.name}); files=files(order);
V=zeros(m,m,count); sig2h=zeros(count,m); sig2r=zeros(count,nr); offset=0;
for k=1:numel(files)
    loaded=load(fullfile(files(k).folder,files(k).name),'chunk');
    if ~isfield(loaded,'chunk')||~isstruct(loaded.chunk)||~isscalar(loaded.chunk)
        error('dsc:ParametersDraws','Invalid posterior chunk: %s',files(k).name);
    end
    chunk=loaded.chunk; n=numel(chunk.iterations);
    if size(chunk.V,1)~=m||size(chunk.V,2)~=m||size(chunk.V,3)~=n|| ...
            ~isequal(size(chunk.sig2h),[n m])||~isequal(size(chunk.sig2r),[n nr])
        error('dsc:ParametersDraws','Posterior chunk %s has invalid parameter-draw dimensions.',files(k).name);
    end
    rows=offset+(1:n);
    if rows(end)>count
        error('dsc:ParametersDraws','Posterior chunks contain more parameter draws than the sampler summary.');
    end
    V(:,:,rows)=chunk.V; sig2h(rows,:)=chunk.sig2h; sig2r(rows,:)=chunk.sig2r;
    offset=offset+n;
end
if offset~=count
    error('dsc:ParametersDraws','Posterior chunks contain %d parameter draws, expected %d.',offset,count);
end
draws=struct('V',V,'sig2h',sig2h,'sig2r',sig2r);
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

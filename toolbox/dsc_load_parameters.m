function bundle = dsc_load_parameters(parameter_dir,panel,explicit_file)
%DSC_LOAD_PARAMETERS Load the latest compatible saved parameter estimate.
%   Latest means the saved UTC estimation time, never file modification time.
%   An explicit file bypasses discovery but receives the same validation.
%   PANEL is the full current panel, including its excluded calibration prefix.
%   New dates may be appended; calibration and likelihood histories, series
%   order and pair order must match. Original priors are returned unchanged.
if nargin<3, explicit_file=''; end
if ~is_text_scalar(parameter_dir)
    error('dsc:ParametersDirectory','parameter_dir must be a nonempty directory path.');
end
parameter_dir=char(parameter_dir);
if isempty(explicit_file)||(isstring(explicit_file)&&isscalar(explicit_file)&&strlength(explicit_file)==0)
    files=dir(fullfile(parameter_dir,'dsc_parameters_*.mat'));
    files=files(~[files.isdir]);
    if isempty(files)
        error('dsc:ParametersMissing','No dsc_parameters_*.mat files exist in %s. Estimate parameters first or choose an explicit parameter file.',parameter_dir);
    end
    latest=-Inf; bundle=[]; path='';
    for k=1:numel(files)
        candidate_path=fullfile(files(k).folder,files(k).name);
        candidate=read_candidate(candidate_path);
        if ~isfield(candidate,'estimated_at_posix')||~finite_scalar(candidate.estimated_at_posix)
            error('dsc:ParametersCorrupt','Cannot determine estimation time of parameter file %s. Choose an explicit valid file.',candidate_path);
        end
        % Directory names break exact timestamp ties deterministically.
        if candidate.estimated_at_posix>=latest
            latest=candidate.estimated_at_posix; bundle=candidate; path=candidate_path;
        end
    end
else
    if ~is_text_scalar(explicit_file)
        error('dsc:ParametersFile','explicit_file must be a scalar file path.');
    end
    path=char(explicit_file);
    if ~isfile(path), path=fullfile(parameter_dir,path); end
    if ~isfile(path), error('dsc:ParametersMissing','Parameter file does not exist: %s',path); end
    bundle=read_candidate(path);
end
path=char(java.io.File(path).getCanonicalPath());
validate_bundle(bundle,path);
validate_panel(bundle,panel,path);
if bundle.schema_version==3
    bundle.convergence=struct('status','not_checked','passed',false, ...
        'reason','Legacy schema 3 did not compute convergence diagnostics.');
end
bundle.parameter_estimate.sig2h=bundle.parameter_estimate.sig2h(:)';
bundle.parameter_estimate.sig2r=bundle.parameter_estimate.sig2r(:)';
bundle.parameter_draws.sig2h=double(bundle.parameter_draws.sig2h);
bundle.parameter_draws.sig2r=double(bundle.parameter_draws.sig2r);
for d=1:size(bundle.parameter_draws.V,3)
    bundle.parameter_draws.V(:,:,d)=(double(bundle.parameter_draws.V(:,:,d))+double(bundle.parameter_draws.V(:,:,d))')/2;
end
bundle.path=path;
end

function value=read_candidate(path)
try
    loaded=load(path,'parameters');
catch problem
    error('dsc:ParametersCorrupt','Cannot read parameter file %s: %s',path,problem.message);
end
if ~isfield(loaded,'parameters')||~isstruct(loaded.parameters)||~isscalar(loaded.parameters)
    error('dsc:ParametersCorrupt','File %s has no scalar parameters bundle. Checkpoints and legacy MAT files cannot be reused as parameter estimates.',path);
end
value=loaded.parameters;
end

function validate_bundle(b,path)
if isfield(b,'schema_version')&&isequal(b.schema_version,1)
    error('dsc:ParametersLegacyCalibration', ...
        ['Parameter file %s uses schema 1, whose model included calibration weeks in the likelihood. ' ...
         'Re-estimate parameters with the calibration weeks excluded before smoothing.'],path);
end
if isfield(b,'schema_version')&&isequal(b.schema_version,2)
    error('dsc:ParametersLegacyDraws', ...
        ['Parameter file %s stores posterior means only. Re-estimate parameters to save retained ' ...
         'parameter draws before smoothing with this interface.'],path);
end
required={'schema_version','model','estimator','estimated_at','estimated_at_posix', ...
    'estimation_start','estimation_end','tickers','pair_i','pair_j','source_hash', ...
    'calibration_dates','calibration_returns','calibration_mask','calibration_weeks','calibration_start','calibration_end', ...
    'training_dates','training_returns','training_mask','priors','parameter_draws','parameter_estimate', ...
    'retained_draws','convergence_established','run_dir','estimation_config'};
if ~all(isfield(b,required))
    error('dsc:ParametersSchema','Parameter file %s is missing required metadata.',path);
end
if ~ismember(b.schema_version,[3 4])||~strcmp(b.model,'joint_dsc_p0')||~strcmp(b.estimator,'posterior_draws')
    error('dsc:ParametersSchema','Unsupported parameter schema, model or estimator in %s.',path);
end
if ~finite_scalar(b.estimated_at_posix)||~is_text_scalar(b.estimated_at)|| ...
        isempty(regexp(char(b.estimated_at),'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z$','once'))
    error('dsc:ParametersSchema','Invalid UTC estimation timestamp in %s.',path);
end
try
    date=datetime(b.estimated_at,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC');
catch
    error('dsc:ParametersSchema','Invalid UTC estimation timestamp in %s.',path);
end
if isnat(date)||abs(posixtime(date)-b.estimated_at_posix)>5e-4
    error('dsc:ParametersSchema','Inconsistent UTC estimation timestamps in %s.',path);
end
if ~weekly_dates(b.training_dates)||~weekly_dates(b.calibration_dates)|| ...
        b.training_dates(1)~=b.calibration_dates(end)+calweeks(1)
    error('dsc:ParametersCalibration','Calibration and training dates must form consecutive nonoverlapping weekly blocks in %s.',path);
end
if ~isequal(b.estimation_start,b.training_dates(1))||~isequal(b.estimation_end,b.training_dates(end))
    error('dsc:ParametersSchema','Estimation dates do not match the saved training history in %s.',path);
end
if ~finite_scalar(b.calibration_weeks)||b.calibration_weeks<1||b.calibration_weeks~=floor(b.calibration_weeks)|| ...
        b.calibration_weeks~=numel(b.calibration_dates)|| ...
        ~isequal(b.calibration_start,b.calibration_dates(1))||~isequal(b.calibration_end,b.calibration_dates(end))
    error('dsc:ParametersCalibration','Calibration dates and week count disagree in %s.',path);
end
if ~(isstring(b.tickers)||iscellstr(b.tickers))||~isvector(b.tickers)
    error('dsc:ParametersSchema','Invalid series names in %s.',path);
end
names=string(b.tickers(:)); m=numel(names); T=numel(b.training_dates); nr=m*(m-1)/2;
if m<2||any(ismissing(names))||any(strlength(names)==0)||numel(unique(names))~=m
    error('dsc:ParametersSchema','Invalid series names in %s.',path);
end
[ii,jj]=find(tril(true(m),-1));
if ~isequal(b.pair_i(:),ii)||~isequal(b.pair_j(:),jj)
    error('dsc:ParametersSchema','Invalid column-major pair ordering in %s.',path);
end
if ~valid_data_block(b.training_returns,b.training_mask,T,m)|| ...
        ~valid_data_block(b.calibration_returns,b.calibration_mask,b.calibration_weeks,m)
    error('dsc:ParametersSchema','Invalid calibration/training returns or observation masks in %s.',path);
end
if ~is_text_scalar(b.source_hash)||~is_text_scalar(b.run_dir)|| ...
        ~isstruct(b.estimation_config)||~isscalar(b.estimation_config)|| ...
        ~isfield(b.estimation_config,'p')||~isequal(b.estimation_config.p,0)|| ...
        ~isfield(b.estimation_config,'prior_weeks')||~isequal(b.estimation_config.prior_weeks,b.calibration_weeks)
    error('dsc:ParametersSchema','Invalid source or configuration metadata in %s.',path);
end
if ~finite_scalar(b.retained_draws)||b.retained_draws<2||b.retained_draws~=floor(b.retained_draws)
    error('dsc:ParametersDraws','Parameter file %s does not contain an estimate from at least two retained draws.',path);
end
if b.schema_version==3&&~isequal(b.convergence_established,false)
    error('dsc:ParametersSchema','Schema 3 does not establish convergence; invalid convergence metadata in %s.',path);
end
v=b.parameter_estimate;
if ~isstruct(v)||~isscalar(v)||~all(isfield(v,{'V','sig2h','sig2r'}))|| ...
        ~positive_covariance(v.V,m)||~positive_vector(v.sig2h,m)||~positive_vector(v.sig2r,nr)
    error('dsc:ParametersValues','Parameter file %s must contain a positive-definite V and positive finite innovation variances.',path);
end
draws=b.parameter_draws;
if ~valid_parameter_draws(draws,m,nr,b.retained_draws)
    error('dsc:ParametersValues','Parameter file %s must contain retained positive parameter draws.',path);
end
if b.schema_version==4
    validate_chain_provenance(b,path);
end
draw_mean=struct('V',mean(draws.V,3),'sig2h',mean(draws.sig2h,1),'sig2r',mean(draws.sig2r,1));
if max(abs(draw_mean.V(:)-v.V(:)))>1e-10|| ...
        max(abs(draw_mean.sig2h(:)-v.sig2h(:)))>1e-10|| ...
        max(abs(draw_mean.sig2r(:)-v.sig2r(:)))>1e-10
    error('dsc:ParametersValues','Parameter means must equal the average of retained draws in %s.',path);
end
p=b.priors;
prior_fields={'Bbar','VBbar','nuB','V0B','mh0','mr0','h0_scale','r0_scale', ...
    'sig2h_mean','sig2r_mean','ig_shape','ig_scale_h','ig_scale_r','source_hash','tickers'};
if ~isstruct(p)||~isscalar(p)||~all(isfield(p,prior_fields))|| ...
        ~finite_vector(p.Bbar,m)||~positive_covariance(p.VBbar,m)|| ...
        ~positive_covariance(p.V0B,m)||~finite_scalar(p.nuB)||p.nuB<=m+1|| ...
        ~finite_vector(p.mh0,m)||~finite_vector(p.mr0,nr)|| ...
        ~positive_vector(p.h0_scale,1)||~positive_vector(p.r0_scale,1)|| ...
        ~positive_vector(p.sig2h_mean,1)||~positive_vector(p.sig2r_mean,1)|| ...
        ~finite_scalar(p.ig_shape)||p.ig_shape<=1|| ...
        ~positive_vector(p.ig_scale_h,1)||~positive_vector(p.ig_scale_r,1)|| ...
        ~isequal(p.tickers,b.tickers)||~isequal(p.source_hash,b.source_hash)
    error('dsc:ParametersPriors','Saved priors are invalid or inconsistent with the estimation sample in %s.',path);
end
if (isfield(p,'initial_weeks')&&~isequal(p.initial_weeks,b.calibration_weeks))|| ...
        (isfield(p,'initial_dates')&&~isequal(p.initial_dates(:),reshape(b.calibration_dates([1 end]),[],1)))
    error('dsc:ParametersCalibration','Initial-prior window does not match the excluded calibration dates in %s.',path);
end
end

function validate_chain_provenance(b,path)
fields={'chain_ids','chain_seeds','chain_dirs','retained_draws_per_chain','convergence'};
if ~all(isfield(b,fields))||~all(isfield(b.parameter_draws,{'chain_id','iteration'}))
    error('dsc:ParametersSchema','Missing chain or diagnostic provenance in %s.',path);
end
ids=b.chain_ids(:); counts=b.retained_draws_per_chain(:); seeds=b.chain_seeds(:);
if isempty(ids)||~finite_vector(ids,numel(ids))||any(ids<1)||any(ids~=floor(ids))|| ...
        numel(unique(ids))~=numel(ids)||~finite_vector(counts,numel(ids))|| ...
        any(counts<1)||any(counts~=floor(counts))||sum(counts)~=b.retained_draws|| ...
        ~finite_vector(seeds,numel(ids))||any(seeds<0)||any(seeds>2^32-1)|| ...
        any(seeds~=floor(seeds))||numel(unique(seeds))~=numel(seeds)|| ...
        numel(b.chain_dirs)~=numel(ids)||numel(unique(string(b.chain_dirs)))~=numel(ids)
    error('dsc:ParametersIdentity','Invalid original chain identities/counts in %s.',path);
end
cfg=b.estimation_config;
if ~all(isfield(cfg,{'seed','burnin','thin'}))|| ...
        ~finite_scalar(cfg.seed)||cfg.seed<0||cfg.seed>2^32-1||cfg.seed~=floor(cfg.seed)|| ...
        ~finite_scalar(cfg.burnin)||cfg.burnin<0||cfg.burnin~=floor(cfg.burnin)|| ...
        ~finite_scalar(cfg.thin)||cfg.thin<1||cfg.thin~=floor(cfg.thin)|| ...
        any(seeds~=cfg.seed+ids-1)
    error('dsc:ParametersIdentity','Invalid chain seeds, warm-up or thinning configuration in %s.',path);
end
if isfield(cfg,'chain_id')&&(~finite_scalar(cfg.chain_id)||cfg.chain_id<1|| ...
        cfg.chain_id~=floor(cfg.chain_id)||cfg.chain_id>2^32)
    error('dsc:ParametersIdentity','Invalid base chain identity in %s.',path);
end
if isfield(cfg,'num_chains')
    if ~finite_scalar(cfg.num_chains)||cfg.num_chains<1||cfg.num_chains~=floor(cfg.num_chains)
        error('dsc:ParametersIdentity','Invalid configured chain count in %s.',path);
    end
    % Explicit coordinator metadata promises consecutive IDs from its base.
    % Legacy pooling may combine an arbitrary subset of independent chains.
    if cfg.num_chains>1&&(~isfield(cfg,'chain_id')||cfg.num_chains~=numel(ids)|| ...
            ~isequal(ids,cfg.chain_id+(0:cfg.num_chains-1)'))
        error('dsc:ParametersIdentity','Configured chain base/count disagree with original chain identities in %s.',path);
    end
end
draw_ids=b.parameter_draws.chain_id(:); iterations=b.parameter_draws.iteration(:);
if ~finite_vector(draw_ids,b.retained_draws)||~all(ismember(draw_ids,ids))|| ...
        ~finite_vector(iterations,b.retained_draws)||any(iterations~=floor(iterations))|| ...
        ~all(isfield(b.estimation_config,{'burnin','thin'}))
    error('dsc:ParametersIdentity','Invalid per-draw provenance in %s.',path);
end
for c=1:numel(ids)
    ix=iterations(draw_ids==ids(c));
    if numel(ix)~=counts(c)||any(diff([b.estimation_config.burnin;ix])~=b.estimation_config.thin)
        error('dsc:ParametersIdentity','Invalid per-chain retention sequence in %s.',path);
    end
end
m=numel(b.tickers); nr=m*(m-1)/2;
expected_quantities=numel(b.training_dates)*(m+nr);
dsc_validate_convergence(b.convergence,ids,counts,expected_quantities);
if ~isequal(b.convergence_established,b.convergence.passed)
    error('dsc:ConvergenceMetadata','Parameter convergence flag disagrees with its diagnostic report.');
end
end

function validate_panel(b,panel,path)
fields={'dates','returns','observation_mask','tickers','pair_i','pair_j','source_hash'};
if ~isstruct(panel)||~isscalar(panel)||~all(isfield(panel,fields))
    error('dsc:ParametersPanel','A complete weekly panel is required to validate %s.',path);
end
if ~isequal(string(panel.tickers(:)),string(b.tickers(:)))|| ...
        ~isequal(panel.pair_i(:),b.pair_i(:))||~isequal(panel.pair_j(:),b.pair_j(:))
    error('dsc:ParametersIdentity','Series names/order or pair order differ from the model estimated at %s.',b.estimated_at);
end
T=numel(panel.dates); m=numel(b.tickers); n=numel(b.training_dates); c=b.calibration_weeks;
if ~isdatetime(panel.dates)||~isvector(panel.dates)||T<c+n||any(isnat(panel.dates))|| ...
        any(diff(panel.dates)<=seconds(0))||~isnumeric(panel.returns)||~isreal(panel.returns)|| ...
        ~isequal(size(panel.returns),[T m])||~isequal(size(panel.observation_mask),[T m])|| ...
        ~(islogical(panel.observation_mask)||(isnumeric(panel.observation_mask)&& ...
        all(ismember(panel.observation_mask(:),[0 1]))))|| ...
        any(~isfinite(panel.returns(logical(panel.observation_mask))))||~is_text_scalar(panel.source_hash)
    error('dsc:ParametersPanel','New data must be a valid full panel containing both calibration and estimation histories.');
end
if ~isequal(panel.dates(1:c),b.calibration_dates)|| ...
        ~isequaln(panel.returns(1:c,:),b.calibration_returns)|| ...
        ~isequal(logical(panel.observation_mask(1:c,:)),b.calibration_mask)|| ...
        ~isequal(panel.dates(c+(1:n)),b.training_dates)|| ...
        ~isequaln(panel.returns(c+(1:n),:),b.training_returns)|| ...
        ~isequal(logical(panel.observation_mask(c+(1:n),:)),b.training_mask)
    error('dsc:ParametersHistory','Historical dates, returns or observation mask have changed since parameter estimation at %s. Re-estimate parameters for revised history.',b.estimated_at);
end
if ~weekly_dates(panel.dates)
    error('dsc:ParametersPanel','The full current panel must retain consecutive weekly dates.');
end
% A CSV hash changes when dates are appended. Exact history matching above,
% not equality of full-file hashes, prevents accidental reuse on revised data.
end

function tf=weekly_dates(dates)
tf=isdatetime(dates)&&isvector(dates)&&~isempty(dates)&&~any(isnat(dates))&& ...
    all(dates(2:end)==dates(1:end-1)+calweeks(1));
end

function tf=valid_data_block(returns,mask,T,m)
tf=isnumeric(returns)&&isreal(returns)&&isequal(size(returns),[T m])&& ...
    ~any(isinf(returns(:)))&&islogical(mask)&&isequal(size(mask),[T m])&& ...
    all(isfinite(returns(mask)));
end

function tf=finite_scalar(value)
tf=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value);
end

function tf=finite_vector(value,n)
tf=isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==n&&all(isfinite(value(:)));
end

function tf=positive_vector(value,n)
tf=finite_vector(value,n)&&all(value(:)>0);
end

function tf=positive_covariance(value,n)
tf=isnumeric(value)&&isreal(value)&&isequal(size(value),[n n])&&all(isfinite(value(:)));
if ~tf, return; end
tf=norm(value-value','fro')<=1e-10*max(1,norm(value,'fro'));
if tf, [~,bad]=chol((value+value')/2); tf=bad==0; end
end

function tf=valid_parameter_draws(draws,m,nr,n)
tf=isstruct(draws)&&isscalar(draws)&&all(isfield(draws,{'V','sig2h','sig2r'}))&& ...
    isnumeric(draws.V)&&isreal(draws.V)&&isequal(size(draws.V),[m m n])&& ...
    isnumeric(draws.sig2h)&&isreal(draws.sig2h)&&isequal(size(draws.sig2h),[n m])&& ...
    isnumeric(draws.sig2r)&&isreal(draws.sig2r)&&isequal(size(draws.sig2r),[n nr])&& ...
    all(isfinite(draws.V(:)))&&all(isfinite(draws.sig2h(:)))&&all(draws.sig2h(:)>0)&& ...
    all(isfinite(draws.sig2r(:)))&&all(draws.sig2r(:)>0);
if ~tf, return; end
for d=1:n
    if ~positive_covariance(draws.V(:,:,d),m), tf=false; return; end
end
end

function tf=is_text_scalar(value)
tf=(ischar(value)&&isrow(value)&&~isempty(value))|| ...
    (isstring(value)&&isscalar(value)&&~ismissing(value)&&strlength(value)>0);
end

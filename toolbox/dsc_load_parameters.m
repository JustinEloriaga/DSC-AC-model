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
    files=dir(fullfile(parameter_dir,'*_parameters.mat'));
    files=files(~[files.isdir]);
    files=files(~contains({files.name},'pending'));
    if isempty(files)
        error('dsc:ParametersMissing','No prefixed *_parameters.mat files exist in %s. Estimate parameters first or choose an explicit parameter file.',parameter_dir);
    end
    latest=-Inf; bundle=[]; path='';
    for k=1:numel(files)
        candidate_path=fullfile(files(k).folder,files(k).name);
        candidate=read_candidate(candidate_path);
        if isfield(candidate,'storage')&&strcmp(candidate.storage,'conditional_parameter_summary_v1'), continue; end
        % Historical dynamic-mean bundles remain untouched and cannot be
        % selected for this constant-mean model. Other panels are separate fits.
        if legacy_mean_model(candidate), continue; end
        if isfield(candidate,'tickers')&&isfield(panel,'tickers')&& ...
                ~isequal(string(candidate.tickers(:)),string(panel.tickers(:))), continue; end
        try
            candidate=expand_report(candidate,panel,candidate_path);
        catch problem
            % Discovery may encounter earlier fits of a revised history.
            % Only a history mismatch makes a healthy report ineligible;
            % malformed reports and missing artifacts must remain visible.
            if strcmp(problem.identifier,'dsc:ParametersHistory'), continue; end
            rethrow(problem);
        end
        validate_bundle(candidate,candidate_path);
        if ~isfield(candidate,'estimated_at_posix')||~finite_scalar(candidate.estimated_at_posix)
            error('dsc:ParametersCorrupt','Cannot determine estimation time of parameter file %s. Choose an explicit valid file.',candidate_path);
        end
        % Directory names break exact timestamp ties deterministically.
        if candidate.estimated_at_posix>=latest
            latest=candidate.estimated_at_posix; bundle=candidate; path=candidate_path;
        end
    end
    if isempty(bundle)
        error('dsc:ParametersMissing','No compatible constant-mean parameter file exists for the current series in %s. Estimate this model first.',parameter_dir);
    end
else
    if ~is_text_scalar(explicit_file)
        error('dsc:ParametersFile','explicit_file must be a scalar file path.');
    end
    path=char(explicit_file);
    if ~isfile(path), path=fullfile(parameter_dir,path); end
    if ~isfile(path), error('dsc:ParametersMissing','Parameter file does not exist: %s',path); end
    bundle=read_candidate(path);
    if isfield(bundle,'storage')&&strcmp(bundle.storage,'conditional_parameter_summary_v1')
        error('dsc:ParametersFixedMode','File %s reports fixed-parameter smoothing; reuse its source_parameter_file instead.',path);
    end
    bundle=expand_report(bundle,panel,path);
end
path=char(java.io.File(path).getCanonicalPath());
validate_bundle(bundle,path);
validate_panel(bundle,panel,path);
bundle.parameter_estimate.mu=double(bundle.parameter_estimate.mu(:)');
bundle.parameter_estimate.sig2h=double(bundle.parameter_estimate.sig2h(:)');
bundle.parameter_estimate.sig2r=double(bundle.parameter_estimate.sig2r(:)');
bundle.parameter_estimate.delta=bundle.parameter_estimate.sig2h;
bundle.parameter_estimate.gamma=bundle.parameter_estimate.sig2r;
bundle.parameter_draws.mu=double(bundle.parameter_draws.mu);
bundle.parameter_draws.sig2h=double(bundle.parameter_draws.sig2h);
bundle.parameter_draws.sig2r=double(bundle.parameter_draws.sig2r);
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
value=without_legacy_chain_statistics(loaded.parameters);
end

function b=expand_report(b,panel,path)
% Final reports store parameter draws once, in their original joint chunks.
if ~isfield(b,'storage')||~strcmp(b.storage,'parameter_report_v1'), return; end
required={'manifest_file','parameter_chunk_files', ...
    'calibration_sample_sha256', ...
    'training_sample_sha256','calibration_weeks','training_dates','tickers', ...
    'pair_i','pair_j','retained_draws','chain_ids','retained_draws_per_chain'};
if ~all(isfield(b,required))
    error('dsc:ParametersSchema','Incomplete parameter report in %s.',path);
end
if ~isstruct(panel)||~all(isfield(panel,{'dates','returns','observation_mask','tickers','pair_i','pair_j'}))|| ...
        ~isequal(string(panel.tickers(:)),string(b.tickers(:)))|| ...
        ~isequal(panel.pair_i(:),b.pair_i(:))||~isequal(panel.pair_j(:),b.pair_j(:))
    error('dsc:ParametersIdentity','A matching full weekly panel is required to expand %s.',path);
end
c=b.calibration_weeks; n=numel(b.training_dates);
if ~finite_scalar(c)||c<1||c~=floor(c)||numel(panel.dates)<c+n
    error('dsc:ParametersHistory','Current panel is shorter than the saved estimation history in %s.',path);
end
cal=sample_block(panel,1:c); train=sample_block(panel,c+(1:n));
if ~strcmp(dsc_sample_fingerprint(cal),b.calibration_sample_sha256)|| ...
        ~strcmp(dsc_sample_fingerprint(train),b.training_sample_sha256)
    error('dsc:ParametersHistory','Calibration or estimation history has changed since %s.',path);
end
b.calibration_returns=cal.returns; b.calibration_mask=logical(cal.observation_mask);
b.training_returns=train.returns; b.training_mask=logical(train.observation_mask);
if ~isfield(b,'priors')
    % Legacy reports stored their original priors in a manifest section.
    % Never reconstruct or recalibrate missing saved priors from current data.
    try
        if ~isfield(b,'prior_manifest_section')
            error('dsc:ParametersSchema','Saved calibrated priors are missing.');
        end
        manifest=jsondecode(fileread(b.manifest_file));
        section=char(b.prior_manifest_section);
        if strcmp(section,'priors')
            p=manifest.priors;
        elseif strcmp(section,'parameter_estimation')
            p=manifest.parameter_estimation.priors;
        else
            error('dsc:ParametersSchema','Unknown prior manifest section.');
        end
        p.tickers=b.tickers;
        for name={'mu0','mh0','mr0'}, p.(name{1})=p.(name{1})(:); end
        for name={'initial_dates','calibration_dates','likelihood_start_date'}
            if isfield(p,name{1})
                p.(name{1})=datetime(string(p.(name{1})),'InputFormat','yyyy-MM-dd');
                p.(name{1})=p.(name{1})(:);
            end
        end
        b.priors=p;
    catch problem
        error('dsc:ParametersCorrupt','Cannot expand parameter report %s: %s',path,problem.message);
    end
end
b.parameter_draws=report_draws(b,path);
end

function block=sample_block(panel,rows)
block=struct('dates',panel.dates(rows),'returns',panel.returns(rows,:), ...
    'observation_mask',panel.observation_mask(rows,:),'tickers',panel.tickers, ...
    'pair_i',panel.pair_i,'pair_j',panel.pair_j);
end

function draws=report_draws(b,path)
m=numel(b.tickers); q=numel(b.pair_i); N=b.retained_draws; T=numel(b.training_dates);
if ~finite_scalar(N)||N<2||N~=floor(N)|| ...
        ~(iscellstr(b.parameter_chunk_files)||isstring(b.parameter_chunk_files))
    error('dsc:ParametersDraws','Invalid retained draw metadata in %s.',path);
end
files=cellstr(b.parameter_chunk_files);
if isempty(files)||numel(unique(files))~=numel(files)
    error('dsc:ParametersDraws','Missing or repeated chunk references in %s.',path);
end
draws=struct('mu',zeros(N,m),'sig2h',zeros(N,m),'sig2r',zeros(N,q), ...
    'chain_id',zeros(N,1),'iteration',zeros(N,1));
offset=0;
for k=1:numel(files)
    try
        loaded=load(files{k},'chunk'); ch=loaded.chunk;
    catch problem
        error('dsc:ParametersDraws','Cannot read joint chunk %s: %s',files{k},problem.message);
    end
    required={'schema','mu','sig2h','sig2r','h','P_pairs','iterations','chain_id', ...
        'dates','tickers','pair_i','pair_j'};
    if ~isstruct(ch)||~isscalar(ch)||~all(isfield(ch,required))
        error('dsc:ParametersDraws','Incomplete joint posterior chunk %s.',files{k});
    end
    count=numel(ch.iterations);
    if ~strcmp(ch.schema,'constant_mu_v1')||count<1||offset+count>N|| ...
            ~finite_scalar(ch.chain_id)||~ismember(ch.chain_id,b.chain_ids)|| ...
            ~isequal(ch.dates,b.training_dates)|| ...
            ~isequal(string(ch.tickers(:)),string(b.tickers(:)))|| ...
            ~isequal(ch.pair_i(:),b.pair_i(:))||~isequal(ch.pair_j(:),b.pair_j(:))|| ...
            ~isequal(size(ch.mu),[count m])||~isequal(size(ch.sig2h),[count m])|| ...
            ~isequal(size(ch.sig2r),[count q])|| ...
            size(ch.h,1)~=T||size(ch.h,2)~=m||size(ch.h,3)~=count|| ...
            size(ch.P_pairs,1)~=T||size(ch.P_pairs,2)~=q||size(ch.P_pairs,3)~=count
        error('dsc:ParametersDraws','Invalid joint draw dimensions or identity in %s.',files{k});
    end
    rows=offset+(1:count);
    draws.mu(rows,:)=ch.mu; draws.sig2h(rows,:)=ch.sig2h; draws.sig2r(rows,:)=ch.sig2r;
    draws.chain_id(rows)=ch.chain_id; draws.iteration(rows)=ch.iterations(:);
    offset=offset+count;
end
if offset~=N
    error('dsc:ParametersDraws','Joint chunk count differs from report %s.',path);
end
end

function validate_bundle(b,path)
if legacy_mean_model(b)
    error('dsc:ParametersMeanModel', ...
        'Parameter file %s belongs to the historical dynamic-mean model. Re-estimate the constant-mean model before smoothing.',path);
end
required={'schema_version','model','estimator','estimated_at','estimated_at_posix', ...
    'estimation_start','estimation_end','tickers','pair_i','pair_j','source_hash', ...
    'calibration_dates','calibration_returns','calibration_mask','calibration_weeks','calibration_start','calibration_end', ...
    'training_dates','training_returns','training_mask','priors','parameter_draws','parameter_estimate', ...
    'retained_draws','run_dir','file_prefix','estimation_config', ...
    'mean_model','innovation_variance_aliases'};
if ~all(isfield(b,required))
    error('dsc:ParametersSchema','Parameter file %s is missing required metadata.',path);
end
if ~isequal(b.schema_version,5)||~strcmp(b.model,'joint_dsc_constant_mean')|| ...
        ~strcmp(b.mean_model,'constant')||~strcmp(b.estimator,'posterior_draws')
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
        ~isfield(b.estimation_config,'prior_weeks')||~isequal(b.estimation_config.prior_weeks,b.calibration_weeks)
    error('dsc:ParametersSchema','Invalid source or configuration metadata in %s.',path);
end
if ~((ischar(b.file_prefix)&&isrow(b.file_prefix))||isempty(b.file_prefix))|| ...
        (~isempty(b.file_prefix)&&isempty(regexp(b.file_prefix,'^[A-Za-z0-9_-]+$','once')))
    error('dsc:ParametersSchema','Invalid run or filename prefix metadata in %s.',path);
end
[run_parent,run_name]=fileparts(char(b.run_dir));
if strcmp(run_name,'estimation'), [~,run_name]=fileparts(run_parent); end
if ~isempty(b.file_prefix), run_name=regexprep(b.file_prefix,'_+$',''); end
if ~finite_scalar(b.retained_draws)||b.retained_draws<2||b.retained_draws~=floor(b.retained_draws)
    error('dsc:ParametersDraws','Parameter file %s does not contain an estimate from at least two retained draws.',path);
end
v=b.parameter_estimate;
if ~isstruct(v)||~isscalar(v)||~all(isfield(v,{'mu','sig2h','sig2r','delta','gamma'}))|| ...
        isfield(v,'B')||isfield(v,'V')|| ...
        ~finite_vector(v.mu,m)||~positive_vector(v.sig2h,m)||~positive_vector(v.sig2r,nr)|| ...
        ~isequal(v.delta(:),v.sig2h(:))||~isequal(v.gamma(:),v.sig2r(:))|| ...
        ~isequal(b.innovation_variance_aliases,struct('delta','sig2h','gamma','sig2r'))
    error('dsc:ParametersValues','Parameter file %s must contain a finite constant mean and positive innovation variances with matching aliases.',path);
end
draws=b.parameter_draws;
if ~valid_parameter_draws(draws,m,nr,b.retained_draws)
    error('dsc:ParametersValues','Parameter file %s must contain retained finite mean and positive innovation-variance draws.',path);
end
validate_chain_provenance(b,path);
draw_mean=struct('mu',mean(draws.mu,1),'sig2h',mean(draws.sig2h,1),'sig2r',mean(draws.sig2r,1));
if max(abs(draw_mean.mu(:)-v.mu(:)))>1e-10|| ...
        max(abs(draw_mean.sig2h(:)-v.sig2h(:)))>1e-10|| ...
        max(abs(draw_mean.sig2r(:)-v.sig2r(:)))>1e-10
    error('dsc:ParametersValues','Parameter means must equal the average of retained draws in %s.',path);
end
p=b.priors;
prior_fields={'mu0','mu_cov0','mh0','mr0','h0_scale','r0_scale', ...
    'sig2h_mean','sig2r_mean','ig_shape','ig_scale_h','ig_scale_r','source_hash','tickers'};
if ~isstruct(p)||~isscalar(p)||~all(isfield(p,prior_fields))|| ...
        ~finite_vector(p.mu0,m)||~positive_covariance(p.mu_cov0,m)|| ...
        any(isfield(p,{'Bbar','VBbar','V0B','nuB'}))|| ...
        ~finite_vector(p.mh0,m)||~finite_vector(p.mr0,nr)|| ...
        ~positive_vector(p.h0_scale,1)||~positive_vector(p.r0_scale,1)|| ...
        ~positive_prior_vector(p.sig2h_mean,m)||~positive_prior_vector(p.sig2r_mean,nr)|| ...
        ~finite_scalar(p.ig_shape)||p.ig_shape<=1|| ...
        ~positive_prior_vector(p.ig_scale_h,m)||~positive_prior_vector(p.ig_scale_r,nr)|| ...
        ~isequal(p.tickers,b.tickers)||~isequal(p.source_hash,b.source_hash)
    error('dsc:ParametersPriors','Saved priors are invalid or inconsistent with the estimation sample in %s.',path);
end
if (isfield(p,'initial_weeks')&&~isequal(p.initial_weeks,b.calibration_weeks))|| ...
        (isfield(p,'initial_dates')&&~isequal(p.initial_dates(:),reshape(b.calibration_dates([1 end]),[],1)))
    error('dsc:ParametersCalibration','Initial-prior window does not match the excluded calibration dates in %s.',path);
end
end

function validate_chain_provenance(b,path)
fields={'chain_ids','chain_seeds','chain_dirs','chain_file_prefixes','retained_draws_per_chain'};
if ~all(isfield(b,fields))||~all(isfield(b.parameter_draws,{'chain_id','iteration'}))
    error('dsc:ParametersSchema','Missing chain or draw provenance in %s.',path);
end
if ~(iscellstr(b.chain_dirs)||(isstring(b.chain_dirs)&&isvector(b.chain_dirs)))|| ...
        ~(iscellstr(b.chain_file_prefixes)||(isstring(b.chain_file_prefixes)&&isvector(b.chain_file_prefixes)))|| ...
        any(ismissing(string(b.chain_dirs)))||any(strlength(string(b.chain_dirs))==0)|| ...
        any(ismissing(string(b.chain_file_prefixes)))
    error('dsc:ParametersIdentity','Invalid source chain directory/prefix metadata in %s.',path);
end
for prefix=string(b.chain_file_prefixes(:))'
    if strlength(prefix)>0&&isempty(regexp(char(prefix),'^[A-Za-z0-9_-]+$','once'))
        error('dsc:ParametersIdentity','Invalid source chain filename prefix in %s.',path);
    end
end
ids=b.chain_ids(:); counts=b.retained_draws_per_chain(:); seeds=b.chain_seeds(:);
if isempty(ids)||~finite_vector(ids,numel(ids))||any(ids<1)||any(ids~=floor(ids))|| ...
        numel(unique(ids))~=numel(ids)||~finite_vector(counts,numel(ids))|| ...
        any(counts<1)||any(counts~=floor(counts))||sum(counts)~=b.retained_draws|| ...
        ~finite_vector(seeds,numel(ids))||any(seeds<0)||any(seeds>2^32-1)|| ...
        any(seeds~=floor(seeds))||numel(unique(seeds))~=numel(seeds)|| ...
        numel(b.chain_dirs)~=numel(ids)||numel(b.chain_file_prefixes)~=numel(ids)|| ...
        numel(unique(strcat(string(b.chain_dirs(:)),filesep,string(b.chain_file_prefixes(:)))))~=numel(ids)
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
end

function b=without_legacy_chain_statistics(b)
% Older parameter files remain reusable without their diagnostic sidecars.
b=remove_fields(b,{'convergence','convergence_established','convergence_file','convergence_section'});
if isfield(b,'estimation_config')
    b.estimation_config=remove_fields(b.estimation_config, ...
        {'convergence_mode','auto_extend','check_every','max_retained_draws', ...
        'rhat_threshold','min_ess','max_mcse_ratio','min_diagnostic_draws'});
end
if isfield(b,'parameter_summary_metadata')
    b.parameter_summary_metadata=remove_fields(b.parameter_summary_metadata,{'convergence_established'});
    if isfield(b.parameter_summary_metadata,'posterior_std_convention')
        b.parameter_summary_metadata.posterior_std_convention='Sample standard deviation with n-1 denominator.';
    end
end
end

function value=remove_fields(value,names)
if isstruct(value)
    present=intersect(fieldnames(value),names);
    if ~isempty(present), value=rmfield(value,present); end
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

function tf=positive_prior_vector(value,n)
% A scalar is the historical shared prior; preserve its stored representation.
tf=positive_vector(value,n)||positive_vector(value,1);
end

function tf=positive_covariance(value,n)
tf=isnumeric(value)&&isreal(value)&&isequal(size(value),[n n])&&all(isfinite(value(:)));
if ~tf, return; end
tf=norm(value-value','fro')<=1e-10*max(1,norm(value,'fro'));
if tf, [~,bad]=chol((value+value')/2); tf=bad==0; end
end

function tf=valid_parameter_draws(draws,m,nr,n)
tf=isstruct(draws)&&isscalar(draws)&&all(isfield(draws,{'mu','sig2h','sig2r'}))&& ...
    ~isfield(draws,'B')&&~isfield(draws,'V')&& ...
    isnumeric(draws.mu)&&isreal(draws.mu)&&isequal(size(draws.mu),[n m])&& ...
    isnumeric(draws.sig2h)&&isreal(draws.sig2h)&&isequal(size(draws.sig2h),[n m])&& ...
    isnumeric(draws.sig2r)&&isreal(draws.sig2r)&&isequal(size(draws.sig2r),[n nr])&& ...
    all(isfinite(draws.mu(:)))&&all(isfinite(draws.sig2h(:)))&&all(draws.sig2h(:)>0)&& ...
    all(isfinite(draws.sig2r(:)))&&all(draws.sig2r(:)>0);
end

function tf=legacy_mean_model(b)
tf=(isfield(b,'schema_version')&&isnumeric(b.schema_version)&&isscalar(b.schema_version)&& ...
    ismember(b.schema_version,1:4))|| ...
    (isfield(b,'model')&&strcmp(b.model,'joint_dsc_p0'))|| ...
    (isfield(b,'parameter_estimate')&&isstruct(b.parameter_estimate)&& ...
    (isfield(b.parameter_estimate,'V')||isfield(b.parameter_estimate,'B')));
end

function tf=is_text_scalar(value)
tf=(ischar(value)&&isrow(value)&&~isempty(value))|| ...
    (isstring(value)&&isscalar(value)&&~ismissing(value)&&strlength(value)>0);
end

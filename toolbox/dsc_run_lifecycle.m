function manifest=dsc_run_lifecycle(action,root,prefix,manifest)
%DSC_RUN_LIFECYCLE Track ownership privately and publish the public manifest.
root=char(root); prefix=char(prefix); action=char(action);
assert(isfolder(root)&&~isempty(regexp(prefix,'^constant_mu_\d{8}T\d{9}_$','once')), ...
    'dsc:RunIdentity','An existing output directory and production prefix are required.');
path=fullfile(root,[prefix 'run_manifest.json']);
state_path=fullfile(root,[prefix 'run_state.json']);
previous=struct();
if isfile(path), previous=jsondecode(fileread(path)); end
if nargin<4, manifest=previous; end
if strcmp(action,'start')
    owner=process_owner();
elseif isfile(state_path)
    owner=jsondecode(fileread(state_path));
elseif isfield(previous,'lifecycle')
    owner=previous.lifecycle; % Compatibility with older unfinished runs.
elseif isfield(manifest,'lifecycle')
    owner=manifest.lifecycle;
elseif strcmp(action,'stopped')
    return % Normal completion already removed the private ownership record.
else
    owner=process_owner();
end
remove_manifest_fields=intersect({'model','mean_model','source_file','tickers','lifecycle'},fieldnames(manifest));
if ~isempty(remove_manifest_fields), manifest=rmfield(manifest,remove_manifest_fields); end
if isfield(manifest,'return_ar1')
    remove_fields=intersect({'calendar_weeks','hac_lag'},fieldnames(manifest.return_ar1));
    if ~isempty(remove_fields), manifest.return_ar1=rmfield(manifest.return_ar1,remove_fields); end
end
if isfield(manifest,'config')
    manifest.config=public_options(manifest.config);
end
defaults=weekly_config();
prior_values=struct('innovation_prior_mean_multiplier',defaults.innovation_prior_mean_multiplier, ...
    'ig_shape',defaults.ig_shape,'mu_prior_std_multiplier',defaults.mu_prior_std_multiplier);
if isfield(manifest,'priors')
    prior_fields=fieldnames(prior_values);
    for k=1:numel(prior_fields)
        field=prior_fields{k};
        if isfield(manifest.priors,field), prior_values.(field)=manifest.priors.(field); end
    end
    % Older saved priors encoded this multiplier only in their covariance.
    p=manifest.priors;
    if ~isfield(p,'mu_prior_std_multiplier')&& ...
            all(isfield(p,{'mu_cov0','mh0','initial_complete_observations'}))
        factors=sqrt(p.initial_complete_observations*diag(p.mu_cov0)./exp(p.mh0(:)));
        assert(all(isfinite(factors))&&all(factors>0)&&max(abs(factors-factors(1)))<1e-10, ...
            'dsc:MeanPriorMultiplier','The saved mu covariance has no common SD multiplier.');
        value=factors(1);
        if abs(value-round(value))<1e-12, value=round(value); end
        prior_values.mu_prior_std_multiplier=value;
    end
end
manifest.prior_parameter_notes=struct( ...
    'innovation_prior_mean_multiplier',struct('value',prior_values.innovation_prior_mean_multiplier, ...
        'description','Multiplies the sample variance of weekly changes in each rolling h/r path to set its inverse-gamma prior mean; calibrated separately for every coordinate.'), ...
    'ig_shape',struct('value',prior_values.ig_shape, ...
        'description','Shape of each inverse-gamma prior for sig2h and sig2r; for shape a > 2, prior SD = prior mean / sqrt(a - 2).'), ...
    'mu_prior_std_multiplier',struct('value',prior_values.mu_prior_std_multiplier, ...
        'description','Normal mu prior SD = this multiplier * sqrt(S0_jj / n0), using the initial-calibration covariance S0 and complete-observation count n0.'));
% Register dependencies before the public manifest exists, so retention can
% protect parameter inputs during the initial data/calibration setup.
if isfield(manifest,'config')
    if isfield(manifest.config,'EstimateParameters')
        owner.estimate_parameters=manifest.config.EstimateParameters;
    end
end
if isfield(manifest,'inference')&&isfield(manifest.inference,'parameter_file')
    owner.parameter_file=manifest.inference.parameter_file;
end
switch action
    case {'start','publishing'}
        owner.status='active';
    case 'finish'
        owner.status='incomplete';
        if isfield(manifest,'run_summary')&&isfield(manifest.run_summary,'status')
            switch char(manifest.run_summary.status)
                case 'iteration_limit', owner.status='completed';
                case 'time_limit', owner.status='stopped';
                case 'failed', owner.status='failed';
            end
        elseif isfile(fullfile(root,[prefix 'parameters.mat'])) && ...
                ~isempty(dir(fullfile(root,[prefix '*posterior_chunk_*.mat'])))
            owner.status='completed';
        end
        owner.finished_at_utc=utc_now();
    case 'failed'
        owner.status='failed';
        if isfield(manifest,'run_summary')&&isfield(manifest.run_summary,'status')&& ...
                strcmp(manifest.run_summary.status,'time_limit')
            owner.status='stopped';
        end
        summary_path=fullfile(root,[prefix 'pilot_summary.json']);
        if isfile(summary_path)
            summary=jsondecode(fileread(summary_path));
            if isfield(summary,'status')&&strcmp(summary.status,'time_limit')
                owner.status='stopped';
            end
        end
        owner.finished_at_utc=utc_now();
    case 'stopped'
        % onCleanup also runs on normal returns. Only an unfinished invocation
        % owned by this process is a manual stop, never a completed publication.
        if ~isfield(owner,'status')||~strcmp(owner.status,'active')|| ...
                ~isfield(owner,'process_id')||owner.process_id~=feature('getpid')
            return
        end
        owner.status='stopped'; owner.finished_at_utc=utc_now();
    otherwise
        error('dsc:RunLifecycle','Unknown lifecycle action %s.',action);
end
owner.updated_at_utc=utc_now();
if strcmp(owner.status,'completed')
    % Keep the active lease until the final manifest has committed. Failure
    % to remove private tracking must not downgrade a completed publication.
    atomic_json(path,manifest,root,prefix);
    try
        atomic_json(state_path,owner,root,prefix);
        delete(state_path);
    catch problem
        warning('dsc:RunStateCleanup','Run published; private state cleanup failed: %s',problem.message);
    end
else
    % Register ownership before the first public snapshot is visible.
    atomic_json(state_path,owner,root,prefix);
    atomic_json(path,manifest,root,prefix);
end
end

function atomic_json(path,value,root,prefix)
temporary=fullfile(root,[prefix 'publication_' char(java.util.UUID.randomUUID()) '.json']);
guard=onCleanup(@()remove_temporary(temporary));
fid=fopen(temporary,'w'); assert(fid>=0,'dsc:RunManifest','Cannot write the run manifest.');
file_guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true)); clear file_guard
[ok,message]=movefile(temporary,path,'f');
assert(ok,'dsc:RunManifest','Cannot publish the run manifest: %s',message);
clear guard
end

function owner=process_owner()
owner=struct('process_id',feature('getpid'), ...
    'host',char(java.net.InetAddress.getLocalHost.getHostName), ...
    'process_started','', 'started_at_utc',utc_now());
if isunix
    [status,value]=system(sprintf('ps -p %.0f -o lstart= 2>&1',owner.process_id));
    if status==0, owner.process_started=strtrim(value); end
end
end

function options=public_options(settings)
% A manifest records callable inputs, never the private sampler configuration.
% Also normalize older in-flight callers when they publish their final files.
names={'WarmupIterations','NumChains','PerDraw','Thin','MaxHours','ChunkSize', ...
    'Seed','EstimateParameters','ParallelChains','CorrelationBackend', ...
    'CorrelationThreads','MeanBackend','CreatePDF','OutputDirectory', ...
    'SourceFile','VariableNames'};
fields={'burnin','num_chains','','thin','','chunk_size','seed','estimate_parameters', ...
    'parallel_chains','correlation_backend','correlation_threads', ...
    'mean_backend','create_pdf','output_root','source_file', ...
    'variable_names'};
options=struct();
for k=1:numel(names)
    name=names{k};
    if isfield(settings,name)
        options.(name)=settings.(name);
    elseif ~isempty(fields{k})&&isfield(settings,fields{k})
        options.(name)=settings.(fields{k});
    elseif strcmp(name,'PerDraw')&&all(isfield(settings,{'max_iterations','burnin','thin'}))
        options.PerDraw=(settings.max_iterations-settings.burnin)/settings.thin;
    elseif strcmp(name,'MaxHours')&&isfield(settings,'max_seconds')
        options.MaxHours=settings.max_seconds/3600;
    end
end
end

function value=utc_now()
value=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end

function remove_temporary(path)
if isfile(path), delete(path); end
end

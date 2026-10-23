function result = run_weekly_research(mode,cfg)
%RUN_WEEKLY_RESEARCH Analyze or pilot the compact constant-mean weekly model.
% 'analysis' never starts MCMC. 'pilot' is capped at 20 warm-up iterations
% and 30 minutes. All outputs use existing directories; pilots use prefixes.
if nargin<1 || isempty(mode), mode='analysis'; end
if nargin<2 || isempty(cfg), cfg=weekly_config(); end
mode=validatestring(mode,{'analysis','pilot'});
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'toolbox'));
assert(isfolder(cfg.output_root),'dsc:OutputDirectory', ...
    'The output directory must already exist: %s',cfg.output_root);
if ~isfield(cfg,'mean_model'), cfg.mean_model='constant'; end
assert(strcmp(cfg.mean_model,'constant'),'dsc:MeanModel', ...
    'Current research runs use a constant mean; historical B/V runs are preserved separately.');
if isfield(cfg,'data_output_dir') && strlength(string(cfg.data_output_dir)) > 0
    outdir = char(cfg.data_output_dir);
else
    outdir = cfg.output_root;
end
assert(isfolder(outdir),'dsc:OutputDirectory', ...
    'The data output directory must already exist: %s',outdir);
% The analysis helper receives an existing directory, so it cannot take its
% historical branch that creates output_root/data.
cfg.data_output_dir=outdir;
if strcmp(mode,'pilot')
    cfg.flat_outputs=true;
    if ~isfield(cfg,'file_prefix')||isempty(cfg.file_prefix)
        assert(~isfield(cfg,'resume')||~cfg.resume,'dsc:ResumePrefix', ...
            'Resume requires the original pilot filename prefix.');
        stamp=char(datetime('now','TimeZone','Europe/Madrid','Format','yyyyMMdd''T''HHmmssSSS'));
        cfg.file_prefix=['pilot_constant_mu_' stamp '_'];
    end
    assert(~isempty(regexp(char(cfg.file_prefix),'^[A-Za-z0-9_-]+$','once')), ...
        'dsc:OutputPrefix','Pilot filename prefixes may contain only letters, numbers, underscores and hyphens.');
end
panel=prepare_weekly_panel(cfg);
[priors,inference]=dsc_prepare_inference(panel,cfg);
if strcmp(mode,'pilot')
    % This bounded pilot invokes a single sampler, not the chain coordinator.
    % Internal continuation uses the already resolved numeric allocation.
    if isfield(cfg,'resume')&&cfg.resume
        requested=cfg.correlation_threads;
        if (ischar(requested)||isstring(requested))&&isequal(string(requested),"auto")
            checkpoint_path=fullfile(cfg.output_root,[cfg.file_prefix 'checkpoint.mat']);
            if ~isfile(checkpoint_path), error('dsc:ResumeMissing','Checkpoint does not exist.'); end
            loaded=load(checkpoint_path,'checkpoint');
            requested=loaded.checkpoint.identity.cfg.correlation_threads;
        end
        cfg.correlation_threads=dsc_resolve_correlation_threads(requested,1,false);
    else
        cfg.correlation_threads=dsc_resolve_correlation_threads( ...
            cfg.correlation_threads,1,false,inference.training_panel,priors,cfg.correlation_backend);
    end
end
summary=analyze_weekly_panel(panel,cfg);
save(fullfile(outdir,'calibrated_priors.mat'),'priors','cfg','-v7.3');
write_json(fullfile(outdir,'prior_summary.json'),priors);
% Manifest schema version, independent of the saved-parameter schema.
manifest.schema_version=1;
manifest.model='joint_dsc_constant_mean';
manifest.mean_model='constant';
manifest.innovation_variances='delta = sig2h, gamma = sig2r; both estimated';
manifest.source_file=cfg.source_file;
manifest.source_hash=panel.source_hash;
manifest.raw_start=char(string(panel.raw_dates(1),'yyyy-MM-dd'));
manifest.raw_end=char(string(panel.raw_dates(end),'yyyy-MM-dd'));
manifest.return_start=char(string(panel.dates(1),'yyyy-MM-dd'));
manifest.return_end=char(string(panel.dates(end),'yyyy-MM-dd'));
manifest.n_returns=size(panel.returns,1);
manifest.n_variables=size(panel.returns,2);
manifest.n_pairs=numel(panel.pair_i);
manifest.tickers=panel.tickers;
manifest.config=cfg;
if isfield(manifest.config,'source_file')
    manifest.config=rmfield(manifest.config,'source_file');
end
manifest.matlab_version=version;
[status,revision]=system(['git -C "' root '" rev-parse HEAD']);
if status==0, manifest.git_revision=strtrim(revision); end
[status,changes]=system(['git -C "' root '" status --porcelain']);
if status==0, manifest.git_dirty=~isempty(strtrim(changes)); end
manifest.matlab_source_hash=code_hash(root);
manifest.created_at=char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
manifest.inference=struct('estimate_parameters',inference.estimate_parameters, ...
    'parameter_estimation_end',inference.estimation_end, ...
    'calibration_weeks',inference.calibration_weeks, ...
    'model_start',char(string(inference.smoothing_panel.dates(1),'yyyy-MM-dd')), ...
    'state_scope','historical smoothing after the reserved initial calibration weeks');
manifest.production_authorized=false;
result=struct('panel',panel,'summary',summary,'priors',priors,'manifest',manifest,'inference',inference);
if strcmp(mode,'pilot')
    if ~inference.estimate_parameters||inference.training_panel.dates(end)~=panel.dates(end)
        error('dsc:PilotMode','Use run_weekly_model_report for saved-parameter smoothing or cutoff estimation.');
    end
    run_dir=cfg.output_root;
    write_json(fullfile(run_dir,[cfg.file_prefix 'run_manifest.json']),manifest);
    result.pilot=run_weekly_pilot(inference.training_panel,priors,cfg,run_dir);
    write_json(fullfile(run_dir,[cfg.file_prefix 'prior_summary.json']),priors);
end
end

function hash=code_hash(root)
% Include tracked and newly authored MATLAB code; Git HEAD alone is not
% sufficient provenance for an uncommitted research implementation.
md=java.security.MessageDigest.getInstance('SHA-256');
folders={'','toolbox'};
for k=1:numel(folders)
    files=dir(fullfile(root,folders{k},'*.m'));
    [~,order]=sort({files.name}); files=files(order);
    for j=1:numel(files)
        relative=fullfile(folders{k},files(j).name);
        md.update(uint8(unicode2native(relative,'UTF-8')));
        md.update(uint8(0));
        fid=fopen(fullfile(root,relative),'rb');
        assert(fid>=0,'DSC:ReadFailed','Cannot hash %s',relative);
        bytes=fread(fid,Inf,'*uint8'); fclose(fid);
        md.update(bytes); md.update(uint8(0));
    end
end
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function write_json(path,value)
fid=fopen(path,'w');
assert(fid>=0,'DSC:WriteFailed','Cannot open %s',path);
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

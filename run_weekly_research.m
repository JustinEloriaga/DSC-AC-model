function result = run_weekly_research(mode,cfg)
%RUN_WEEKLY_RESEARCH Prepare/analyze, and optionally pilot the weekly model.
% 'analysis' never starts MCMC. 'pilot' uses the bounded pilot defaults.
if nargin<1 || isempty(mode), mode='analysis'; end
if nargin<2 || isempty(cfg), cfg=weekly_config(); end
mode=validatestring(mode,{'analysis','pilot'});
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'toolbox'),fullfile(root,'core'));
if ~exist(cfg.output_root,'dir'), mkdir(cfg.output_root); end
panel=prepare_weekly_panel(cfg);
summary=analyze_weekly_panel(panel,cfg);
priors=dsc_calibrate_priors(panel,cfg);
save(fullfile(cfg.output_root,'data','calibrated_priors.mat'),'priors','cfg','-v7.3');
write_json(fullfile(cfg.output_root,'data','prior_summary.json'),priors);
manifest.schema_version=1;
manifest.source_file=cfg.source_file;
manifest.source_hash=panel.source_hash;
manifest.return_start=char(string(panel.dates(1),'yyyy-MM-dd'));
manifest.return_end=char(string(panel.dates(end),'yyyy-MM-dd'));
manifest.n_returns=size(panel.returns,1);
manifest.n_variables=size(panel.returns,2);
manifest.n_pairs=numel(panel.pair_i);
manifest.tickers=panel.tickers;
manifest.config=cfg;
manifest.matlab_version=version;
[status,revision]=system(['git -C "' root '" rev-parse HEAD']);
if status==0, manifest.git_revision=strtrim(revision); end
[status,changes]=system(['git -C "' root '" status --porcelain']);
if status==0, manifest.git_dirty=~isempty(strtrim(changes)); end
manifest.matlab_source_hash=code_hash(root);
manifest.created_at=char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
manifest.inference='historical full-sample smoothing; empirical-Bayes priors';
manifest.production_authorized=false;
write_json(fullfile(cfg.output_root,'run_manifest.json'),manifest);
result=struct('panel',panel,'summary',summary,'priors',priors,'manifest',manifest);
if strcmp(mode,'pilot')
    stamp=char(datetime('now','Format','yyyyMMdd-HHmmss-SSS'));
    run_dir=fullfile(cfg.output_root,'runs',[stamp '-chain' num2str(cfg.chain_id)]);
    if ~exist(run_dir,'dir'), mkdir(run_dir); end
    write_json(fullfile(run_dir,'run_manifest.json'),manifest);
    result.pilot=run_weekly_pilot(panel,priors,cfg,run_dir);
    write_json(fullfile(cfg.output_root,'latest_pilot.json'),result.pilot);
end
end

function hash=code_hash(root)
% Include tracked and newly authored MATLAB code; Git HEAD alone is not
% sufficient provenance for an uncommitted research implementation.
md=java.security.MessageDigest.getInstance('SHA-256');
folders={'','core','toolbox'};
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

function result = run_weekly_sampler_100x10()
%RUN_WEEKLY_SAMPLER_100X10 User-requested short posterior pilot.
% Runs 100 completed sweeps with 10 warm-up sweeps, retaining post-warm-up
% packed P(i,j,t) draws. This is not a convergence assessment.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root,'toolbox'), fullfile(root,'core'));
cfg = weekly_config();
cfg.max_iterations = 100;
cfg.burnin = 10;
cfg.max_seconds = 7200;
cfg.thin = 1;
cfg.chunk_size = 10;
cfg.resume = false;
cfg.chain_id = 1;
cfg.seed = 20260914;
cfg.correlation_backend = 'auto';
cfg.correlation_threads = 4;

if ~exist(cfg.output_root,'dir'), mkdir(cfg.output_root); end
panel = prepare_weekly_panel(cfg);
summary = analyze_weekly_panel(panel,cfg);
priors = dsc_calibrate_priors(panel,cfg);
data_dir = fullfile(cfg.output_root,'data');
if ~exist(data_dir,'dir'), mkdir(data_dir); end
save(fullfile(data_dir,'calibrated_priors.mat'),'priors','cfg','-v7.3');
write_json_local(fullfile(data_dir,'prior_summary.json'),priors);

stamp = char(datetime('now','Format','yyyyMMdd-HHmmss-SSS'));
run_dir = fullfile(cfg.output_root,'runs',[stamp '-chain' num2str(cfg.chain_id) '-100x10']);
mkdir(run_dir);
run_manifest = struct();
run_manifest.schema_version = 1;
run_manifest.description = 'Short user-requested Bayesian sampler run: 100 sweeps, 10 warm-up; not convergence evidence.';
run_manifest.source_file = cfg.source_file;
run_manifest.source_hash = panel.source_hash;
run_manifest.return_start = char(string(panel.dates(1),'yyyy-MM-dd'));
run_manifest.return_end = char(string(panel.dates(end),'yyyy-MM-dd'));
run_manifest.n_returns = size(panel.returns,1);
run_manifest.n_variables = size(panel.returns,2);
run_manifest.n_pairs = numel(panel.pair_i);
run_manifest.tickers = panel.tickers;
run_manifest.config = cfg;
run_manifest.panel_summary = summary;
run_manifest.created_at = char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
run_manifest.convergence_established = false;
write_json_local(fullfile(run_dir,'run_manifest.json'),run_manifest);

result = dsc_sample(panel,priors,cfg,run_dir);
write_json_local(fullfile(cfg.output_root,'latest_100x10.json'),result);
end

function write_json_local(path,value)
fid = fopen(path,'w');
assert(fid>=0,'DSC:WriteFailed','Cannot open %s',path);
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

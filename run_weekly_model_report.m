function result = run_weekly_model_report(options)
%RUN_WEEKLY_MODEL_REPORT Run the weekly DSC model and build its PDF report.
% Example:
%   result = run_weekly_model_report(WarmupIterations=10,RetainedDraws=100);
%
% RetainedDraws is the number of draws per chain after warm-up and thinning.
% The sampler therefore performs WarmupIterations + RetainedDraws*Thin
% completed sweeps unless MaxHours stops it first.
arguments
    options.WarmupIterations (1,1) double {mustBeInteger,mustBeNonnegative} = 10
    options.RetainedDraws (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.RetainedDraws,2)} = 100
    options.Thin (1,1) double {mustBeInteger,mustBePositive} = 1
    options.MaxHours (1,1) double {mustBePositive,mustBeFinite} = 12
    options.ChunkSize (1,1) double {mustBeInteger,mustBePositive} = 1000
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260914
    options.ChainID (1,1) double {mustBeInteger,mustBePositive} = 1
    options.NumChains (1,1) double {mustBeInteger,mustBePositive} = 1
    options.ParallelChains (1,1) logical = false
    options.ConvergenceMode (1,1) string {mustBeMember(options.ConvergenceMode,["off","report"])} = "report"
    options.AutoExtend (1,1) logical = false
    options.MaxRetainedDraws (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.MaxRetainedDraws,2)} = 2000
    options.RhatThreshold (1,1) double {mustBeGreaterThan(options.RhatThreshold,1),mustBeFinite} = 1.01
    options.MinESS (1,1) double {mustBePositive,mustBeFinite} = 400
    options.MaxMCSERatio (1,1) double {mustBePositive,mustBeFinite} = 0.05
    options.CorrelationBackend (1,1) string {mustBeMember(options.CorrelationBackend,["auto","mex","matlab"])} = "mex"
    options.CorrelationThreads = "auto"
    options.SourceFile (1,1) string = ""
    options.VariableNames (1,:) string = "all"
    options.EstimateParameters (1,1) logical = true
    options.EstimationEndDate (1,1) string = ""
    options.KeepNewestParameterFiles (1,1) double {mustBeInteger,mustBeNonnegative} = 3
    options.ParameterSmoothing (1,1) string {mustBeMember(options.ParameterSmoothing,["draws","mean"])} = "mean"
    options.RunTests (1,1) logical = true
    options.BuildReport (1,1) logical = false
    options.VerifyReport (1,1) logical = false
    options.PythonExecutable (1,1) string = "python3"
end

if options.VerifyReport && ~options.BuildReport
    error('DSC:PipelineOptions','VerifyReport requires BuildReport=true.');
end
if isstring(options.CorrelationThreads) || ischar(options.CorrelationThreads)
    thread_text=string(options.CorrelationThreads);
    if thread_text=="auto"
        if ~options.RunTests
            warning('DSC:AutoThreadsNeedsTests', ...
                'CorrelationThreads="auto" requires RunTests=true; enabling RunTests automatically.');
            options.RunTests=true;
        end
        options.CorrelationThreads=min(4,feature('numcores'));
    else
        options.CorrelationThreads=str2double(thread_text);
    end
end
if ~isscalar(options.CorrelationThreads)||~isfinite(options.CorrelationThreads)|| ...
        options.CorrelationThreads<1||options.CorrelationThreads~=floor(options.CorrelationThreads)|| ...
        options.CorrelationThreads>64
    error('DSC:CorrelationThreads','CorrelationThreads must be "auto" or an integer from 1 to 64.');
end
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'toolbox'),fullfile(root,'scripts'));

cfg = weekly_config();
if strlength(options.SourceFile)>0, cfg.source_file=char(options.SourceFile); end
cfg.variable_names = options.VariableNames;
if ~exist(cfg.source_file,'file')
    error('DSC:SourceFile','Current-data source file does not exist: %s',cfg.source_file);
end
cfg.burnin = options.WarmupIterations;
cfg.thin = options.Thin;
cfg.max_iterations = options.WarmupIterations + options.RetainedDraws*options.Thin;
cfg.max_seconds = options.MaxHours*3600;
cfg.chunk_size = options.ChunkSize;
cfg.seed = options.Seed;
cfg.chain_id = options.ChainID;
cfg.num_chains = options.NumChains;
cfg.parallel_chains = options.ParallelChains;
cfg.convergence_mode = char(options.ConvergenceMode);
cfg.auto_extend = options.AutoExtend;
cfg.check_every = cfg.chunk_size;
cfg.max_retained_draws = options.MaxRetainedDraws;
cfg.min_diagnostic_draws = cfg.chunk_size;
cfg.rhat_threshold = options.RhatThreshold;
cfg.min_ess = options.MinESS;
cfg.max_mcse_ratio = options.MaxMCSERatio;
cfg.resume = false;
cfg.correlation_backend = char(options.CorrelationBackend);
cfg.correlation_threads = options.CorrelationThreads;
cfg.estimate_parameters = options.EstimateParameters;
cfg.estimation_end_date = char(options.EstimationEndDate);
cfg.parameter_dir = fullfile(cfg.output_root,'parameters');
cfg.parameter_file = '';
cfg.parameter_smoothing = char(options.ParameterSmoothing);
cfg.write_publication_exports = options.BuildReport;
cfg.keep_newest_parameter_files = options.KeepNewestParameterFiles;
parameter_cleanup = onCleanup(@()delete_old_parameter_files( ...
    cfg.parameter_dir,cfg.keep_newest_parameter_files)); %#ok<NASGU>
cfg.data_output_dir = dsc_create_timestamped_directory( ...
    fullfile(cfg.output_root,'data'),'run-',['-chain' num2str(cfg.chain_id)]);
data_workspace_cleanup = onCleanup(@()delete_data_workspace(cfg.data_output_dir)); %#ok<NASGU>

% Reject unsupported execution settings before data preparation or sampling.
if cfg.auto_extend&&(strcmp(cfg.convergence_mode,'off')||cfg.max_retained_draws<options.RetainedDraws)
    error('dsc:ConvergenceConfig','AutoExtend needs convergence checking and MaxRetainedDraws >= RetainedDraws.');
end
if cfg.parallel_chains&&~dsc_parallel_available()
    error('dsc:ParallelUnavailable','ParallelChains=true requires Parallel Computing Toolbox and an available license. Use ParallelChains=false for sequential chains.');
end

if options.RunTests
    run_weekly_acceptance(cfg);
end
analysis = run_weekly_research('analysis',cfg);

run_name = sprintf('chain%d-w%d-r%d-thin%d',cfg.chain_id,cfg.burnin,options.RetainedDraws,cfg.thin);
run_dir = dsc_create_timestamped_directory(fullfile(cfg.output_root,'runs'),'', ...
    ['-' run_name]);
[~,run_name]=fileparts(run_dir);
run_marker=fullfile(run_dir,'.run_in_progress');
marker_fid=fopen(run_marker,'w');
if marker_fid<0, error('DSC:RunMarker','Cannot mark active run: %s',run_dir); end
fclose(marker_fid);
run_marker_cleanup=onCleanup(@()delete_file_if_present(run_marker)); %#ok<NASGU>
intermediate_cleanup = onCleanup(@()delete_intermediate_run_files(run_dir)); %#ok<NASGU>
run_manifest = analysis.manifest;
run_manifest.run_kind = 'configurable sampler and report pipeline';
run_manifest.requested_warmup_iterations = options.WarmupIterations;
run_manifest.requested_retained_draws = options.RetainedDraws; % legacy alias: per chain
run_manifest.requested_retained_draws_per_chain = options.RetainedDraws;
run_manifest.requested_total_sweeps = cfg.max_iterations;
run_manifest.requested_num_chains = cfg.num_chains;
run_manifest.parallel_chains = cfg.parallel_chains;
run_manifest.convergence_mode = cfg.convergence_mode;
run_manifest.requested_max_hours = options.MaxHours;
run_manifest.convergence_established = false;
write_json_local(fullfile(run_dir,'run_manifest.json'),run_manifest);

result = dsc_run_inference(analysis.panel,analysis.priors,cfg,run_dir,analysis.inference);
if result.saved_draws < 2
    error('DSC:TooFewRetainedDraws', ...
        ['The run was preserved at %s, but it retained only %d draw(s). ' ...
         'Increase MaxHours or reduce WarmupIterations.'],run_dir,result.saved_draws);
end
if any(result.saved_draws_per_chain < options.RetainedDraws)
    warning('DSC:IncompleteDrawTarget', ...
        'Some chains retained fewer than %d requested draws each; %d draws were retained in total.', ...
        options.RetainedDraws,result.saved_draws);
end

result.figure_paths = plot_bayes_correlation_paths(run_dir);
result.convergence_artifacts = write_convergence_artifacts(run_dir,result);
result.report_pdf = '';
result.figure_pdf = '';
result.publication_checked = false;
if options.BuildReport
    publication_exports_cleanup = onCleanup(@()delete_publication_exports(cfg.data_output_dir)); %#ok<NASGU>
    builder = fullfile(root,'scripts','build_publication.py');
    summary_path = fullfile(run_dir,'pilot_summary.json');
    command = sprintf('%s %s --data %s --pilot %s --posterior-run %s --run-id %s', ...
        shell_arg(options.PythonExecutable),shell_arg(builder),shell_arg(cfg.data_output_dir), ...
        shell_arg(summary_path),shell_arg(run_dir),shell_arg(run_name));
    [status,output] = system(command,'-echo');
    if status~=0
        error('DSC:PublicationBuild','Report builder failed with status %d:\n%s',status,output);
    end
    result.report_pdf = fullfile(root,'report',['report_' run_name '.pdf']);
    result.figure_pdf = fullfile(root,'report',['figure_' run_name '.pdf']);
    if options.VerifyReport
        checker = fullfile(root,'scripts','check_publication.py');
        command = sprintf('%s %s --run-id %s',shell_arg(options.PythonExecutable), ...
            shell_arg(checker),shell_arg(run_name));
        [status,output] = system(command,'-echo');
        if status~=0
            error('DSC:PublicationCheck','Report verification failed with status %d:\n%s',status,output);
        end
        result.publication_checked = true;
    end
end

latest = struct('run_dir',run_dir,'summary',fullfile(cfg.output_root,'latest_model_report.json'), ...
    'figures',fullfile(run_dir,'figures'),'report_pdf',result.report_pdf,'figure_pdf',result.figure_pdf, ...
    'requested_retained_draws',options.RetainedDraws, ...
    'requested_retained_draws_per_chain',options.RetainedDraws, ...
    'actual_retained_draws',result.saved_draws, ...
    'actual_retained_draws_total',result.saved_draws, ...
    'warmup_iterations',cfg.burnin,'thin',cfg.thin,'completed_iterations',result.completed_iterations, ...
    'status',result.status,'created_at',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')));
latest.inference=result.inference;
latest.convergence_artifacts=result.convergence_artifacts;
latest.manifest=run_manifest;
latest.num_chains=cfg.num_chains;
latest.actual_retained_draws_per_chain=result.saved_draws_per_chain;
latest.convergence=result.convergence;
latest.run_result=struct('status',result.status,'reason',result.reason, ...
    'completed_iterations',result.completed_iterations, ...
    'completed_iterations_per_chain',result.completed_iterations_per_chain, ...
    'saved_draws_total',result.saved_draws, ...
    'saved_draws_per_chain',result.saved_draws_per_chain, ...
    'chain_ids',result.chain_ids,'seeds',result.seeds, ...
    'total_seconds',result.total_seconds,'stage_seconds',result.stage_seconds);
write_json_local(fullfile(cfg.output_root,'latest_model_report.json'),latest);
delete_file_if_present(run_marker);
dsc_prune_outputs(cfg.output_root,fullfile(root,'report'),run_dir);
end

function delete_file_if_present(path)
if isfile(path), delete(path); end
end

function delete_publication_exports(data_dir)
files = {'weekly_levels.csv','weekly_returns.csv','weekly_panel.mat', ...
    'weekly_quality.csv','data_summary.json'};
for k=1:numel(files)
    path = fullfile(data_dir,files{k});
    if isfile(path), delete(path); end
end
end

function delete_data_workspace(data_dir)
if isfolder(data_dir), rmdir(data_dir,'s'); end
data_parent=fileparts(data_dir);
if isfolder(data_parent)
    entries=dir(data_parent);
    entries=entries(~ismember({entries.name},{'.','..'}));
    if isempty(entries)
        [removed,message]=rmdir(data_parent);
        if ~removed
            remaining=dir(data_parent);
            remaining=remaining(~ismember({remaining.name},{'.','..'}));
            if isempty(remaining)
                warning('DSC:DataWorkspaceCleanup', ...
                    'Could not remove the empty temporary data folder %s: %s',data_parent,message);
            end
        end
    end
end
end

function delete_old_parameter_files(parameter_dir,keep_count)
if keep_count <= 0 || ~isfolder(parameter_dir), return; end
files = dir(fullfile(char(parameter_dir),'parameters_*.mat'));
if numel(files) <= keep_count, return; end
names=string({files.name});
[~,order] = sort(names,'descend');
for k=keep_count+1:numel(order)
    delete(fullfile(files(order(k)).folder,files(order(k)).name));
end
end

function delete_intermediate_run_files(run_dir)
% Remove files that can be regenerated from retained posterior chunks.
patterns = {'posterior_correlation_bands.mat','checkpoint.mat','diagnostics.mat','result.mat','*.DS_Store'};
for p=1:numel(patterns)
    files=dir(fullfile(run_dir,'**',patterns{p}));
    for k=1:numel(files)
        if ~files(k).isdir, delete(fullfile(files(k).folder,files(k).name)); end
    end
end
% These JSON files are runtime handoffs and are reproducible from the
% retained chunks, result, manifest, and convergence diagnostics.
for name={'pilot_summary.json','inference_metadata.json','prior_summary.json', ...
        'bayes_correlation_paths_manifest.json'}
    files=dir(fullfile(run_dir,'**',name{1}));
    for k=1:numel(files)
        if ~files(k).isdir, delete(fullfile(files(k).folder,files(k).name)); end
    end
end
chain_dirs=dir(fullfile(run_dir,'chains','chain_*'));
for c=1:numel(chain_dirs)
    if ~chain_dirs(c).isdir, continue; end
    folder=fullfile(chain_dirs(c).folder,chain_dirs(c).name);
    for name={'diagnostics.mat','result.mat','pilot_summary.json','checkpoint.mat'}
        path=fullfile(folder,name{1});
        if isfile(path), delete(path); end
    end
end
root_result=fullfile(run_dir,'result.mat');
if isfile(root_result), delete(root_result); end
end

function value = shell_arg(value)
value = char(value);
quote = char(39);
value = strrep(value,quote,[quote '"' quote '"' quote]);
value = [quote value quote];
end

function write_json_local(path,value)
fid=fopen(path,'w');
assert(fid>=0,'DSC:WriteFailed','Cannot write %s.',path);
guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

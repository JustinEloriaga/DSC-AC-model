function result = dsc_run_inference(panel,priors,cfg,run_dir,inference)
%DSC_RUN_INFERENCE Estimate/save static parameters or condition on saved ones.
% With an earlier cutoff, estimation and full-history conditional smoothing
% are separate stages. Each has the configured draw target and time budget.
estimation_result=[];
modeled=inference.smoothing_panel;
parameter_smoothing=validatestring(inference.parameter_smoothing,{'draws','mean'});
if inference.estimate_parameters
    training=inference.training_panel;
    extends_history=training.dates(end)<modeled.dates(end);
    needs_conditional_smoothing=extends_history||strcmp(parameter_smoothing,'mean');
    if needs_conditional_smoothing
        estimation_dir=fullfile(run_dir,'estimation');
    else
        estimation_dir=run_dir;
    end
    estimation_cfg=cfg; estimation_cfg.estimate_parameters=true;
    fprintf('Estimating parameters through %s (%d weekly returns).\n', ...
        inference.estimation_end,numel(training.dates));
    estimation_result=dsc_sample(training,priors,estimation_cfg,estimation_dir);
    require_retained(estimation_result);
    requested=floor((cfg.max_iterations-cfg.burnin)/cfg.thin);
    if estimation_result.saved_draws<requested
        warning('dsc:IncompleteParameterEstimation', ...
            'Parameter estimation retained %d of %d requested draws. Saving the partial estimate with its actual count.', ...
            estimation_result.saved_draws,requested);
    end
    parameter_path=dsc_save_parameters(estimation_result,training,priors, ...
        estimation_cfg,inference.parameter_dir,inference.calibration_panel);
    parameters=dsc_load_parameters(inference.parameter_dir,panel,parameter_path);
    if ~needs_conditional_smoothing
        result=estimation_result;
    end
else
    parameters=inference.parameters;
    needs_conditional_smoothing=true; % Run conditional smoothing even with no new dates.
end
if needs_conditional_smoothing
    smoothing_cfg=cfg;
    smoothing_cfg.estimate_parameters=false;
    smoothing_cfg.parameter_smoothing=parameter_smoothing;
    if strcmp(parameter_smoothing,'draws')
        smoothing_cfg.fixed_parameter_draws=parameters.parameter_draws;
        smoothing_cfg.fixed_parameters=parameters.parameter_estimate;
    else
        smoothing_cfg.fixed_parameters=parameters.parameter_estimate;
    end
    smoothing_cfg.parameter_file=parameters.path;
    smoothing_priors=parameters.priors;
    smoothing_priors.estimation_source_hash=parameters.source_hash;
    smoothing_priors.source_hash=panel.source_hash;
    fprintf('Smoothing states through %s with fixed parameters estimated at %s.\n', ...
        char(string(modeled.dates(end),'yyyy-MM-dd')),parameters.estimated_at);
    result=dsc_sample(modeled,smoothing_priors,smoothing_cfg,run_dir);
    require_retained(result);
end
metadata=struct('estimate_parameters_requested',inference.estimate_parameters, ...
    'inference_mode',result.inference_mode,'parameter_file',parameters.path, ...
    'parameters_estimated_at',parameters.estimated_at, ...
    'parameter_estimation_start',char(string(parameters.estimation_start,'yyyy-MM-dd')), ...
    'parameter_estimation_end',char(string(parameters.estimation_end,'yyyy-MM-dd')), ...
    'parameter_estimation_draws',parameters.retained_draws, ...
    'parameter_estimation_run',parameters.run_dir, ...
    'calibration_weeks',parameters.calibration_weeks, ...
    'calibration_start',char(string(parameters.calibration_start,'yyyy-MM-dd')), ...
    'calibration_end',char(string(parameters.calibration_end,'yyyy-MM-dd')), ...
    'full_data_start',char(string(panel.dates(1),'yyyy-MM-dd')), ...
    'smoothing_start',char(string(modeled.dates(1),'yyyy-MM-dd')), ...
    'smoothing_end',char(string(panel.dates(end),'yyyy-MM-dd')), ...
    'parameter_smoothing',parameter_smoothing, ...
    'fixed_parameter_summary',parameter_summary_text(parameter_smoothing), ...
    'state_scope','B, h and r jointly smoothed from smoothing_start through smoothing_end; calibration weeks excluded; historical modeled states may revise.', ...
    'parameter_uncertainty_in_bands',strcmp(result.inference_mode,'parameter_estimation')||strcmp(parameter_smoothing,'draws'));
if ~isempty(estimation_result)
    metadata.parameter_estimation_status=estimation_result.status;
elseif isfield(parameters,'estimation_status')
    metadata.parameter_estimation_status=parameters.estimation_status;
end
if isfield(parameters,'requested_retained_draws')
    metadata.parameter_estimation_requested_draws=parameters.requested_retained_draws;
end
result.inference=metadata;
result.paths.parameters=parameters.path;
% Snapshot the priors used by this run: another cutoff can overwrite the
% shared analysis directory even when the CSV source hash is unchanged.
priors=result.priors;
write_json(fullfile(run_dir,'prior_summary.json'),priors);
temporary=[tempname(run_dir),'.mat'];
save(temporary,'priors','cfg','-v7.3');
movefile(temporary,fullfile(run_dir,'calibrated_priors.mat'),'f');
write_json(fullfile(run_dir,'inference_metadata.json'),metadata);
% Keep the selected run's manifest and result authoritative for rebuilding.
manifest_path=fullfile(run_dir,'run_manifest.json');
if isfile(manifest_path), manifest=jsondecode(fileread(manifest_path));
else, manifest=struct('source_hash',panel.source_hash); end
manifest.inference=metadata;
write_json(manifest_path,manifest);
temporary=[tempname(run_dir),'.mat'];
save(temporary,'result','-v7.3');
movefile(temporary,fullfile(run_dir,'result.mat'),'f');
end

function text=parameter_summary_text(parameter_smoothing)
if strcmp(parameter_smoothing,'draws')
    text='Retained parameter draws of V, sig2h and sig2r.';
else
    text='Retained-draw posterior means of V, sig2h and sig2r.';
end
end

function require_retained(result)
if strcmp(result.status,'failed')
    error('dsc:InferenceFailed','Sampler failed; completed work is preserved at %s. %s',result.run_dir,result.reason);
end
if result.saved_draws<2
    error('dsc:TooFewParameterDraws', ...
        'Only %d retained draw(s) at %s. Increase the time budget or reduce warm-up.',result.saved_draws,result.run_dir);
end
end

function write_json(path,value)
temporary=[tempname(fileparts(path)),'.json'];
fid=fopen(temporary,'w');
assert(fid>=0,'dsc:OutputFile','Cannot write %s.',path);
guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true)); clear guard
movefile(temporary,path,'f');
end

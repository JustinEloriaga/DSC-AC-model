function [priors,inference] = dsc_prepare_inference(panel,cfg)
%DSC_PREPARE_INFERENCE Calibrate through the cutoff or load frozen parameters.
% No sampler is launched here. Appended observations never recalibrate a
% loaded model. The first prior_weeks are excluded from every likelihood.
if ~isfield(cfg,'estimate_parameters'), cfg.estimate_parameters=true; end
if ~islogical(cfg.estimate_parameters)||~isscalar(cfg.estimate_parameters)
    error('dsc:InferenceConfig','estimate_parameters must be a scalar logical.');
end
if ~isfield(cfg,'estimation_end_date'), cfg.estimation_end_date=''; end
if ~isfield(cfg,'parameter_dir')||strlength(string(cfg.parameter_dir))==0
    cfg.parameter_dir=fullfile(cfg.output_root,'parameters');
end
if ~isfield(cfg,'parameter_file'), cfg.parameter_file=''; end
if ~isfield(cfg,'parameter_smoothing')||strlength(string(cfg.parameter_smoothing))==0
    cfg.parameter_smoothing='mean';
end
cfg.parameter_smoothing=validatestring(cfg.parameter_smoothing,{'draws','mean'});
inference=struct('estimate_parameters',cfg.estimate_parameters, ...
    'parameter_dir',char(cfg.parameter_dir),'parameters',[],'training_panel',[], ...
    'calibration_panel',[],'smoothing_panel',[], ...
    'parameter_smoothing',char(cfg.parameter_smoothing));
if cfg.estimate_parameters
    if ~isfield(cfg,'prior_weeks')||~isscalar(cfg.prior_weeks)|| ...
            ~isfinite(cfg.prior_weeks)||cfg.prior_weeks<1||cfg.prior_weeks~=floor(cfg.prior_weeks)
        error('dsc:InferenceConfig','prior_weeks must be a positive integer.');
    end
    excluded=cfg.prior_weeks;
    if strlength(string(cfg.parameter_file))>0
        error('dsc:InferenceConfig','parameter_file applies only when estimate_parameters=false.');
    end
    last=numel(panel.dates);
    if strlength(string(cfg.estimation_end_date))>0
        cutoff=datetime(cfg.estimation_end_date,'InputFormat','yyyy-MM-dd');
        if ~isscalar(cutoff)||isnat(cutoff)||cutoff<panel.dates(1)||cutoff>panel.dates(end)
            error('dsc:EstimationDate','Estimation cutoff must lie within the available return dates.');
        end
        last=find(panel.dates<=cutoff,1,'last');
    end
    if last<=excluded
        error('dsc:EstimationDate','The estimation cutoff must follow the %d initial calibration weeks.',excluded);
    end
    % Initial moments use only the reserved prefix. Preserve the existing
    % rolling evolution calibration through the cutoff, then remove that
    % initial prefix physically from the sampler's observation arrays.
    priors=dsc_calibrate_priors(dsc_panel_slice(panel,1,last),cfg);
    training=dsc_panel_slice(panel,excluded+1,last);
    priors.likelihood_start_date=training.dates(1);
    priors.excluded_initial_weeks=excluded;
    priors.calibration_dates=panel.dates([1 excluded]);
    priors.scope='Initial prior_weeks excluded from likelihood; evolution priors use rolling windows through the estimation cutoff.';
    inference.training_panel=training;
    inference.estimation_end=char(string(training.dates(end),'yyyy-MM-dd'));
else
    if strlength(string(cfg.estimation_end_date))>0
        error('dsc:InferenceConfig','The saved model determines the estimation cutoff; leave estimation_end_date empty when loading.');
    end
    parameters=dsc_load_parameters(cfg.parameter_dir,panel,cfg.parameter_file);
    excluded=parameters.calibration_weeks;
    priors=parameters.priors;
    % The loader has verified the exact training prefix. Preserve its source
    % identity separately while binding this copy to the extended panel.
    priors.estimation_source_hash=parameters.source_hash;
    priors.source_hash=panel.source_hash;
    inference.parameters=parameters;
    inference.estimation_end=char(string(parameters.estimation_end,'yyyy-MM-dd'));
    fprintf('Loaded parameters: %s\nEstimated at: %s\nEstimated using returns through: %s\n', ...
        parameters.path,parameters.estimated_at,inference.estimation_end);
end
inference.calibration_weeks=excluded;
inference.calibration_panel=dsc_panel_slice(panel,1,excluded);
inference.smoothing_panel=dsc_panel_slice(panel,excluded+1,numel(panel.dates));
fprintf('Initial-prior calibration: %d weeks (%s to %s), excluded from the likelihood.\n', ...
    excluded,char(string(panel.dates(1),'yyyy-MM-dd')),char(string(panel.dates(excluded),'yyyy-MM-dd')));
fprintf('Model sample: %s to %s (%d weeks available for smoothing).\n', ...
    char(string(panel.dates(excluded+1),'yyyy-MM-dd')), ...
    char(string(panel.dates(end),'yyyy-MM-dd')),numel(panel.dates)-excluded);
end

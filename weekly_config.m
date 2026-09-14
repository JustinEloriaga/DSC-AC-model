function cfg = weekly_config()
%WEEKLY_CONFIG Reproducible defaults for the staged weekly research workflow.
root = fileparts(mfilename('fullpath'));
cfg.source_file = fullfile(root,'data','yad_tickers_no_usdcnh_from_20030804.csv');
cfg.start_date = '2003-08-04';
cfg.closure_policy = 'asof';
cfg.output_root = fullfile(root,'outputs','weekly_research');
cfg.freq = 'weekly';
cfg.p = 0;
cfg.seed = 20260911;
cfg.chain_id = 1;
cfg.prior_weeks = 104;
cfg.shrinkage = 0.05;
cfg.calibration_windows = [52 104 156];
cfg.calibration_window = 104;
cfg.ig_shape = 10;
cfg.kB = 0.01;
cfg.max_iterations = 20;
cfg.max_seconds = 1800;
cfg.burnin = 20; % Every draw in the initial pilot is warm-up.
cfg.thin = 1;
cfg.chunk_size = 10;
cfg.resume = false;
cfg.correlation_backend = 'auto'; % Use the validated MEX kernel if built; otherwise MATLAB.
cfg.correlation_threads = 4; % Independent date blocks; no parallel chains in the pilot.
end

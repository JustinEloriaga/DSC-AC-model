function result = run_weekly_pilot(panel,priors,cfg,run_dir)
% Explicit caller-controlled pilot; never claims posterior convergence.
cfg.max_iterations=min(cfg.max_iterations,20);
cfg.max_seconds=min(cfg.max_seconds,1800);
cfg.burnin=max(cfg.burnin,cfg.max_iterations);
result=dsc_sample(panel,priors,cfg,run_dir);
end

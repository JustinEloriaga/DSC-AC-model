% Weekly research entrypoint. Default action does not run MCMC.
% See weekly_config.m and README.md for the explicit bounded pilot command.
%result = run_weekly_research('analysis');

result = run_weekly_model_report(WarmupIterations=10,RetainedDraws=100,MaxHours=25);

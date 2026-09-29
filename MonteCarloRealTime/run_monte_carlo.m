function result = run_monte_carlo(N,repetitions,output_dir)
% Two-series DSC: compare all four dates with posterior(1:3) updated by date 4.
% Paired importance sampling preserves identity without intermediate resampling.
if nargin<1, N=20000; end
if nargin<2, repetitions=5; end
validateattributes(N,{'numeric'},{'scalar','integer','>=',100});
validateattributes(repetitions,{'numeric'},{'scalar','integer','>=',1});
folder=fileparts(mfilename('fullpath'));
if nargin<3
    output_dir=fullfile(folder,'results',['run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
end
if exist(output_dir,'dir')&&~isempty(dir(fullfile(output_dir,'REPORT.md')))
    error('rt:ExistingRun','Choose a new output directory to preserve an earlier exercise.');
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
original_rng=rng; cleanup=onCleanup(@()rng(original_rng));
started=tic;

prior=struct('mB',[0 0],'CB',[.25 .05;.05 .25], ...
    'Psi',[.15 .03;.03 .20],'nu',8,'mh',log([.7 1.2]),'mr',atanh(.2), ...
    'ch',10,'cr',10,'a',6,'bh',.15,'br',.05);
data_seed=20260929;
truth=dsc_rt_draw_prior(1,4,prior,data_seed);
Y=zeros(4,2);
for t=1:4
    rho=tanh(truth.r(t)); sd=exp(truth.h(:,:,t)/2);
    sigma=(sd'*sd).*[1 rho;rho 1];
    Y(t,:)=truth.B(:,:,t)+randn(1,2)*chol(sigma);
end
data_table=array2table([(1:4)' Y],'VariableNames',{'t','y1','y2'});
writetable(data_table,fullfile(output_dir,'observations.csv'));
write_json(fullfile(output_dir,'prior.json'),prior);
save(fullfile(output_dir,'synthetic_data.mat'),'Y','truth','prior','data_seed');

diagnostics=zeros(repetitions,10); parameter_means=zeros(repetitions,6);
for rep=1:repetitions
    seed=6400+rep;
    batch_particles=dsc_rt_draw_prior(N,4,prior,seed);
    batch_loglik=dsc_rt_loglik(batch_particles,Y,1:4);
    [batch_weights,batch_logweights,batch_logtotal]=dsc_rt_normalize(sum(batch_loglik,2));
    [batch_hyper,hyper_names]=dsc_rt_hyperparameters(batch_particles,prior);
    batch=struct('particles',batch_particles,'weights',batch_weights, ...
        'log_weights',batch_logweights,'hyper',batch_hyper, ...
        'log_evidence',batch_logtotal-log(N));

    old_particles=dsc_rt_draw_prior(N,3,prior,seed);
    old_rng=rng;
    old_loglik=dsc_rt_loglik(old_particles,Y(1:3,:),1:3);
    [old_weights,old_logweights,old_logtotal]=dsc_rt_normalize(sum(old_loglik,2));
    old_hyper=dsc_rt_hyperparameters(old_particles,prior);
    previous=struct('particles',old_particles,'weights',old_weights, ...
        'log_weights',old_logweights,'hyper',old_hyper,'prior',prior, ...
        'rng_state',old_rng,'log_evidence',old_logtotal-log(N));
    if rep==1
        save(fullfile(output_dir,'posterior_T3.mat'),'previous','-v7');
        loaded=load(fullfile(output_dir,'posterior_T3.mat'),'previous');
        previous=loaded.previous;
    end
    sequential=dsc_rt_update(previous,Y(4,:));
    assert(isequal(batch_particles,sequential.particles),'rt:ParticlePairing', ...
        'Common random numbers must produce exactly the same candidate paths.');
    hyper_error=max(abs(batch.hyper-sequential.hyper),[],'all');
    weight_error=max(abs(batch.weights-sequential.weights));
    evidence_error=abs(batch.log_evidence-sequential.log_evidence);
    [batch_mean,batch_sd,parameter_names]=parameter_moments(batch.hyper,batch.weights);
    [seq_mean,seq_sd]=parameter_moments(sequential.hyper,sequential.weights);
    mean_error=max(abs(batch_mean-seq_mean));
    sd_error=max(abs(batch_sd-seq_sd));
    assert(max([hyper_error weight_error evidence_error mean_error sd_error])<1e-11, ...
        'rt:Equality','The two posterior representations must agree to floating-point precision.');
    ess3=1/sum(previous.weights.^2); ess4=1/sum(batch.weights.^2);
    diagnostics(rep,:)=[rep,seed,ess3,ess4,hyper_error,weight_error,evidence_error,mean_error,sd_error,max(batch.weights)];
    parameter_means(rep,:)=batch_mean;
    fprintf('Repeat %d/%d: ESS %.0f -> %.0f of %d; hyper diff %.3g, weight diff %.3g\n', ...
        rep,repetitions,ess3,ess4,N,hyper_error,weight_error);
    if rep==1
        first=struct('batch',batch,'sequential',sequential,'previous',previous);
        save(fullfile(output_dir,'posterior_T4.mat'),'batch','sequential','-v7');
        prior_values=[prior.nu prior.Psi(1,1) prior.Psi(1,2) prior.Psi(2,2) prior.a prior.bh prior.bh prior.a prior.br];
        hyper_table=table(hyper_names',prior_values', ...
            (previous.weights'*previous.hyper)',(batch.weights'*batch.hyper)', ...
            (sequential.weights'*sequential.hyper)', ...
            max(abs(batch.hyper-sequential.hyper),[],1)', ...
            'VariableNames',{'Hyperparameter','Prior','After_1_3','Batch_1_4','Sequential_3_then_4','Max_component_difference'});
        [old_mean,old_sd]=parameter_moments(previous.hyper,previous.weights);
        parameter_table=table(parameter_names',old_mean',batch_mean',seq_mean', ...
            batch_sd',seq_sd',abs(batch_mean-seq_mean)', ...
            'VariableNames',{'Parameter','Mean_1_3','Mean_batch','Mean_sequential','SD_batch','SD_sequential','Mean_difference'});
        state_table=state_summary(first);
        writetable(hyper_table,fullfile(output_dir,'posterior_hyperparameters.csv'));
        writetable(parameter_table,fullfile(output_dir,'posterior_parameters.csv'));
        writetable(state_table,fullfile(output_dir,'smoothed_states.csv'));
        selected=(1:min(N,12))';
        component_table=array2table([selected batch.weights(selected) sequential.weights(selected) ...
            batch.hyper(selected,:) sequential.hyper(selected,:)], ...
            'VariableNames',["Particle","Weight_batch","Weight_sequential", ...
            "Batch_"+hyper_names,"Sequential_"+hyper_names]);
        writetable(component_table,fullfile(output_dir,'example_mixture_components.csv'));
    end
end
diagnostic_names={'Repetition','Seed','ESS_T3','ESS_T4','Hyper_difference', ...
    'Weight_difference','Log_evidence_difference','Mean_difference','SD_difference','Largest_weight_T4'};
diagnostic_table=array2table(diagnostics,'VariableNames',diagnostic_names);
for k=1:6
    diagnostic_table.(['Mean_' char(parameter_names(k))])=parameter_means(:,k);
end
writetable(diagnostic_table,fullfile(output_dir,'replications.csv'));
result=struct('output_dir',output_dir,'particles',N,'repetitions',repetitions, ...
    'data_seed',data_seed,'prior',prior,'observations',Y, ...
    'hyperparameter_table',hyper_table,'parameter_table',parameter_table,'state_table',state_table, ...
    'diagnostic_table',diagnostic_table,'max_hyperparameter_difference',max(diagnostics(:,5)), ...
    'max_weight_difference',max(diagnostics(:,6)), ...
    'max_parameter_mean_difference',max(diagnostics(:,8)), ...
    'parameter_mean_sd_across_repetitions',std(parameter_means,0,1), ...
    'elapsed_seconds',toc(started));
metadata=rmfield(result,{'hyperparameter_table','parameter_table','state_table','diagnostic_table'});
write_json(fullfile(output_dir,'summary.json'),metadata);
write_report(result,first);
fprintf('\n%s\n',fullfile(output_dir,'REPORT.md'));
disp(hyper_table);
end

function [means,sds,names] = parameter_moments(hyper,weights)
% Integrate theta conditional on each latent path (Rao-Blackwellization).
nu=hyper(:,1); scale=hyper(:,2:4); denom=nu-3;
conditional_mean=scale./denom;
conditional_var=zeros(size(scale));
conditional_var(:,[1 3])=2*scale(:,[1 3]).^2./(denom.^2.*(nu-5));
conditional_var(:,2)=((nu-1).*scale(:,2).^2+(nu-3).*scale(:,1).*scale(:,3)) ...
    ./((nu-2).*(nu-3).^2.*(nu-5));
shapes=hyper(:,[5 5 8]); scales=hyper(:,[6 7 9]);
conditional_mean=[conditional_mean scales./(shapes-1)];
conditional_var=[conditional_var scales.^2./((shapes-1).^2.*(shapes-2))];
means=weights'*conditional_mean;
second=weights'*(conditional_var+conditional_mean.^2);
sds=sqrt(max(second-means.^2,0));
names=["V11","V12","V22","sig2h1","sig2h2","sig2r"];
end

function output = state_summary(first)
labels=["B1","B2","h1","h2","r","rho"];
names=strings(24,1); time=zeros(24,1); values=zeros(24,5);
for t=1:4
    xb=[first.batch.particles.B(:,:,t),first.batch.particles.h(:,:,t), ...
        first.batch.particles.r(:,t),tanh(first.batch.particles.r(:,t))];
    xs=[first.sequential.particles.B(:,:,t),first.sequential.particles.h(:,:,t), ...
        first.sequential.particles.r(:,t),tanh(first.sequential.particles.r(:,t))];
    mb=first.batch.weights'*xb; ms=first.sequential.weights'*xs;
    sb=sqrt(first.batch.weights'*(xb-mb).^2); ss=sqrt(first.sequential.weights'*(xs-ms).^2);
    before=first.previous.weights'*xs;
    rows=(t-1)*6+(1:6);
    names(rows)=labels'; time(rows)=t;
    values(rows,:)=[before' mb' ms' sb' ss'];
end
output=table(time,names,values(:,1),values(:,2),values(:,3),values(:,4),values(:,5), ...
    'VariableNames',{'Time','State','Mean_before_y4','Mean_batch','Mean_sequential','SD_batch','SD_sequential'});
end

function write_json(path,value)
fid=fopen(path,'w'); assert(fid>=0,'rt:Output','Cannot open output file.');
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

function write_report(result,first)
fid=fopen(fullfile(result.output_dir,'REPORT.md'),'w');
assert(fid>=0,'rt:Output','Cannot open report.'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Two-Series DSC: Batch Versus Sequential Bayes\n\n');
fprintf(fid,'Both routes agree to floating-point precision. This run used %d particles in each of %d paired repetitions, with the same four observations.\n\n',result.particles,result.repetitions);
fprintf(fid,'Maximum component-hyperparameter difference: `%.3g`. Maximum normalized-weight difference: `%.3g`. Maximum parameter-mean difference: `%.3g`.\n\n', ...
    result.max_hyperparameter_difference,result.max_weight_difference,result.max_parameter_mean_difference);
fprintf(fid,'## Data and Starting Prior\n\n');
fprintf(fid,'Data were simulated from the DSC prior with seed %d. Inference receives only the observations; the simulated latent truth is saved separately for inspection and is never treated as observed.\n\n',result.data_seed);
fprintf(fid,'| Date | y1 | y2 |\n|---|---:|---:|\n');
for t=1:4, fprintf(fid,'| %d | %.10f | %.10f |\n',t,result.observations(t,:)); end
fprintf(fid,'\n`B1 ~ N(mB, CB)`, independently of V. `V ~ IW(Psi, nu)`.\n\n');
fprintf(fid,'`sig2h_j ~ IG(a, bh)` and `sig2r ~ IG(a, br)`, with density proportional to `q^(-a-1) exp(-b/q)`.\n\n');
fprintf(fid,'`h1_j | sig2h_j ~ N(mh_j, (ch+1)*sig2h_j)`; `r1 | sig2r ~ N(mr, (cr+1)*sig2r)`. Later states follow random walks. `rho_t = tanh(r_t)` and `Sigma_t = diag(exp(h_t/2))*R_t*diag(exp(h_t/2))`.\n\n');
fprintf(fid,'All given hyperparameters, identical in both routes:\n\n```json\n%s\n```\n\n',jsonencode(result.prior,PrettyPrint=true));
fprintf(fid,'## What Is Compared\n\n');
fprintf(fid,'1. **Batch:** draw complete candidate paths from the original prior and weight them by `L1*L2*L3*L4`.\n');
fprintf(fid,'2. **Sequential:** draw paths through date 3, weight by `L1*L2*L3`, save and reload that joint posterior, append date-4 states from their transition distributions, and multiply the old weights by **L4 only**. `dsc_rt_update` accepts no old observations.\n');
fprintf(fid,'3. Use matched random numbers and no intermediate resampling. Thus the two routes have the same candidate paths and the same final weights, up to floating-point arithmetic. Old path weights change, so dates 1--3 are smoothed using observation 4.\n\n');
fprintf(fid,'The common prior particles are a controlled comparison, not independent Monte Carlo runs. The original continuous posterior is approximated. Exact equality of the paired approximations does not prove that finite-sample Monte Carlo error is zero. Different repetition seeds illustrate that uncertainty.\n\n');
fprintf(fid,'## Posterior Hyperparameters\n\n');
fprintf(fid,'Because states are latent, the observed-data posterior has **no single IW/IG hyperparameter vector**. We represent it by a weighted mixture: for each candidate path s,\n\n');
fprintf(fid,'- `V | path_s, y ~ IW(Psi_s, nu_s)`;\n- `sig2h_j | path_s, y ~ IG(a_h, b_hj_s)`;\n- `sig2r | path_s, y ~ IG(a_r, b_r_s)`.\n\n');
fprintf(fid,'Every component hyperparameter is checked between routes. The table reports the shared degrees of freedom/shapes and the **posterior-weighted mean of each scale**, using repetition 1. Scale averages do not turn the mixture into one IW/IG distribution.\n\n');
fprintf(fid,'| Hyperparameter | Prior | After 1:3 | Batch 1:4 | Sequential 3+1 | Max component difference |\n|---|---:|---:|---:|---:|---:|\n');
for k=1:height(result.hyperparameter_table)
    row=result.hyperparameter_table(k,:);
    fprintf(fid,'| %s | %.12f | %.12f | %.12f | %.12f | %.3g |\n', ...
        row.Hyperparameter,row.Prior,row.After_1_3,row.Batch_1_4,row.Sequential_3_then_4,row.Max_component_difference);
end
fprintf(fid,'\nIn particular, `nu_V: 8 -> 10 -> 11`, and `a_h = a_r: 6 -> 7.5 -> 8`. Each path-specific scale is updated by its new squared increment (or outer product for V). The h/r initial quadratic term is included.\n\n');
fprintf(fid,'Complete mixture components and weights are in `posterior_T4.mat`; twelve illustrative components are in `example_mixture_components.csv`. These are inferred candidate paths, not the simulated true path.\n\n');
fprintf(fid,'## Static Parameter Posterior Summaries\n\n');
fprintf(fid,'Means and standard deviations integrate each conditional IW/IG component analytically and then average over its posterior path weight. This reduces Monte Carlo noise and retains parameter uncertainty.\n\n');
fprintf(fid,'| Parameter | Mean after 1:3 | Batch mean | Sequential mean | Batch SD | Sequential SD |\n|---|---:|---:|---:|---:|---:|\n');
for k=1:height(result.parameter_table)
    row=result.parameter_table(k,:);
    fprintf(fid,'| %s | %.12f | %.12f | %.12f | %.12f | %.12f |\n', ...
        row.Parameter,row.Mean_1_3,row.Mean_batch,row.Mean_sequential,row.SD_batch,row.SD_sequential);
end
fprintf(fid,'\n## The Old History Also Updates\n\n');
fprintf(fid,'The first-date states below change after observation 4. Both routes give the same new smoothed means. All dates, including SDs, are in `smoothed_states.csv`.\n\n');
fprintf(fid,'| First-date state | Before y4 | Batch after y4 | Sequential after y4 |\n|---|---:|---:|---:|\n');
for k=1:6
    row=result.state_table(k,:);
    fprintf(fid,'| %s | %.12f | %.12f | %.12f |\n',row.State,row.Mean_before_y4,row.Mean_batch,row.Mean_sequential);
end
fprintf(fid,'\nFor date 4, the `Mean_before_y4` column in the CSV is a predictive mean; for dates 1--3 it is the old smoothed mean. State labels B1/B2 and h1/h2 identify series; the Time column identifies the date.\n\n');
fprintf(fid,'## Monte Carlo Diagnostics\n\n');
fprintf(fid,'ESS here means importance-weight ESS, `1/sum(w.^2)`, not MCMC autocorrelation ESS.\n\n');
fprintf(fid,'| Repeat | ESS after 1:3 | ESS after 1:4 | Max hyperparameter difference | Max weight difference |\n|---|---:|---:|---:|---:|\n');
for k=1:height(result.diagnostic_table)
    row=result.diagnostic_table(k,:);
    fprintf(fid,'| %d | %.1f | %.1f | %.3g | %.3g |\n',row.Repetition,row.ESS_T3,row.ESS_T4,row.Hyper_difference,row.Weight_difference);
end
fprintf(fid,'\nAcross-repetition SD of the estimated means (V11, V12, V22, sig2h1, sig2h2, sig2r):\n\n`%s`\n\n',mat2str(result.parameter_mean_sd_across_repetitions,6));
fprintf(fid,'This is empirical Monte Carlo variation for one fixed data set, not a posterior SD and not a study of repeated data sets. Each pair of routes shares this numerical error.\n\n');
fprintf(fid,'Elapsed time: %.2f seconds. Both evidence calculations agree: direct `mean(L1*L2*L3*L4)` and old evidence times the sequential predictive density.\n\n',result.elapsed_seconds);
fprintf(fid,'The final weighted posterior is the same, with no refitting of hyperpriors, no use of latent truth, and no repeated old-data likelihood in the sequential update. This is a small teaching exercise, not an alteration of the production sampler.\n');
assert(max(abs(first.batch.weights-first.sequential.weights))<1e-11);
end

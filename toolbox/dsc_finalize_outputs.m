function result=dsc_finalize_outputs(result,cfg,manifest,calibration_panel)
%DSC_FINALIZE_OUTPUTS Publish the run manifest, parameters, draws and path PDFs.
root=canonical(cfg.output_root); prefix=char(cfg.file_prefix);
assert(isfolder(root)&&~isempty(regexp(prefix,'^[A-Za-z0-9_-]+_$','once')), ...
    'dsc:FinalizeIdentity','An existing output directory and valid run prefix are required.');
manifest_path=fullfile(root,[prefix 'run_manifest.json']);
if nargin<3, manifest=jsondecode(fileread(manifest_path)); end
if nargin>=4&&~isempty(calibration_panel)
    validate_calibration_panel(calibration_panel,result,cfg);
    manifest.return_ar1=dsc_return_ar1(calibration_panel);
end
parameter_path=fullfile(root,[prefix 'parameters.mat']);
source_parameter_file=char(result.paths.parameters);
has_own_bundle=isfile(source_parameter_file)&&owned(source_parameter_file,root,prefix);
stage=[];
stage_path=fullfile(root,[prefix 'estimation_result.mat']);
if isfile(stage_path)
    loaded=load(stage_path,'result'); stage=loaded.result;
end

% Validate every retained draw against its chain and iteration before any
% file is replaced. The chunks remain unchanged: parameters and paths form
% a single joint posterior draw and have one authoritative storage location.
if has_own_bundle
    loaded=load(source_parameter_file,'parameters'); parameters=loaded.parameters;
    parameter_chunks=validate_chunks(parameters,root,prefix);
    report=parameters;
    report.storage='parameter_report_v1';
    report.calibration_sample_sha256=sample_fingerprint(parameters,'calibration');
    report.training_sample_sha256=sample_fingerprint(parameters,'training');
    report.posterior_std=parameter_std(parameters.parameter_draws);
    if isempty(stage)
        report.prior_manifest_section='priors';
    else
        report.prior_manifest_section='parameter_estimation';
        manifest.parameter_estimation=struct('priors',json_dates(parameters.priors), ...
            'run_summary',run_summary(stage));
    end
    report=remove_fields(report,{'priors','parameter_draws', ...
        'training_returns','training_mask','calibration_returns','calibration_mask'});
    report.parameter_chunk_files=parameter_chunks;
else
    % Conditional smoothing produces new mu/state draws. These summaries are
    % explicitly not reusable parameter-estimation bundles.
    [draws,parameter_chunks]=current_parameter_draws(result,root,prefix);
    report=struct('storage','conditional_parameter_summary_v1', ...
        'model','joint_dsc_constant_mean','mean_model','constant', ...
        'source_parameter_file',source_parameter_file, ...
        'parameter_estimate',parameter_mean(draws), ...
        'posterior_std',parameter_std(draws), ...
        'retained_draws',numel(draws.chain_id), ...
        'parameter_chunk_files',{parameter_chunks}, ...
        'prior_manifest_section','priors');
    report.prior_posterior_summary=dsc_parameter_summary(draws,result.priors, ...
        result.tickers,result.pair_i,result.pair_j);
    report.parameter_summary_metadata=struct( ...
        'posterior_source','Retained conditional-smoothing draws.', ...
        'posterior_std_convention','Sample standard deviation with n-1 denominator.');
end
report.manifest_file=manifest_path;
manifest.run_summary=run_summary(result);
if isfield(result,'inference')
    manifest.inference=remove_fields(result.inference,{'estimation_convergence','state_convergence','parameter_smoothing'});
    if has_own_bundle, manifest.inference.parameter_file=parameter_path; end
end
draw_chunks=owned_files(root,[prefix '*posterior_chunk_*.mat'],prefix);
% Complete all publications before deleting staging files. Each replacement
% is atomic; failure leaves the original joint chunks and remaining sidecars.
manifest=dsc_run_lifecycle('publishing',root,prefix,manifest);
atomic_mat(parameter_path,'parameters',report,root,prefix);

patterns={'*checkpoint.mat','*pilot_summary.json','*prior_summary.json', ...
    '*calibrated_priors.mat','*inference_metadata.json', ...
    '*posterior_correlation_bands.mat', ...
    '*chain_*diagnostics.mat','*chain_*result.mat','estimation_result.mat', ...
    'result.mat','*.pdf'};
for p=1:numel(patterns)
    files=owned_files(root,[prefix patterns{p}],prefix);
    for k=1:numel(files), delete(files{k}); end
end
% A cutoff bundle was saved with an estimation_ prefix; publish its report
% under the invocation's single prefix and remove the now redundant bundle.
if has_own_bundle&&~strcmp(canonical(source_parameter_file),canonical(parameter_path))
    delete(source_parameter_file);
end
result.paths=struct('manifest',manifest_path,'parameters',parameter_path, ...
    'draw_chunks',{draw_chunks});
result.figure_paths=strings(0,1); result.figure_pdf=''; result.report_pdf='';
create_pdf=true;
if isfield(cfg,'create_pdf'), create_pdf=cfg.create_pdf; end
% Plot from the published joint chunks and manifest, without retaining a
% duplicate bands MAT file. Cleanup has finished, so these PDFs are retained.
if create_pdf
    result.figure_paths=plot_bayes_correlation_paths(root,prefix);
    result.figure_pdf=char(result.figure_paths(1));
    result.paths.figures=result.figure_paths;
end
dsc_run_lifecycle('finish',root,prefix,manifest);
% Only fully published invocations can trigger retention. Failed estimation
% and publication leave evidence for the next successful run to classify.
try
    dsc_prune_outputs(root);
catch problem
    warning('dsc:OutputRetention','Run saved, but output retention failed: %s',problem.message);
end
end

function validate_calibration_panel(panel,result,cfg)
required={'dates','tickers','returns','observation_mask','pair_i','pair_j','source_hash'};
assert(isstruct(panel)&&isscalar(panel)&&all(isfield(panel,required))&& ...
    all(isfield(result,{'dates','tickers','m','pair_i','pair_j','priors','inference'}))&& ...
    isstruct(result.inference)&&isscalar(result.inference)&& ...
    all(isfield(result.inference,{'calibration_weeks','calibration_start','calibration_end'}))&& ...
    isfield(cfg,'prior_weeks'), ...
    'dsc:ReturnAR1Identity','The initial prior-calibration panel and its run identity are required for the return AR(1) appendix.');
weeks=cfg.prior_weeks; dates=panel.dates(:);
assert(isnumeric(weeks)&&isreal(weeks)&&isscalar(weeks)&&isfinite(weeks)&& ...
    weeks>=1&&weeks==floor(weeks)&&isequal(result.inference.calibration_weeks,weeks)&& ...
    isdatetime(dates)&&numel(dates)==weeks&&~any(isnat(dates))&& ...
    all(dates(2:end)==dates(1:end-1)+calweeks(1))&& ...
    isdatetime(result.dates)&&~isempty(result.dates)&& ...
    dates(end)+calweeks(1)==result.dates(1)&& ...
    isequal(string(panel.tickers(:)),string(result.tickers(:)))&& ...
    isequal(size(panel.returns),[weeks result.m])&& ...
    isequal(size(panel.observation_mask),size(panel.returns))&& ...
    isequal(panel.pair_i(:),result.pair_i(:))&& ...
    isequal(panel.pair_j(:),result.pair_j(:))&& ...
    isfield(result.priors,'source_hash')&&isequal(panel.source_hash,result.priors.source_hash), ...
    'dsc:ReturnAR1Identity','The return AR(1) appendix must use exactly the initial prior-calibration weeks preceding the model, with matching series, pair order, mask dimensions and source.');
expected_dates=string(dates([1 end]),'yyyy-MM-dd');
assert(isequal(expected_dates,[string(result.inference.calibration_start);string(result.inference.calibration_end)]), ...
    'dsc:ReturnAR1Identity','The return AR(1) dates differ from the run''s initial prior-calibration dates.');
for field={'calibration_dates','initial_dates'}
    name=field{1};
    if isfield(result.priors,name)
        value=result.priors.(name);
        if isdatetime(value), value=string(value(:),'yyyy-MM-dd');
        else, value=string(value(:)); end
        assert(isequal(value,expected_dates),'dsc:ReturnAR1Identity', ...
            'The return AR(1) dates differ from the saved %s.',name);
    end
end
if isfield(result.priors,'initial_weeks')
    assert(isequal(result.priors.initial_weeks,weeks),'dsc:ReturnAR1Identity', ...
        'The return AR(1) week count differs from the initial prior-calibration window.');
end
end

function paths=validate_chunks(parameters,root,prefix)
assert(all(isfield(parameters,{'parameter_draws','chain_ids','chain_file_prefixes','retained_draws_per_chain'})), ...
    'dsc:FinalizeParameters','Parameter bundle lacks chain/draw identity.');
draws=parameters.parameter_draws;
assert(all(isfield(draws,{'chain_id','iteration','mu','sig2h','sig2r'})), ...
    'dsc:FinalizeParameters','Incomplete parameter draws.');
keys=[draws.chain_id(:) draws.iteration(:)]; n=size(keys,1);
assert(size(unique(keys,'rows'),1)==n&&all(isfinite(keys),'all')&& ...
    all(keys==floor(keys),'all'),'dsc:FinalizeParameters','Parameter draw identities must be unique integers.');
for field={'mu','sig2h','sig2r'}
    assert(size(draws.(field{1}),1)==n,'dsc:FinalizeParameters','Parameter draw counts disagree.');
end
ids=parameters.chain_ids(:); prefixes=cellstr(string(parameters.chain_file_prefixes));
counts=parameters.retained_draws_per_chain(:); paths={}; seen=false(n,1);
assert(numel(ids)==numel(prefixes)&&numel(counts)==numel(ids)&&numel(unique(ids))==numel(ids), ...
    'dsc:FinalizeParameters','Invalid parameter chain identities/counts.');
for c=1:numel(ids)
    chain_prefix=prefixes{c};
    assert(startsWith(chain_prefix,prefix)&&~isempty(regexp(chain_prefix,'^[A-Za-z0-9_-]+_$','once')), ...
        'dsc:FinalizeIdentity','A parameter source chain is outside this run.');
    files=owned_files(root,[chain_prefix 'posterior_chunk_*.mat'],prefix);
    expected=find(draws.chain_id(:)==ids(c)); total=0;
    assert(numel(expected)==counts(c),'dsc:FinalizeParameters','Saved chain count disagrees with parameter draws.');
    for k=1:numel(files)
        loaded=load(files{k},'chunk');
        assert(isfield(loaded,'chunk'),'dsc:FinalizeChunks','Missing chunk in %s.',files{k});
        chunk=loaded.chunk;
        assert(all(isfield(chunk,{'chain_id','iterations','mu','sig2h','sig2r','h','P_pairs'}))&& ...
            isequal(chunk.chain_id,ids(c)),'dsc:FinalizeChunks','Invalid chunk chain identity.');
        iterations=chunk.iterations(:); m=numel(iterations);
        assert(m>0&&numel(unique(iterations))==m,'dsc:FinalizeChunks','Invalid chunk iteration identities.');
        [found,rows]=ismember([repmat(ids(c),m,1) iterations],keys,'rows');
        assert(all(found)&&~any(seen(rows)),'dsc:FinalizeChunks','Chunk iterations are missing or duplicated in the parameter bundle.');
        for field={'mu','sig2h','sig2r'}
            key=field{1};
            assert(isequaln(double(chunk.(key)),double(draws.(key)(rows,:))), ...
                'dsc:FinalizeChunks','Chunk %s differs from saved parameter draws (%s).',files{k},key);
        end
        seen(rows)=true; total=total+m; paths{end+1}=files{k}; %#ok<AGROW>
    end
    assert(total==numel(expected),'dsc:FinalizeChunks','Chunks do not cover all saved chain draws.');
end
assert(all(seen),'dsc:FinalizeChunks','Chunks do not cover the parameter bundle.');
end

function hash=sample_fingerprint(parameters,kind)
block=struct('dates',parameters.([kind '_dates']), ...
    'returns',parameters.([kind '_returns']), ...
    'observation_mask',parameters.([kind '_mask']), ...
    'tickers',parameters.tickers,'pair_i',parameters.pair_i,'pair_j',parameters.pair_j);
hash=dsc_sample_fingerprint(block);
end

function [draws,paths]=current_parameter_draws(result,root,prefix)
paths=owned_files(root,[prefix 'chain_*posterior_chunk_*.mat'],prefix);
draws=struct('mu',[],'sig2h',[],'sig2r',[],'chain_id',[],'iteration',[]);
for k=1:numel(paths)
    loaded=load(paths{k},'chunk'); chunk=loaded.chunk; n=numel(chunk.iterations);
    assert(all(isfield(chunk,{'mu','sig2h','sig2r','chain_id','iterations'})), ...
        'dsc:FinalizeChunks','Incomplete conditional draw chunk.');
    for field={'mu','sig2h','sig2r'}
        key=field{1};
        assert(size(chunk.(key),1)==n,'dsc:FinalizeChunks','Conditional draw dimensions disagree.');
        draws.(key)=[draws.(key);double(chunk.(key))];
    end
    draws.chain_id=[draws.chain_id;repmat(chunk.chain_id,n,1)];
    draws.iteration=[draws.iteration;chunk.iterations(:)];
end
keys=[draws.chain_id draws.iteration];
assert(size(keys,1)==result.saved_draws&&size(unique(keys,'rows'),1)==size(keys,1), ...
    'dsc:FinalizeChunks','Conditional draw count or identity mismatch.');
end

function value=parameter_mean(draws)
value=struct('mu',mean(draws.mu,1),'sig2h',mean(draws.sig2h,1),'sig2r',mean(draws.sig2r,1));
end

function value=parameter_std(draws)
value=struct('mu',std(draws.mu,0,1),'sig2h',std(draws.sig2h,0,1),'sig2r',std(draws.sig2r,0,1));
end

function summary=run_summary(result)
summary=struct();
names={'total_seconds','session_seconds','stage_seconds','stage_elapsed_seconds'};
for k=1:numel(names)
    key=names{k}; if isfield(result,key), summary.(key)=result.(key); end
end
end

function value=json_dates(value)
if ~isstruct(value), return; end
for key={'initial_dates','calibration_dates','likelihood_start_date'}
    name=key{1};
    if isfield(value,name)&&isdatetime(value.(name))
        dates=string(value.(name),'yyyy-MM-dd');
        if isscalar(dates), value.(name)=char(dates);
        else, value.(name)=cellstr(dates); end
    end
end
end

function value=remove_fields(value,names)
if isstruct(value)
    present=intersect(fieldnames(value),names);
    if ~isempty(present), value=rmfield(value,present); end
end
end

function paths=owned_files(root,pattern,prefix)
files=dir(fullfile(root,pattern)); paths={};
[~,order]=sort({files.name}); files=files(order);
for k=1:numel(files)
    if files(k).isdir, continue; end
    path=fullfile(files(k).folder,files(k).name);
    assert(owned(path,root,prefix),'dsc:FinalizeIdentity','Output path escapes the current run.');
    paths{end+1}=path; %#ok<AGROW>
end
end

function yes=owned(path,root,prefix)
path=canonical(path); [directory,name,extension]=fileparts(path);
yes=strcmp(directory,root)&&startsWith([name extension],prefix);
end

function path=canonical(path)
path=char(java.io.File(char(path)).getCanonicalPath());
end

function atomic_mat(path,name,value,root,prefix)
temporary=temporary_path(root,prefix,'.mat');
guard=onCleanup(@()remove_temporary(temporary));
data=struct(); data.(name)=value;
save(temporary,'-struct','data','-v7.3');
publish(temporary,path); clear guard
end

function path=temporary_path(root,prefix,extension)
path=fullfile(root,[prefix 'finalizing_' char(java.util.UUID.randomUUID()) extension]);
end

function publish(temporary,path)
[ok,message]=movefile(temporary,path,'f');
assert(ok,'dsc:FinalizeWrite','Cannot publish %s: %s',path,message);
end

function remove_temporary(path)
if isfile(path), delete(path); end
end

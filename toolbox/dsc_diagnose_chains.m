function [report,details] = dsc_diagnose_chains(chain_dirs,cfg,output_dir)
%DSC_DIAGNOSE_CHAINS Diagnose retained chunks, preserving original chains.
% No burn-in is removed here: posterior_chunk files already contain only
% retained draws. Unequal chains use their earliest common retained length.
% Actual P_pairs and h state paths are checked at EVERY modeled date.
% The latent r coordinates, B states and hyperparameters remain available
% in saved chunks but are excluded from the convergence verdict.
% A pass requires at least four original chains, enough retained draws,
% and finite passing diagnostics for every checked scalar quantity.
%
% Current sampler chunks contain a scalar struct, whose numeric fields
% cannot be sliced directly with matfile. A single-chunk conversion pass
% creates temporary top-level arrays, then reads bounded blocks with matfile.
% Peak memory is one source chunk plus one diagnostic block, not all states
% across all retained draws. Temporary staging is removed on exit.
started=tic;
if ischar(chain_dirs)||isstring(chain_dirs), chain_dirs=cellstr(chain_dirs); end
if ~iscell(chain_dirs)||~all(cellfun(@(x)ischar(x)||(isstring(x)&&isscalar(x)),chain_dirs))
    error('dsc:DiagnosticConfig','chain_dirs must contain directory paths.');
end
chain_dirs=cellfun(@char,chain_dirs(:)','UniformOutput',false);
if nargin<2||isempty(cfg), cfg=struct; end
if nargin<3||isempty(output_dir)
    if isempty(chain_dirs), error('dsc:DiagnosticConfig','An output directory is required.'); end
    output_dir=chain_dirs{1};
end
cfg=diagnostic_config(cfg);
if ~isfolder(output_dir), mkdir(output_dir); end
thresholds=struct('rhat_threshold',cfg.rhat_threshold,'min_ess',cfg.min_ess, ...
    'max_mcse_ratio',cfg.max_mcse_ratio,'min_diagnostic_draws',cfg.min_diagnostic_draws, ...
    'min_original_chains',4,'rhat',cfg.rhat_threshold,'ess',cfg.min_ess, ...
    'mcse_ratio',cfg.max_mcse_ratio,'min_draws',cfg.min_diagnostic_draws);
scope='P_pairs and h state paths at every modeled date; r, B and hyperparameters excluded.';
report=struct('status','not_checked','reason','','passed',false,'chains',numel(chain_dirs), ...
    'draws_per_chain',0,'actual_draws_per_chain',zeros(1,numel(chain_dirs)), ...
    'quantities_checked',0,'max_rhat',NaN,'min_ess_bulk',NaN,'min_ess_tail',NaN, ...
    'max_mcse_sd_ratio',NaN,'thresholds',thresholds,'scope',scope, ...
    'failing_quantities',{{}},'unavailable_quantities',{{}}, ...
    'chain_ids',nan(1,numel(chain_dirs)),'chain_seeds',nan(1,numel(chain_dirs)), ...
    'chain_dirs',{chain_dirs}, ...
    'chunk_metadata',{cell(1,numel(chain_dirs))}, ...
    'retained_iteration_ranges',nan(numel(chain_dirs),2), ...
    'alignment','Earliest common retained draw count; no additional burn-in removal.', ...
    'method','Rank-normalized split/folded Rhat; Geyer positive/monotone bulk/tail ESS; raw-mean MCSE.', ...
    'mode',cfg.convergence_mode);
report.paths=struct('json',fullfile(output_dir,'convergence_diagnostics.json'), ...
    'csv',fullfile(output_dir,'convergence_diagnostics.csv'), ...
    'mat',fullfile(output_dir,'convergence_diagnostics.mat'));
details=empty_details();
if strcmp(cfg.convergence_mode,'off')
    % Counts are unknown because off intentionally does not inspect chunks.
    report=rmfield(report,'actual_draws_per_chain');
    report.reason='Convergence checking is disabled.';
    write_outputs(report,details); return
end
if ~cfg.estimate_parameters&&isfield(cfg,'parameter_smoothing')&&strcmp(cfg.parameter_smoothing,'draws')
    report=rmfield(report,'actual_draws_per_chain');
    report.reason='Cycling saved parameter draws changes the conditional transition kernel; these draws cannot establish convergence to a single fixed target.';
    write_outputs(report,details); return
end
counts_known=false;
try
check_budget();
chain_dirs=cellfun(@(p)char(java.io.File(p).getCanonicalPath()),chain_dirs,'UniformOutput',false);
if numel(unique(chain_dirs))~=numel(chain_dirs)
    error('dsc:DiagnosticIdentity','Original chain directories must be distinct.');
end
report.chain_dirs=chain_dirs;
stage_dir=tempname(output_dir); mkdir(stage_dir);
stage_cleanup=onCleanup(@()remove_staging(stage_dir)); %#ok<NASGU>
staged=cell(1,numel(chain_dirs)); identity=[]; target=[];
for c=1:numel(chain_dirs)
    check_budget();
    chunks=ordered_chunks(chain_dirs{c});
    if isempty(chunks), continue; end
    [current,provenance]=validated_checkpoint(chain_dirs{c},cfg.estimate_parameters);
    check_budget();
    if isempty(target), target=current;
    elseif ~isequaln(current,target)
        error('dsc:DiagnosticIdentity','Original chains must have the same data, priors, fixed parameters and configuration.');
    end
    stage=[];
    count=0; previous=provenance.burnin; chain_id=provenance.chain_id;
    for k=1:numel(chunks)
        check_budget();
        loaded=load(chunks{k},'chunk');
        if ~isfield(loaded,'chunk'), error('dsc:DiagnosticChunk','Missing chunk in %s.',chunks{k}); end
        chunk=loaded.chunk; clear loaded;
        check_budget();
        meta=validate_chunk(chunk,cfg.estimate_parameters,chunks{k});
        if chunk.chain_id~=chain_id||~isequal(chunk.dates(:),provenance.dates(:))|| ...
                ~isequal(string(chunk.tickers(:)),string(provenance.tickers(:)))
            error('dsc:DiagnosticIdentity','Chunk chain identity or sample differs from its checkpoint.');
        end
        if isempty(identity), identity=meta; end
        if ~isequaln(meta,identity)
            error('dsc:DiagnosticIdentity','Dates, tickers, pair order or state dimensions differ in %s.',chunks{k});
        end
        if k==1
            stage=matfile(fullfile(stage_dir,sprintf('chain_%d.mat',c)),'Writable',true);
            stage.h=zeros(provenance.saved,meta.T*meta.m);
            stage.P_pairs=zeros(provenance.saved,meta.T*meta.nr);
        end
        if k==1
            report.retained_iteration_ranges(c,1)=chunk.iterations(1);
        end
        iterations=double(chunk.iterations(:));
        if any(diff([previous;iterations])~=provenance.thin)||iterations(end)>provenance.completed
            error('dsc:DiagnosticIterations','Retained iterations do not match the checkpoint retention schedule in %s.',chunks{k});
        end
        previous=iterations(end); number=numel(iterations); rows=count+(1:number);
        file_info=dir(chunks{k});
        record=struct('file',chunks{k},'retained_draws',number, ...
            'first_iteration',iterations(1),'last_iteration',iterations(end), ...
            'bytes',file_info.bytes,'modified_datenum',file_info.datenum);
        if k==1, report.chunk_metadata{c}=record; else, report.chunk_metadata{c}(k)=record; end
        % The source is already retained, regardless of cfg.burnin/cfg.thin.
        for family={'h','P_pairs'}
            check_budget();
            name=family{1}; values=chunk.(name);
            if strcmp(name,'h')
                stage.h(rows,:)=reshape(values,[],number)';
            else
                stage.P_pairs(rows,:)=reshape(values,[],number)';
            end
        end
        count=count+number;
    end
    if count~=provenance.saved
        error('dsc:DiagnosticIterations','Saved retained chunk counts disagree with the source checkpoint.');
    end
    staged{c}=stage;
    report.actual_draws_per_chain(c)=count;
    report.chain_ids(c)=chain_id;
    report.chain_seeds(c)=provenance.seed;
    report.retained_iteration_ranges(c,2)=previous;
end
counts_known=true;
clear chunk values stage;
ids=report.chain_ids(isfinite(report.chain_ids));
seeds=report.chain_seeds(isfinite(report.chain_seeds));
if numel(unique(ids))~=numel(ids)||numel(unique(seeds))~=numel(seeds)
    error('dsc:DiagnosticIdentity','Original chain IDs and effective seeds must be distinct across chain directories.');
end
if isempty(chain_dirs), aligned=0; else, aligned=min(report.actual_draws_per_chain); end
report.draws_per_chain=aligned;
if aligned==0
    report.status='insufficient_draws';
    report.reason='At least one original chain has no saved retained draws.';
    write_outputs(report,details); return
end
details=quantity_details(identity);
metrics=nan(height(details),5); available=false(height(details),1);
passed=false(height(details),1); reasons=strings(height(details),1);
% Limit the draw cube to approximately 32 MiB (and at most 64 scalars).
block_size=max(1,min(64,floor(32*2^20/(8*aligned*numel(chain_dirs)))));
families=unique(details.family,'stable');
for f=1:numel(families)
    family=char(families(f)); indices=find(details.family==family);
    for first=1:block_size:numel(indices)
        check_budget();
        block=first:min(first+block_size-1,numel(indices));
        draws=zeros(aligned,numel(chain_dirs),numel(block));
        for c=1:numel(chain_dirs)
            check_budget();
            if strcmp(family,'h')
                values=staged{c}.h(1:aligned,block);
            else
                values=staged{c}.P_pairs(1:aligned,block);
            end
            draws(:,c,:)=reshape(values,aligned,1,numel(block));
        end
        for b=1:numel(block)
            check_budget();
            row=indices(block(b)); d=dsc_mcmc_diagnostics(draws(:,:,b));
            metrics(row,:)=[d.rhat d.ess_bulk d.ess_tail d.mcse_mean d.mcse_sd_ratio];
            available(row)=d.available; reasons(row)=string(d.reason);
            passed(row)=d.available&&d.rhat<cfg.rhat_threshold&& ...
                d.ess_bulk>=cfg.min_ess&&d.ess_tail>=cfg.min_ess&& ...
                d.mcse_sd_ratio<=cfg.max_mcse_ratio;
        end
    end
end
check_budget();
details{:,{'rhat','ess_bulk','ess_tail','mcse_mean','mcse_sd_ratio'}}=metrics;
details.available=available; details.passed=passed; details.reason=reasons;
report.quantities_checked=height(details);
report.max_rhat=finite_extreme(details.rhat,'max');
report.min_ess_bulk=finite_extreme(details.ess_bulk,'min');
report.min_ess_tail=finite_extreme(details.ess_tail,'min');
report.max_mcse_sd_ratio=finite_extreme(details.mcse_sd_ratio,'max');
report.failing_quantities=cellstr(details.quantity(details.available&~details.passed));
report.unavailable_quantities=cellstr(details.quantity(~details.available));
if aligned<cfg.min_diagnostic_draws
    report.status='insufficient_draws';
    report.reason='The aligned retained count is below min_diagnostic_draws.';
elseif numel(chain_dirs)<4
    report.status='insufficient_chains';
    report.reason='At least four distinct original chains are required; split halves are not original chains.';
elseif any(report.actual_draws_per_chain~=aligned)
    report.status='failed';
    report.reason='Unequal retained counts: aligned diagnostics are exploratory because pooled summaries use all retained draws.';
elseif ~isempty(details)&&all(details.passed)
    report.status='passed'; report.passed=true;
    report.reason='Every checked quantity satisfies the configured diagnostic thresholds.';
else
    report.status='failed';
    report.reason='At least one checked quantity has unavailable or failing diagnostics.';
end
write_outputs(report,details);
catch problem
    if ~strcmp(problem.identifier,'dsc:DiagnosticTimeLimit'), rethrow(problem); end
    report.status='not_checked'; report.passed=false;
    report.reason='Diagnostic time budget exhausted; convergence has not been established.';
    if ~counts_known, report=rmfield(report,'actual_draws_per_chain'); end
    details=empty_details(); report.quantities_checked=0;
    write_outputs(report,details);
end

    function check_budget()
        if toc(started)>=cfg.diagnostic_max_seconds
            error('dsc:DiagnosticTimeLimit','Diagnostic time budget exhausted.');
        end
    end
end

function cfg=diagnostic_config(cfg)
defaults=struct('convergence_mode','report','rhat_threshold',1.01, ...
    'min_ess',400,'max_mcse_ratio',.05,'min_diagnostic_draws',100, ...
    'estimate_parameters',true,'diagnostic_max_seconds',Inf);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(cfg,names{k}), cfg.(names{k})=defaults.(names{k}); end
end
cfg.convergence_mode=validatestring(cfg.convergence_mode,{'off','report'});
if ~isnumeric(cfg.diagnostic_max_seconds)||~isreal(cfg.diagnostic_max_seconds)|| ...
        ~isscalar(cfg.diagnostic_max_seconds)||isnan(cfg.diagnostic_max_seconds)||cfg.diagnostic_max_seconds<0
    error('dsc:DiagnosticConfig','diagnostic_max_seconds must be nonnegative (or Inf).');
end
for name={'rhat_threshold','min_ess','max_mcse_ratio','min_diagnostic_draws'}
    x=cfg.(name{1});
    if ~isnumeric(x)||~isreal(x)||~isscalar(x)||~isfinite(x)||x<=0
        error('dsc:DiagnosticConfig','%s must be a positive finite scalar.',name{1});
    end
end
if cfg.min_diagnostic_draws~=floor(cfg.min_diagnostic_draws)||cfg.min_diagnostic_draws<6
    error('dsc:DiagnosticConfig','min_diagnostic_draws must be an integer of at least 6.');
end
if ~isscalar(cfg.estimate_parameters)||~ismember(cfg.estimate_parameters,[false true])
    error('dsc:DiagnosticConfig','estimate_parameters must be true or false.');
end
end

function files=ordered_chunks(folder)
entries=dir(fullfile(folder,'posterior_chunk_*.mat'));
files=cell(1,numel(entries)); numbers=zeros(1,numel(entries));
for k=1:numel(entries)
    token=regexp(entries(k).name,'^posterior_chunk_(\d+)\.mat$','tokens','once');
    if isempty(token), error('dsc:DiagnosticChunk','Invalid posterior chunk filename.'); end
    numbers(k)=str2double(token{1}); files{k}=fullfile(folder,entries(k).name);
end
if numel(unique(numbers))~=numel(numbers)
    error('dsc:DiagnosticChunk','Duplicate numeric posterior chunk indices.');
end
[~,order]=sort(numbers); files=files(order);
end

function [target,p]=validated_checkpoint(folder,estimate)
path=fullfile(folder,'checkpoint.mat');
if ~isfile(path), error('dsc:DiagnosticIdentity','A source checkpoint is required to verify the chain target: %s.',path); end
loaded=load(path,'checkpoint');
if ~isfield(loaded,'checkpoint')||~isstruct(loaded.checkpoint)||~isscalar(loaded.checkpoint)|| ...
        ~all(isfield(loaded.checkpoint,{'identity','saved','completed'}))
    error('dsc:DiagnosticIdentity','Incomplete source checkpoint.');
end
cp=loaded.checkpoint; identity=cp.identity;
if ~isstruct(identity)||~isscalar(identity)|| ...
        ~all(isfield(identity,{'dates','tickers','returns','mask','source_hash','cfg','priors'}))|| ...
        ~isstruct(identity.cfg)||~isscalar(identity.cfg)|| ...
        ~all(isfield(identity.cfg,{'seed','chain_id','burnin','thin','estimate_parameters'}))
    error('dsc:DiagnosticIdentity','Checkpoint target identity or retention settings are missing.');
end
settings=identity.cfg;
numbers={settings.seed,settings.chain_id,settings.burnin,settings.thin,cp.saved,cp.completed};
if ~all(cellfun(@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x==floor(x),numbers))|| ...
        settings.seed<0||settings.chain_id<1||settings.burnin<0||settings.thin<1||cp.saved<0||cp.completed<0|| ...
        ~isequal(logical(settings.estimate_parameters),logical(estimate))|| ...
        ~isequal(size(identity.returns),[numel(identity.dates) numel(identity.tickers)])|| ...
        ~isequal(size(identity.mask),size(identity.returns))
    error('dsc:DiagnosticIdentity','Checkpoint has invalid mode, target dimensions, seeds or retention counts.');
end
p=struct('chain_id',settings.chain_id,'seed',settings.seed+settings.chain_id-1, ...
    'burnin',settings.burnin,'thin',settings.thin,'saved',cp.saved,'completed',cp.completed, ...
    'dates',identity.dates,'tickers',identity.tickers);
target=identity;
for field={'seed','chain_id','max_iterations','max_seconds','resume','output_root', ...
        'resume_validate_only','resume_finalize_only','diagnostic_max_seconds'}
    if isfield(target.cfg,field{1}), target.cfg=rmfield(target.cfg,field{1}); end
end
end

function meta=validate_chunk(chunk,estimate,path)
required={'B','h','P_pairs','dates','tickers','pair_i','pair_j','chain_id','iterations'};
if estimate, required=[required {'V','sig2h','sig2r'}]; end
if ~isstruct(chunk)||~isscalar(chunk)||~all(isfield(chunk,required))
    error('dsc:DiagnosticChunk','Required retained draw fields are missing in %s.',path);
end
iterations=chunk.iterations;
if ~isnumeric(iterations)||~isreal(iterations)||~isvector(iterations)||isempty(iterations)|| ...
        any(~isfinite(iterations))||any(iterations<1)||any(iterations~=floor(iterations))||any(diff(iterations)<=0)
    error('dsc:DiagnosticIterations','Retained iterations must increase strictly in %s.',path);
end
if ~isnumeric(chunk.chain_id)||~isscalar(chunk.chain_id)||~isreal(chunk.chain_id)|| ...
        ~isfinite(chunk.chain_id)||chunk.chain_id<1||chunk.chain_id~=floor(chunk.chain_id)
    error('dsc:DiagnosticIdentity','A positive integer original chain_id is required.');
end
T=numel(chunk.dates); m=numel(chunk.tickers); nr=m*(m-1)/2; n=numel(iterations);
if ~isdatetime(chunk.dates)||T<1||any(isnat(chunk.dates))||any(diff(chunk.dates)<=seconds(0))||m<2
    error('dsc:DiagnosticIdentity','Modeled dates must be valid and strictly increasing.');
end
[pi,pj]=find(tril(true(m),-1));
if ~isequal(chunk.pair_i(:),pi)||~isequal(chunk.pair_j(:),pj)
    error('dsc:DiagnosticIdentity','Correlation pairs must follow strict lower-triangle order.');
end
check_shape(chunk.B,[T m n],path,'B'); check_shape(chunk.h,[T m n],path,'h');
check_shape(chunk.P_pairs,[T nr n],path,'P_pairs');
if estimate
    check_shape(chunk.V,[m m n],path,'V');
    check_shape(chunk.sig2h,[n m],path,'sig2h');
    check_shape(chunk.sig2r,[n nr],path,'sig2r');
end
meta=struct('T',T,'m',m,'nr',nr,'dates',chunk.dates(:), ...
    'tickers',string(chunk.tickers(:)),'pair_i',pi,'pair_j',pj);
end

function check_shape(values,expected,path,name)
observed=size(values); observed(end+1:numel(expected))=1;
expected(end+1:numel(observed))=1;
if ~isnumeric(values)||~isreal(values)||~isequal(observed,expected)
    error('dsc:DiagnosticChunk','Invalid %s draw dimensions in %s.',name,path);
end
end

function details=empty_details()
details=table('Size',[0 13],'VariableTypes', ...
    {'string','string','double','double','string','double','double','double', ...
    'double','double','logical','logical','string'},'VariableNames', ...
    {'quantity','family','index_i','index_j','date','rhat','ess_bulk','ess_tail', ...
    'mcse_mean','mcse_sd_ratio','available','passed','reason'});
end

function details=quantity_details(meta)
count=meta.T*(meta.m+meta.nr);
quantity=strings(count,1); family=strings(count,1); date=strings(count,1);
index_i=nan(count,1); index_j=nan(count,1);
row=0;
date_labels=string(meta.dates,'yyyy-MM-dd');
for name={'h','P_pairs'}
    field=name{1}; width=meta.m; if strcmp(field,'P_pairs'), width=meta.nr; end
    for k=1:width
        ix=row+(1:meta.T); family(ix)=string(field); date(ix)=date_labels;
        if strcmp(field,'P_pairs')
            index_i(ix)=meta.pair_i(k); index_j(ix)=meta.pair_j(k);
            prefix=sprintf('P(%d,%d)[',meta.pair_i(k),meta.pair_j(k));
        else
            index_i(ix)=k; prefix=sprintf('%s(%d)[',field,k);
        end
        quantity(ix)=string(prefix)+date_labels+"]"; row=row+meta.T;
    end
end
details=table(quantity,family,index_i,index_j,date,nan(count,1),nan(count,1), ...
    nan(count,1),nan(count,1),nan(count,1),false(count,1),false(count,1),strings(count,1), ...
    'VariableNames',{'quantity','family','index_i','index_j','date','rhat','ess_bulk', ...
    'ess_tail','mcse_mean','mcse_sd_ratio','available','passed','reason'});
end

function value=finite_extreme(values,which)
values=values(isfinite(values)); value=NaN;
if isempty(values), return; end
if strcmp(which,'max'), value=max(values); else, value=min(values); end
end

function write_outputs(report,details)
temporary=[tempname(fileparts(report.paths.mat)) '.mat'];
save(temporary,'report','details','-v7.3'); movefile(temporary,report.paths.mat,'f');
temporary=[tempname(fileparts(report.paths.csv)) '.csv'];
writetable(details,temporary); movefile(temporary,report.paths.csv,'f');
temporary=[tempname(fileparts(report.paths.json)) '.json'];
fid=fopen(temporary,'w');
if fid<0, error('dsc:DiagnosticOutput','Cannot write diagnostic JSON.'); end
closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true)); clear closer;
movefile(temporary,report.paths.json,'f');
end

function remove_staging(folder)
if isfolder(folder), rmdir(folder,'s'); end
end

function figure_paths = plot_bayes_correlation_paths(run_dir)
%PLOT_BAYES_CORRELATION_PATHS Plot short-run Bayesian P(i,j,t) paths.
% Pools retained post-warm-up draws while preserving original chain identity.
% Uses direct chunks or the coordinator summary's chain_dirs. The colored
% vertical bands use Figure 1's continuous red-white-blue correlation scale.
if nargin<1 || strlength(string(run_dir))==0
    error('DSC:RunDir','Provide a sampler run directory.');
end
run_dir = char(run_dir);
summary_path = fullfile(run_dir,'pilot_summary.json');
if exist(summary_path,'file')
    run_summary = jsondecode(fileread(summary_path));
else
    run_summary=struct();
end
[P,identity,iterations,draw_chain_ids,chain_dirs,chain_ids,retained_draws_per_chain,inference]= ...
    read_chain_draws(run_dir,run_summary);
if size(P,3)<2, error('DSC:TooFewDraws','At least two retained draws are needed.'); end
dates=identity.dates; tickers=string(identity.tickers);
pair_i=identity.pair_i; pair_j=identity.pair_j;
warmup=max(0,min(iterations)-1);
if isfield(run_summary,'burnin'), warmup=run_summary.burnin; end
if ~isnumeric(warmup)||~isscalar(warmup)||~isfinite(warmup)||warmup<0||warmup~=floor(warmup)|| ...
        any(iterations<=warmup)|| ...
        (isfield(run_summary,'saved_draws')&&run_summary.saved_draws~=size(P,3))|| ...
        (isfield(run_summary,'saved_draws_per_chain')&& ...
        ~isequal(run_summary.saved_draws_per_chain(:),retained_draws_per_chain(:)))|| ...
        (isfield(run_summary,'chain_ids')&&~isequal(run_summary.chain_ids(:),chain_ids(:)))|| ...
        (isfield(run_summary,'num_chains')&&run_summary.num_chains~=numel(chain_ids))
    error('DSC:PathDraws','Retained iterations or draw count disagree with the sampler summary.');
end
warmup_per_chain=repmat(warmup,1,numel(chain_ids));
convergence_established=false;
if isfield(run_summary,'convergence_established')
    flag=run_summary.convergence_established;
    if ~isscalar(flag)||~(islogical(flag)||isnumeric(flag))||~ismember(flag,[0 1])
        error('DSC:ConvergenceMetadata','Invalid convergence-established flag in the run summary.');
    end
    convergence_established=logical(run_summary.convergence_established);
end
validate_plot_convergence(run_dir,run_summary,chain_ids,retained_draws_per_chain, ...
    numel(dates),numel(tickers),convergence_established);
q=prctile(P,[16 50 84],3);
lower=q(:,:,1); median_path=q(:,:,2); upper=q(:,:,3);
prior_sampling=[]; draw_parameter_ids=[];
if ~isempty(inference)&&strcmp(inference.inference_mode,'prior_parameter_smoothing')
    prior_sampling=run_summary.prior_sampling;
    for c=1:numel(chain_ids)
        draw_parameter_ids=[draw_parameter_ids repmat( ...
            prior_sampling.parameter_draw_ids_per_chain(c),1,retained_draws_per_chain(c))]; %#ok<AGROW>
    end
end
draw_note=sprintf('%d pooled retained draws from %d chain(s) | %d warm-up iterations removed per chain', ...
    size(P,3),numel(chain_ids),warmup);

fig_dir = fullfile(run_dir,'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
save(fullfile(run_dir,'posterior_correlation_bands.mat'), ...
    'lower','median_path','upper','dates','tickers','pair_i','pair_j','iterations','inference', ...
    'draw_chain_ids','chain_ids','chain_dirs','retained_draws_per_chain','warmup_per_chain', ...
    'prior_sampling','draw_parameter_ids','-v7.3');

pair_labels = erase(tickers(pair_i)," Index") + " / " + erase(tickers(pair_j)," Index");
pages = ceil(numel(pair_labels)/9);
figure_paths = strings(pages,1);

for page = 1:pages
    ix = (page-1)*9 + (1:9);
    ix = ix(ix<=numel(pair_labels));
    figure_paths(page) = plot_page(fullfile(fig_dir,sprintf('bayes_correlation_paths_all_%02d.pdf',page)), ...
        dates, median_path, lower, upper, pair_labels, ix, ...
        sprintf('All Bayesian correlation paths, page %d of %d',page,pages),inference,convergence_established,draw_note,false);
end
end

function validate_plot_convergence(run_dir,summary,ids,counts,T,m,established)
if ~isfield(summary,'convergence')
    if established
        error('DSC:ConvergenceMetadata','A converged plot requires a saved diagnostic report.');
    end
    return
end
report=summary.convergence;
if ~isstruct(report)||~isscalar(report)||~isfield(report,'passed')|| ...
        ~isequal(established,report.passed)
    error('DSC:ConvergenceMetadata','The plot convergence flag disagrees with its diagnostic report.');
end
if isfield(summary,'inference_mode')&&strcmp(summary.inference_mode,'prior_parameter_smoothing')&& ...
        (established||~isfield(report,'status')||~strcmp(report.status,'not_checked'))
    error('DSC:ConvergenceMetadata','Prior-mixture diagnostics must be reported separately within each fixed parameter group.');
end
nr=m*(m-1)/2; expected=T*(m+nr);
if established&&(~isfield(summary,'inference_mode')|| ...
        ~ismember(summary.inference_mode,{'parameter_estimation','fixed_parameter_smoothing'}))
    error('DSC:ConvergenceMetadata','A passing diagnostic report must identify its inference target.');
end
dsc_validate_convergence(report,ids,counts,expected);
if established
    path=fullfile(run_dir,'convergence','convergence_diagnostics.json');
    if ~isfile(path)||~isequaln(report,jsondecode(fileread(path)))
        error('DSC:ConvergenceMetadata','Passing run-summary diagnostics disagree with the saved diagnostic JSON.');
    end
end
end

function [P,identity,iterations,draw_ids,dirs,ids,counts,inference]=read_chain_draws(run_dir,summary)
dirs={run_dir};
if isfield(summary,'chain_dirs')
    dirs=summary.chain_dirs;
    if ischar(dirs)||isstring(dirs), dirs=cellstr(dirs); end
end
if ~iscell(dirs)||isempty(dirs)||~all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),dirs))
    error('DSC:PathIdentity','The sampler summary must contain nonempty chain directory paths.');
end
dirs=dirs(:)'; canonical=cell(size(dirs));
for c=1:numel(dirs)
    file=java.io.File(dirs{c}); canonical{c}=char(file.getCanonicalPath());
end
if numel(unique(canonical))~=numel(dirs)
    error('DSC:PathIdentity','Each posterior chain directory must be distinct.');
end
blocks={}; identity=[]; iterations=[]; draw_ids=[]; ids=nan(1,numel(dirs)); counts=zeros(1,numel(dirs));
for c=1:numel(dirs)
    files=dir(fullfile(dirs{c},'posterior_chunk_*.mat'));
    if isempty(files), error('DSC:NoChunks','No posterior chunks found in %s',dirs{c}); end
    numbers=zeros(1,numel(files));
    for n=1:numel(files)
        token=regexp(files(n).name,'^posterior_chunk_(\d+)\.mat$','tokens','once');
        if isempty(token), error('DSC:PathDraws','Invalid posterior chunk filename.'); end
        numbers(n)=str2double(token{1});
    end
    if numel(unique(numbers))~=numel(numbers), error('DSC:PathDraws','Duplicate chunk indices.'); end
    [~,order]=sort(numbers); files=files(order); last_iteration=0;
    for n=1:numel(files)
        loaded=load(fullfile(dirs{c},files(n).name),'chunk');
        if ~isfield(loaded,'chunk'), error('DSC:PathIdentity','Missing posterior chunk.'); end
        chunk=loaded.chunk;
        fields={'dates','tickers','pair_i','pair_j','chain_id','iterations','P_pairs'};
        if ~isstruct(chunk)||~isscalar(chunk)||~all(isfield(chunk,fields))
            error('DSC:PathIdentity','A posterior chunk is missing required identity or draw fields.');
        end
        current=struct('dates',chunk.dates,'tickers',chunk.tickers, ...
            'pair_i',chunk.pair_i,'pair_j',chunk.pair_j);
        if isempty(identity), identity=current;
        elseif ~isequaln(identity,current)
            error('DSC:PathIdentity','Posterior chunks have different dates, series or pair order.');
        end
        if ~isnumeric(chunk.chain_id)||~isscalar(chunk.chain_id)||~isfinite(chunk.chain_id)|| ...
                chunk.chain_id<1||chunk.chain_id~=floor(chunk.chain_id)
            error('DSC:PathIdentity','Invalid original chain identity.');
        end
        if n==1, ids(c)=chunk.chain_id;
        elseif ids(c)~=chunk.chain_id
            error('DSC:PathIdentity','One directory contains multiple original chains.');
        end
        index=chunk.iterations(:)'; draws=numel(index);
        if ~isnumeric(index)||~isreal(index)||isempty(index)||any(~isfinite(index))|| ...
                any(index~=floor(index))||index(1)<=last_iteration||any(diff(index)<=0)
            error('DSC:PathDraws','Retained iterations must increase strictly within each original chain.');
        end
        T=numel(chunk.dates); m=numel(chunk.tickers); nr=m*(m-1)/2;
        [ii,jj]=find(tril(true(m),-1)); values=chunk.P_pairs;
        if ~isdatetime(chunk.dates)||~iscolumn(chunk.dates)||T<2||any(isnat(chunk.dates))|| ...
                any(diff(chunk.dates)<=seconds(0))||numel(unique(string(chunk.tickers)))~=m|| ...
                ~isequal(chunk.pair_i(:),ii)||~isequal(chunk.pair_j(:),jj)|| ...
                ~isnumeric(values)||~isreal(values)||size(values,1)~=T||size(values,2)~=nr|| ...
                size(values,3)~=draws||numel(values)~=T*nr*draws|| ...
                any(~isfinite(values(:)))||any(abs(values(:))>1)
            error('DSC:PathIdentity','Invalid dates, pairs or finite correlation draw dimensions.');
        end
        blocks{end+1}=values; iterations=[iterations index]; %#ok<AGROW>
        draw_ids=[draw_ids repmat(ids(c),1,draws)]; %#ok<AGROW>
        counts(c)=counts(c)+draws; last_iteration=index(end);
    end
end
if numel(unique(ids))~=numel(ids)
    error('DSC:PathIdentity','Posterior directories must have distinct original chain IDs.');
end
inference=read_inference(run_dir,identity.dates,summary);
if ~isempty(inference)&& ...
        ((isfield(inference,'num_chains')&&inference.num_chains~=numel(ids))|| ...
        (isfield(inference,'retained_draws_per_chain')&& ...
        ~isequal(inference.retained_draws_per_chain(:),counts(:))))
    error('DSC:InferenceMetadata','Inference chain counts disagree with the saved draws.');
end
if ~isempty(inference)&&strcmp(inference.inference_mode,'prior_parameter_smoothing')
    validate_prior_mixture(inference,summary,ids,counts,numel(identity.dates),numel(identity.tickers));
end
% Only a complete, equally weighted prior mixture may enter pooled bands.
P=cat(3,blocks{:});
end

function path = plot_page(path,dates,median_path,lower,upper,pair_labels,indices,title_text,inference,convergence_established,draw_note,show_header)
fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1200 1500]);
layout = tiledlayout(numel(indices),1,'TileSpacing','compact','Padding','compact');
if show_header && ~isempty(inference), layout.Padding='loose'; end
x = datenum(dates);
cmap = red_blue_map(257);
colormap(fig,cmap);
for pos = 1:numel(indices)
    k = indices(pos);
    ax = nexttile; hold(ax,'on');
    draw_bars(ax,x,median_path(:,k));
    fill(ax,[x; flipud(x)],[lower(:,k); flipud(upper(:,k))],[0.08 0.19 0.31], ...
        'FaceAlpha',0.16,'EdgeColor','none');
    plot(ax,x,median_path(:,k),'Color',[0.02 0.10 0.20],'LineWidth',1.1);
    yline(ax,0,'Color',[0.55 0.55 0.55],'LineWidth',0.5);
    cutoff=NaN;
    if ~isempty(inference)&&~strcmp(inference.inference_mode,'prior_parameter_smoothing')
        cutoff=datenum(datetime(inference.parameter_estimation_end,'InputFormat','yyyy-MM-dd'));
        if cutoff>x(1)&&cutoff<x(end)
            xline(ax,cutoff,'--','Color',[.25 .25 .25],'LineWidth',.9,'HandleVisibility','off');
        end
    end
    ylim(ax,[-1 1]); xlim(ax,[x(1) x(end)]);
    clim(ax,[-1 1]);
    datetick(ax,'x','yyyy','keeplimits');
    ylabel(ax,'P_{ij,t}');
    title(ax,char(pair_labels(k)),'Interpreter','none','FontWeight','normal');
    box(ax,'off'); grid(ax,'on');
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = 'Median correlation';
cb.Ticks = -1:0.25:1;
if show_header
    if convergence_established, status='Convergence diagnostics passed';
    else, status='Convergence not established'; end
    if isempty(inference)
        header={title_text,'Median and 68% pointwise band | Vertical colors show the median correlation value', ...
            draw_note,status};
    else
        if strcmp(inference.inference_mode,'prior_parameter_smoothing')
            mode='Prior-parameter conditional smoothing | Prior-weighted parameter uncertainty';
        elseif strcmp(inference.inference_mode,'fixed_parameter_smoothing')
            if isfield(inference,'parameter_smoothing')&&strcmp(inference.parameter_smoothing,'draws')
                mode='Conditional smoothing with saved parameter draws | Bands include saved-parameter variation';
            else
                mode='Conditional on posterior-mean parameters | Bands exclude parameter uncertainty';
            end
        else
            mode='Joint parameter and state draws | Bands include parameter uncertainty';
        end
        if strcmp(inference.inference_mode,'prior_parameter_smoothing')
            sample=sprintf('Priors calibrated: %s to %s | Smoothing: %s to %s', ...
                inference.prior_calibration_start,inference.prior_calibration_end, ...
                inference.smoothing_start,inference.smoothing_end);
            estimated=sprintf('%d equally weighted prior parameter draws | Conditional-state diagnostics recorded by parameter draw', ...
                inference.prior_parameter_draws);
        else
            sample=sprintf('Parameters: %s to %s | Smoothing through %s', ...
                inference.parameter_estimation_start,inference.parameter_estimation_end,inference.smoothing_end);
            estimated=['Model estimated at ' inference.parameters_estimated_at ' (UTC)'];
            if cutoff>x(1)&&cutoff<x(end), estimated=[estimated ' | Dashed line: estimation cutoff']; end
        end
        if isfield(inference,'calibration_weeks')
            sample=sprintf('%s | First %d weeks excluded',sample,inference.calibration_weeks);
        end
        header={title_text,mode,sample,estimated, ...
            'Median and 68% pointwise band | Figure 1 colors | Historical states may revise',draw_note,status};
    end
    title(layout,header,'FontWeight','bold','Interpreter','none','FontSize',11);
end
if isempty(inference)
    exportgraphics(fig,path,'ContentType','vector');
else
    % Fixed paper margins keep long provenance headers and axis labels clear
    % of the PDF edge; exportgraphics tightly crops figure text in R2024b.
    set(fig,'PaperUnits','inches','PaperSize',[8.5 10.5], ...
        'PaperPosition',[.35 .35 7.8 9.8],'PaperPositionMode','manual');
    print(fig,path,'-dpdf','-painters');
end
close(fig);
end

function inference=read_inference(run_dir,dates,summary)
inference=[];
metadata_path=fullfile(run_dir,'inference_metadata.json');
manifest_path=fullfile(run_dir,'run_manifest.json');
if isfile(manifest_path), manifest=jsondecode(fileread(manifest_path));
else, manifest=struct(); end
if ~isfile(metadata_path)
    if (isfield(manifest,'inference')&&isstruct(manifest.inference))|| ...
            (isfield(summary,'inference_mode')&&ismember(summary.inference_mode, ...
            {'fixed_parameter_smoothing','prior_parameter_smoothing'}))
        error('DSC:InferenceMetadata','Inference metadata is missing; cannot label conditional correlation paths.');
    end
    return
end
inference=jsondecode(fileread(metadata_path));
if isstruct(inference)&&isscalar(inference)&&isfield(inference,'inference_mode')&& ...
        strcmp(inference.inference_mode,'prior_parameter_smoothing')
    if ~isfield(manifest,'inference')||~isequaln(inference,manifest.inference)
        error('DSC:InferenceMetadata','Prior inference metadata and run manifest disagree.');
    end
    validate_prior_inference(inference,dates,summary);
    return
end
required={'inference_mode','parameter_file','parameters_estimated_at','parameter_estimation_start', ...
    'parameter_estimation_end','parameter_estimation_draws','parameter_estimation_run', ...
    'smoothing_end','parameter_uncertainty_in_bands','parameter_smoothing'};
if ~isstruct(inference)||~isscalar(inference)||~all(isfield(inference,required))|| ...
        ~isfield(manifest,'inference')||~isequaln(inference,manifest.inference)
    error('DSC:InferenceMetadata','Inference metadata and run manifest are missing fields or disagree.');
end
for k=[1 2 3 4 5 7 8]
    value=inference.(required{k});
    if ~(ischar(value)&&isrow(value)&&~isempty(value))
        error('DSC:InferenceMetadata','Inference date, mode and source fields must be nonempty text.');
    end
end
mode=inference.inference_mode;
parameter_smoothing=inference.parameter_smoothing;
if ~ismember(mode,{'parameter_estimation','fixed_parameter_smoothing'})|| ...
        ~ismember(parameter_smoothing,{'draws','mean'})|| ...
        ~islogical(inference.parameter_uncertainty_in_bands)||~isscalar(inference.parameter_uncertainty_in_bands)|| ...
        inference.parameter_uncertainty_in_bands~=(strcmp(mode,'parameter_estimation')||strcmp(parameter_smoothing,'draws'))|| ...
        ~isscalar(inference.parameter_estimation_draws)||inference.parameter_estimation_draws<2|| ...
        ~isfinite(inference.parameter_estimation_draws)|| ...
        inference.parameter_estimation_draws~=floor(inference.parameter_estimation_draws)
    error('DSC:InferenceMetadata','Inference mode, retained parameter count or uncertainty scope is invalid.');
end
try
    start=datetime(inference.parameter_estimation_start,'InputFormat','yyyy-MM-dd');
    cutoff=datetime(inference.parameter_estimation_end,'InputFormat','yyyy-MM-dd');
    last=datetime(inference.smoothing_end,'InputFormat','yyyy-MM-dd');
    estimated=datetime(inference.parameters_estimated_at,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC');
catch
    error('DSC:InferenceMetadata','Invalid dates in inference metadata.');
end
if any(isnat([start cutoff last]))||isnat(estimated)||start~=dates(1)||last~=dates(end)|| ...
        cutoff<start||cutoff>last||~any(dates==cutoff)|| ...
        (strcmp(mode,'parameter_estimation')&&cutoff~=last)|| ...
        (isfield(summary,'inference_mode')&&~strcmp(summary.inference_mode,mode))|| ...
        (isfield(summary,'first_date')&&~strcmp(summary.first_date,char(string(dates(1),'yyyy-MM-dd'))))|| ...
        (isfield(summary,'last_date')&&~strcmp(summary.last_date,char(string(dates(end),'yyyy-MM-dd'))))
    error('DSC:InferenceMetadata','Inference dates or mode disagree with the saved correlation sample.');
end
if isfield(inference,'calibration_weeks')
    n=inference.calibration_weeks;
    if ~isscalar(n)||~isfinite(n)||n<1||n~=floor(n)|| ...
            ~all(isfield(inference,{'calibration_start','calibration_end','smoothing_start','full_data_start'}))|| ...
            ~strcmp(inference.smoothing_start,inference.parameter_estimation_start)|| ...
            ~strcmp(inference.full_data_start,inference.calibration_start)|| ...
            ~isfield(summary,'excluded_initial_weeks')||summary.excluded_initial_weeks~=n
        error('DSC:InferenceMetadata','Calibration exclusion metadata is inconsistent.');
    end
    initial=datetime(inference.calibration_start,'InputFormat','yyyy-MM-dd');
    final=datetime(inference.calibration_end,'InputFormat','yyyy-MM-dd');
    if isnat(initial)||isnat(final)||start~=final+calweeks(1)||final~=initial+calweeks(n-1)
        error('DSC:InferenceMetadata','Calibration and model dates must be adjacent, non-overlapping weekly blocks.');
    end
end
end

function validate_prior_inference(inference,dates,summary)
required={'parameter_source','parameter_file','prior_parameter_file','parameter_smoothing', ...
    'parameter_uncertainty_in_bands','uncertainty_source','prior_parameter_draws', ...
    'prior_calibration_start','prior_calibration_end','calibration_weeks', ...
    'calibration_start','calibration_end','smoothing_start','smoothing_end', ...
    'num_chains','retained_draws_per_chain'};
if ~all(isfield(inference,required))||~strcmp(inference.parameter_source,'prior')|| ...
        ~strcmp(inference.parameter_smoothing,'prior_draws')|| ...
        ~isequal(inference.parameter_uncertainty_in_bands,true)|| ...
        ~strcmp(inference.uncertainty_source,'prior_parameter_and_conditional_state')|| ...
        ~positive_integer(inference.prior_parameter_draws)|| ...
        ~positive_integer(inference.calibration_weeks)||~positive_integer(inference.num_chains)|| ...
        ~ischar(inference.parameter_file)||~isrow(inference.parameter_file)||isempty(inference.parameter_file)|| ...
        ~strcmp(inference.parameter_file,inference.prior_parameter_file)
    error('DSC:InferenceMetadata','Prior smoothing provenance or uncertainty scope is invalid.');
end
% Prior draws have no fitted-parameter provenance, even when legacy keys are saved.
empty_fields={'parameters_estimated_at','parameter_estimation_start','parameter_estimation_end','parameter_estimation_run'};
for k=1:numel(empty_fields)
    if isfield(inference,empty_fields{k})&&~isempty(inference.(empty_fields{k}))
        error('DSC:InferenceMetadata','Prior parameter draws cannot claim parameter estimation provenance.');
    end
end
if isfield(inference,'parameter_estimation_draws')&&~isequal(inference.parameter_estimation_draws,0)
    error('DSC:InferenceMetadata','Prior parameter draws cannot be counted as estimated posterior draws.');
end
date_fields={'prior_calibration_start','prior_calibration_end','calibration_start', ...
    'calibration_end','smoothing_start','smoothing_end'};
parsed=NaT(1,numel(date_fields));
for k=1:numel(date_fields)
    value=inference.(date_fields{k});
    if ~ischar(value)||~isrow(value)||isempty(regexp(value,'^\d{4}-\d{2}-\d{2}$','once'))
        error('DSC:InferenceMetadata','Prior calibration and smoothing dates must be ISO dates.');
    end
    try
        parsed(k)=datetime(value,'InputFormat','yyyy-MM-dd');
    catch
        error('DSC:InferenceMetadata','Invalid prior calibration or smoothing dates.');
    end
    if isnat(parsed(k))||~strcmp(char(string(parsed(k),'yyyy-MM-dd')),value)
        error('DSC:InferenceMetadata','Invalid prior calibration or smoothing dates.');
    end
end
if parsed(1)>parsed(2)||parsed(2)>dates(end)||parsed(1)>parsed(3)||parsed(2)<parsed(4)|| ...
        parsed(5)~=dates(1)||parsed(6)~=dates(end)|| ...
        parsed(4)~=parsed(3)+calweeks(inference.calibration_weeks-1)|| ...
        parsed(5)~=parsed(4)+calweeks(1)|| ...
        ~isfield(summary,'excluded_initial_weeks')||summary.excluded_initial_weeks~=inference.calibration_weeks|| ...
        ~isfield(summary,'inference_mode')||~strcmp(summary.inference_mode,inference.inference_mode)|| ...
        (isfield(summary,'first_date')&&~strcmp(summary.first_date,inference.smoothing_start))|| ...
        (isfield(summary,'last_date')&&~strcmp(summary.last_date,inference.smoothing_end))
    error('DSC:InferenceMetadata','Prior calibration and smoothing provenance disagree with saved correlation dates.');
end
if isfield(inference,'full_data_start')&&~strcmp(inference.full_data_start,inference.calibration_start)
    error('DSC:InferenceMetadata','The full sample must begin with the initial prior calibration block.');
end
end

function validate_prior_mixture(inference,summary,ids,counts,T,m)
if ~isfield(summary,'prior_sampling')||~isstruct(summary.prior_sampling)||~isscalar(summary.prior_sampling)
    error('DSC:InferenceMetadata','Prior mixture provenance is missing.');
end
p=summary.prior_sampling;
required={'requested_parameter_draws','conditional_chains_per_parameter','parameter_draw_ids_per_chain', ...
    'mixture_complete','mixture_weights','requested_retained_draws_per_chain','conditional_convergence'};
if ~all(isfield(p,required))||~positive_integer(p.requested_parameter_draws)|| ...
        ~positive_integer(p.conditional_chains_per_parameter)||~positive_integer(p.requested_retained_draws_per_chain)|| ...
        ~isequal(p.mixture_complete,true)
    error('DSC:InferenceMetadata','Prior mixture is incomplete or its requested quotas are invalid.');
end
K=p.requested_parameter_draws; C=p.conditional_chains_per_parameter;
weights=p.mixture_weights; group_ids=p.parameter_draw_ids_per_chain;
if K~=inference.prior_parameter_draws||numel(ids)~=K*C|| ...
        ~isnumeric(group_ids)||~isequal(group_ids(:)',repelem(1:K,C))|| ...
        ~isnumeric(weights)||~isreal(weights)||numel(weights)~=K||any(~isfinite(weights(:)))|| ...
        any(abs(weights(:)-1/K)>8*eps(1/K))|| ...
        any(counts~=p.requested_retained_draws_per_chain)|| ...
        ~isfield(summary,'num_chains')||summary.num_chains~=K*C|| ...
        ~isfield(summary,'chain_ids')||~isequal(summary.chain_ids(:),ids(:))|| ...
        ~isfield(summary,'saved_draws_per_chain')||~isequal(summary.saved_draws_per_chain(:),counts(:))|| ...
        numel(p.conditional_convergence)~=K|| ...
        ~(iscell(p.conditional_convergence)||isstruct(p.conditional_convergence))
    error('DSC:InferenceMetadata','Prior bands require every parameter group, equal weights and equal completed chain quotas.');
end
if ~isfield(summary,'convergence')||~isstruct(summary.convergence)|| ...
        ~all(isfield(summary.convergence,{'status','passed'}))|| ...
        ~strcmp(summary.convergence.status,'not_checked')||~isequal(summary.convergence.passed,false)
    error('DSC:ConvergenceMetadata','Prior mixtures require within-parameter diagnostics, not a pooled convergence claim.');
end
for k=1:K
    if iscell(p.conditional_convergence), report=p.conditional_convergence{k};
    else, report=p.conditional_convergence(k); end
    selection=(k-1)*C+(1:C);
    if isstruct(report)&&isfield(report,'chain_ids')&& ...
            ~isequal(report.chain_ids(:),reshape(ids(selection),[],1))
        error('DSC:ConvergenceMetadata','Conditional diagnostics must retain their original fixed-parameter chain identities.');
    end
    dsc_validate_convergence(report,ids(selection),counts(selection),T*(m+m*(m-1)/2));
end
end

function yes=positive_integer(value)
yes=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>=1&&value==floor(value);
end

function draw_bars(ax,x,y)
edges = [x(1); (x(1:end-1)+x(2:end))/2; x(end)];
image(ax,'XData',[edges(1) edges(end)],'YData',[-1 1], ...
    'CData',reshape(y,1,[]),'CDataMapping','scaled','AlphaData',0.22);
set(ax,'YDir','normal');
end

function cmap = red_blue_map(n)
if nargin<1, n = 257; end
% ColorBrewer/Matplotlib RdBu_r anchors, matching Figure 1.
anchors = [5 48 97; 33 102 172; 67 147 195; 146 197 222; 209 229 240; ...
    247 247 247; 253 219 199; 244 165 130; 214 96 77; 178 24 43; 103 0 31]/255;
cmap = interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n),'linear');
end

function write_json_local(path,value,raw_inference)
fid = fopen(path,'w');
assert(fid>=0,'DSC:WriteFailed','Cannot open %s',path);
cleanup = onCleanup(@()fclose(fid));
payload=jsonencode(value,PrettyPrint=true);
if nargin>=3&&~isempty(raw_inference)
    payload=[payload(1:end-1) ',"inference":' raw_inference '}'];
end
fprintf(fid,'%s\n',payload);
end

function value = boolean_text(value)
if value, value='true'; else, value='false'; end
end

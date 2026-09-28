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
[P,identity,iterations,draw_chain_ids,chain_dirs,chain_ids,retained_draws_per_chain]= ...
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
q=prctile(P,[16 50 84],3);
lower=q(:,:,1); median_path=q(:,:,2); upper=q(:,:,3);
inference=read_inference(run_dir,dates,run_summary);
if ~isempty(inference)&& ...
        ((isfield(inference,'num_chains')&&inference.num_chains~=numel(chain_ids))|| ...
        (isfield(inference,'retained_draws_per_chain')&& ...
        ~isequal(inference.retained_draws_per_chain(:),retained_draws_per_chain(:))))
    error('DSC:InferenceMetadata','Inference chain counts disagree with the saved posterior draws.');
end
validate_plot_convergence(run_dir,run_summary,chain_ids,retained_draws_per_chain, ...
    numel(dates),numel(tickers),convergence_established);
draw_note=sprintf('%d pooled retained draws from %d chain(s) | %d warm-up iterations removed per chain', ...
    size(P,3),numel(chain_ids),warmup);

fig_dir = fullfile(run_dir,'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
save(fullfile(run_dir,'posterior_correlation_bands.mat'), ...
    'lower','median_path','upper','dates','tickers','pair_i','pair_j','iterations','inference', ...
    'draw_chain_ids','chain_ids','chain_dirs','retained_draws_per_chain','warmup_per_chain','-v7.3');

pair_labels = erase(tickers(pair_i)," Index") + " / " + erase(tickers(pair_j)," Index");
selected = select_pairs(median_path,pair_labels);
pages = ceil(numel(pair_labels)/9);
figure_paths = strings(pages+1,1);
figure_paths(1) = plot_page(fullfile(fig_dir,'bayes_correlation_paths_selected.pdf'), ...
    dates, median_path, lower, upper, pair_labels, selected, ...
    'Selected Bayesian correlation paths',inference,convergence_established,draw_note,true);

for page = 1:pages
    ix = (page-1)*9 + (1:9);
    ix = ix(ix<=numel(pair_labels));
    figure_paths(page+1) = plot_page(fullfile(fig_dir,sprintf('bayes_correlation_paths_all_%02d.pdf',page)), ...
        dates, median_path, lower, upper, pair_labels, ix, ...
        sprintf('All Bayesian correlation paths, page %d of %d',page,pages),inference,convergence_established,draw_note,false);
end
manifest=struct('run_dir',run_dir,'retained_draws',size(P,3),'warmup_removed',warmup, ...
    'num_chains',numel(chain_ids),'chain_ids',chain_ids,'chain_dirs',{chain_dirs}, ...
    'retained_draws_per_chain',retained_draws_per_chain,'warmup_removed_per_chain',warmup_per_chain, ...
    'draw_chain_ids',draw_chain_ids,'retained_iterations',iterations, ...
    'paths',{cellstr(figure_paths)},'color_scale',[-1 1], ...
    'first_date',char(string(dates(1),'yyyy-MM-dd')),'last_date',char(string(dates(end),'yyyy-MM-dd')), ...
    'color_rule','Continuous RdBu_r scale applied to the posterior median correlation.', ...
    'created_at',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')), ...
    'note',sprintf('%d pooled retained draws across %d chains after %d warm-up iterations per chain; pointwise 16/50/84 percentiles; convergence established: %s.', ...
    size(P,3),numel(chain_ids),warmup,boolean_text(convergence_established)));
raw_inference='';
if ~isempty(inference)
    % jsondecode converts JSON null to [] and singleton object arrays to
    % structs. Preserve the verified source JSON when copying provenance.
    raw_inference=fileread(fullfile(run_dir,'inference_metadata.json'));
end
write_json_local(fullfile(fig_dir,'bayes_correlation_paths_manifest.json'),manifest,raw_inference);
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
nr=m*(m-1)/2; expected=T*(2*m+nr);
if isfield(summary,'inference_mode')&&strcmp(summary.inference_mode,'parameter_estimation')
    expected=expected+m*(m+1)/2+m+nr;
elseif established&&(~isfield(summary,'inference_mode')|| ...
        ~strcmp(summary.inference_mode,'fixed_parameter_smoothing'))
    error('DSC:ConvergenceMetadata','A passing diagnostic report must identify its inference target.');
end
dsc_validate_convergence(report,ids,counts,expected);
if established
    path=fullfile(run_dir,'convergence_diagnostics.json');
    if ~isfile(path)||~isequaln(report,jsondecode(fileread(path)))
        error('DSC:ConvergenceMetadata','Passing run-summary diagnostics disagree with the saved diagnostic JSON.');
    end
end
end

function [P,identity,iterations,draw_ids,dirs,ids,counts]=read_chain_draws(run_dir,summary)
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
P=[]; identity=[]; iterations=[]; draw_ids=[]; ids=nan(1,numel(dirs)); counts=zeros(1,numel(dirs));
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
        P=cat(3,P,values); iterations=[iterations index]; %#ok<AGROW>
        draw_ids=[draw_ids repmat(ids(c),1,draws)]; %#ok<AGROW>
        counts(c)=counts(c)+draws; last_iteration=index(end);
    end
end
if numel(unique(ids))~=numel(ids)
    error('DSC:PathIdentity','Posterior directories must have distinct original chain IDs.');
end
end

function selected = select_pairs(median_path,pair_labels)
range = max(median_path,[],1) - min(median_path,[],1);
score = max(abs(median_path),[],1) + range;
[~,order] = sort(score,'descend');
must = find(pair_labels=="TUKXG / SX5T" | pair_labels=="SX5T / BCOMCOT");
selected = unique([must(:)' order(1:min(6,numel(order)))],'stable');
selected = selected(1:min(6,numel(selected)));
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
    if ~isempty(inference)
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
        if strcmp(inference.inference_mode,'fixed_parameter_smoothing')
            if isfield(inference,'parameter_smoothing')&&strcmp(inference.parameter_smoothing,'draws')
                mode='Conditional smoothing with saved parameter draws | Bands include saved-parameter variation';
            else
                mode='Conditional on posterior-mean parameters | Bands exclude parameter uncertainty';
            end
        else
            mode='Joint parameter and state draws | Bands include parameter uncertainty';
        end
        sample=sprintf('Parameters: %s to %s | Smoothing through %s', ...
            inference.parameter_estimation_start,inference.parameter_estimation_end,inference.smoothing_end);
        if isfield(inference,'calibration_weeks')
            sample=sprintf('%s | First %d weeks excluded',sample,inference.calibration_weeks);
        end
        estimated=['Model estimated at ' inference.parameters_estimated_at ' (UTC)'];
        if cutoff>x(1)&&cutoff<x(end), estimated=[estimated ' | Dashed line: estimation cutoff']; end
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
            (isfield(summary,'inference_mode')&&strcmp(summary.inference_mode,'fixed_parameter_smoothing'))
        error('DSC:InferenceMetadata','Inference metadata is missing; cannot label conditional correlation paths.');
    end
    return
end
inference=jsondecode(fileread(metadata_path));
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

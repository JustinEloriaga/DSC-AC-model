function [figure_paths,plot_summary] = plot_bayes_correlation_paths(run_dir,file_prefix)
%PLOT_BAYES_CORRELATION_PATHS Plot pooled and chain-specific P(i,j,t) means.
% Pools retained post-warm-up draws while preserving original chain identity.
% Uses direct chunks, prefixed flat chunks or the coordinator chain_dirs. The colored
% vertical bands use Figure 1's continuous red-white-blue correlation scale.
if nargin<1 || strlength(string(run_dir))==0
    error('DSC:RunDir','Provide a sampler run directory.');
end
run_dir = char(run_dir);
if ~isfolder(run_dir), error('DSC:RunDir','Provide an existing sampler run directory.'); end
if nargin<2, file_prefix=''; end
file_prefix=validate_file_prefix(file_prefix);
manifest=struct();
manifest_path=fullfile(run_dir,[file_prefix 'run_manifest.json']);
if isfile(manifest_path), manifest=jsondecode(fileread(manifest_path)); end
parameters=read_parameter_report(run_dir,file_prefix,manifest);
summary_path = fullfile(run_dir,[file_prefix 'pilot_summary.json']);
result_path=fullfile(run_dir,[file_prefix 'result.mat']);
if isfile(result_path)
    loaded=load(result_path,'result');
    if ~isfield(loaded,'result')||~isstruct(loaded.result)||~isscalar(loaded.result)
        error('DSC:PathIdentity','The saved sampler result is invalid.');
    end
    run_summary=loaded.result;
elseif exist(summary_path,'file')
    run_summary = jsondecode(fileread(summary_path));
else
    run_summary=struct();
    if isfile(manifest_path)
        if isfield(manifest,'run_summary')
            if ~isstruct(manifest.run_summary)||~isscalar(manifest.run_summary)
                error('DSC:PathIdentity','The manifest run summary is invalid.');
            end
            run_summary=manifest.run_summary;
        end
    end
end
[retention,run_summary]=read_retention_settings(manifest,parameters,run_summary);
[P,identity,iterations,draw_chain_ids,chain_dirs,chain_ids,retained_draws_per_chain]= ...
    read_chain_draws(run_dir,run_summary,file_prefix,parameters,retention);
if size(P,3)<2, error('DSC:TooFewDraws','At least two retained draws are needed.'); end
dates=identity.dates; tickers=string(identity.tickers);
pair_i=identity.pair_i; pair_j=identity.pair_j;
warmup=max(0,min(iterations)-1);
if isfield(run_summary,'burnin'), warmup=run_summary.burnin; end
if ~isempty(retention.burnin), warmup=retention.burnin; end
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
q=prctile(P,[16 84],3);
lower=q(:,:,1); mean_path=mean(P,3); upper=q(:,:,2);
chain_means=zeros(size(P,1),size(P,2),numel(chain_ids),'like',P);
for c=1:numel(chain_ids)
    chain_means(:,:,c)=mean(P(:,:,draw_chain_ids==chain_ids(c)),3);
end
inference=read_inference(run_dir,dates,run_summary,file_prefix,manifest,parameters);
if ~isempty(inference)&& ...
        ((isfield(inference,'num_chains')&&inference.num_chains~=numel(chain_ids))|| ...
        (isfield(inference,'retained_draws_per_chain')&& ...
        ~isequal(inference.retained_draws_per_chain(:),retained_draws_per_chain(:))))
    error('DSC:InferenceMetadata','Inference chain counts disagree with the saved posterior draws.');
end
draw_note=sprintf('%d pooled retained draws from %d chain(s) | %d warm-up iterations removed per chain', ...
    size(P,3),numel(chain_ids),warmup);

fig_dir = fullfile(run_dir,'figures');
compact_output=~isempty(file_prefix)||is_final_compact_output(run_summary);
if compact_output||~isfolder(fig_dir), fig_dir=run_dir; end
if ~compact_output
    save(fullfile(run_dir,[file_prefix 'posterior_correlation_bands.mat']), ...
        'lower','mean_path','upper','dates','tickers','pair_i','pair_j','iterations','inference', ...
        'draw_chain_ids','chain_ids','chain_dirs','retained_draws_per_chain','warmup_per_chain','-v7.3');
end

pair_labels = erase(tickers(pair_i)," Index") + " / " + erase(tickers(pair_j)," Index");
pages = ceil(numel(pair_labels)/9);
figure_paths = string(fullfile(fig_dir,[file_prefix 'correlation_paths.pdf']));
if numel(chain_ids)>1
    figure_paths(2,1)=string(fullfile(fig_dir,[file_prefix 'convergence_correlation_paths.pdf']));
end
plot_summary=struct('pooled_mean',mean_path,'lower',lower,'upper',upper, ...
    'chain_means',chain_means,'chain_ids',chain_ids, ...
    'retained_draws_per_chain',retained_draws_per_chain, ...
    'pooled_y_limits',nan(numel(pair_labels),2));
if numel(chain_ids)>1
    plot_summary.chain_y_limits=nan(numel(pair_labels),2);
end

for page = 1:pages
    ix = (page-1)*9 + (1:9);
    ix = ix(ix<=numel(pair_labels));
    pooled_y_limits=plot_page(figure_paths(1),dates,mean_path,lower,upper,pair_labels,ix, ...
        sprintf('Pooled posterior mean correlation paths | Page %d of %d',page,pages), ...
        inference,draw_note,page>1,[],[]);
    plot_summary.pooled_y_limits(ix,:)=pooled_y_limits;
    if numel(chain_ids)>1
        chain_y_limits=plot_page(figure_paths(2),dates,mean_path,lower,upper,pair_labels,ix, ...
            sprintf('Convergence: mean correlation paths by chain | Page %d of %d',page,pages), ...
            inference,draw_note,page>1,chain_means,chain_ids);
        plot_summary.chain_y_limits(ix,:)=chain_y_limits;
    end
end
if numel(chain_ids)>1
    dsc_plot_run_tables(figure_paths(2),run_dir,file_prefix);
end
end

function [P,identity,iterations,draw_ids,dirs,ids,counts]=read_chain_draws(run_dir,summary,default_prefix,parameters,retention)
dirs={run_dir};
if isfield(summary,'chain_dirs')
    dirs=summary.chain_dirs;
    if ischar(dirs)||isstring(dirs), dirs=cellstr(dirs); end
end
if ~iscell(dirs)||isempty(dirs)||~all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),dirs))
    error('DSC:PathIdentity','The sampler summary must contain nonempty chain directory paths.');
end
if ~isfield(summary,'chain_dirs')&&isfield(summary,'chain_file_prefixes')
    prefixes=summary.chain_file_prefixes;
    if ischar(prefixes)||isstring(prefixes), prefixes=cellstr(prefixes); end
    summary.chain_file_prefixes=prefixes;
    dirs=repmat({run_dir},1,numel(prefixes));
elseif ~isfield(summary,'chain_file_prefixes')
    if ~isfield(summary,'chain_dirs')|| ...
            all(cellfun(@(p)strcmp(char(java.io.File(p).getCanonicalPath()), ...
            char(java.io.File(run_dir).getCanonicalPath())),dirs))
        [dirs,summary.chain_file_prefixes]=discover_chain_files(run_dir,default_prefix);
    end
end
if ~iscell(dirs)||isempty(dirs)||~all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),dirs))
    error('DSC:PathIdentity','The sampler summary must contain nonempty chain directory paths.');
end
dirs=dirs(:)'; canonical=cell(size(dirs));
prefixes=repmat({default_prefix},size(dirs));
if isfield(summary,'chain_file_prefixes')
    prefixes=summary.chain_file_prefixes;
    if isstring(prefixes), prefixes=cellstr(prefixes); end
    if ~iscell(prefixes)||numel(prefixes)~=numel(dirs)
        error('DSC:PathIdentity','A filename prefix is required for each original chain.');
    end
    prefixes=cellfun(@validate_file_prefix,prefixes(:)','UniformOutput',false);
end
for c=1:numel(dirs)
    file=java.io.File(dirs{c}); canonical{c}=char(file.getCanonicalPath());
end
keys=cellfun(@(p,f)[p filesep f],canonical,prefixes,'UniformOutput',false);
if numel(unique(keys))~=numel(dirs)
    error('DSC:PathIdentity','Each posterior chain directory and filename-prefix identity must be distinct.');
end
% Use saved summary or MAT counts to allocate the pooled cube once rather
% than copying its accumulated history for every loaded chunk.
declared_counts=[]; declared_total=[];
if isfield(summary,'saved_draws_per_chain')
    declared_counts=summary.saved_draws_per_chain;
    if ~isnumeric(declared_counts)||~isreal(declared_counts)|| ...
            numel(declared_counts)~=numel(dirs)||any(~isfinite(declared_counts(:)))|| ...
            any(declared_counts(:)<1)||any(declared_counts(:)~=floor(declared_counts(:)))
        error('DSC:PathDraws','Invalid declared retained draw counts per chain.');
    end
    declared_counts=declared_counts(:)';
    declared_total=sum(declared_counts);
end
if isfield(summary,'saved_draws')
    value=summary.saved_draws;
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)|| ...
            value<1||value~=floor(value)||(~isempty(declared_total)&&value~=declared_total)
        error('DSC:PathDraws','Invalid or inconsistent declared total retained draw count.');
    end
    declared_total=value;
end
if isfield(parameters,'retained_draws')
    value=parameters.retained_draws;
    if ~valid_count(value,1)||(~isempty(declared_total)&&declared_total~=value)
        error('DSC:PathDraws','Invalid or inconsistent retained draw count in the parameter report.');
    end
    declared_total=value;
end
P=[]; identity=[]; iterations=[]; draw_ids=[]; ids=nan(1,numel(dirs)); counts=zeros(1,numel(dirs));
offset=0; chunk_values={}; chunk_iterations={}; chunk_ids={};
for c=1:numel(dirs)
    files=dir(fullfile(dirs{c},[prefixes{c} 'posterior_chunk_*.mat']));
    if isempty(files), error('DSC:NoChunks','No posterior chunks found in %s',dirs{c}); end
    numbers=zeros(1,numel(files));
    for n=1:numel(files)
        token=regexp(files(n).name,['^' regexptranslate('escape',prefixes{c}) 'posterior_chunk_(\d+)\.mat$'],'tokens','once');
        if isempty(token), error('DSC:PathDraws','Invalid posterior chunk filename.'); end
        numbers(n)=str2double(token{1});
    end
    if ~isequal(sort(numbers),1:numel(numbers))
        error('DSC:PathDraws','Chunk indices must be consecutive positive integers starting at one.');
    end
    [~,order]=sort(numbers); files=files(order); last_iteration=0;
    if ~isempty(retention.burnin), last_iteration=retention.burnin; end
    for n=1:numel(files)
        % A chunk also contains the full h paths and parameter draws. Release
        % the previous struct before load materializes the next chunk.
        clear loaded chunk values
        loaded=load(fullfile(dirs{c},files(n).name),'chunk');
        if ~isfield(loaded,'chunk'), error('DSC:PathIdentity','Missing posterior chunk.'); end
        chunk=loaded.chunk;
        fields={'dates','tickers','pair_i','pair_j','chain_id','iterations','P_pairs'};
        if ~isstruct(chunk)||~isscalar(chunk)||~all(isfield(chunk,fields))
            error('DSC:PathIdentity','A posterior chunk is missing required identity or draw fields.');
        end
        current=struct('dates',chunk.dates,'tickers',{chunk.tickers}, ...
            'pair_i',chunk.pair_i,'pair_j',chunk.pair_j,'compact', ...
            isfield(chunk,'mu')||(isfield(chunk,'schema')&&strcmp(chunk.schema,'constant_mu_v1')));
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
        if ~isempty(retention.thin)&& ...
                ((~isempty(retention.burnin)||n>1)&& ...
                any(diff([last_iteration index])~=retention.thin)|| ...
                any(diff(index)~=retention.thin))
            error('DSC:PathDraws','Saved iterations disagree with per-chain warm-up and thinning.');
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
        next_offset=offset+draws;
        if ~isempty(declared_total)
            if next_offset>declared_total|| ...
                    (~isempty(declared_counts)&&counts(c)+draws>declared_counts(c))
                error('DSC:PathDraws','Loaded posterior draws exceed the declared retained count.');
            end
            if offset==0
                P=zeros(T,nr,declared_total,'like',values);
                iterations=zeros(1,declared_total,'like',index);
                draw_ids=zeros(1,declared_total);
            elseif ~strcmp(class(P),class(values))
                % Determine MATLAB's concatenation class from empty arrays
                % when reading historical mixed-class chunks.
                prototype=cat(1,zeros(0,1,'like',P),zeros(0,1,'like',values));
                if ~strcmp(class(P),class(prototype)), P=cast(P,'like',prototype); end
            end
            range=offset+1:next_offset;
            if offset>0&&~strcmp(class(iterations),class(index))
                prototype=cat(2,zeros(1,0,'like',iterations),zeros(1,0,'like',index));
                if ~strcmp(class(iterations),class(prototype))
                    iterations=cast(iterations,'like',prototype);
                end
            end
            P(:,:,range)=values; iterations(range)=index; draw_ids(range)=ids(c);
        else
            % Historical direct-chunk calls may lack a count summary. Hold
            % references and concatenate once after reading every chunk.
            chunk_values{end+1}=values; %#ok<AGROW>
            chunk_iterations{end+1}=index; %#ok<AGROW>
            chunk_ids{end+1}=repmat(ids(c),1,draws); %#ok<AGROW>
        end
        offset=next_offset;
        counts(c)=counts(c)+draws; last_iteration=index(end);
    end
    if ~isempty(declared_counts)&&counts(c)~=declared_counts(c)
        error('DSC:PathDraws','Loaded posterior draws disagree with the declared chain count.');
    end
end
clear loaded chunk values
if isempty(declared_total)
    P=cat(3,chunk_values{:}); iterations=[chunk_iterations{:}]; draw_ids=[chunk_ids{:}];
elseif offset~=declared_total
    error('DSC:PathDraws','Loaded posterior draws disagree with the declared total count.');
end
if numel(unique(ids))~=numel(ids)
    error('DSC:PathIdentity','Posterior directories must have distinct original chain IDs.');
end
validate_parameter_chain_counts(parameters,ids,counts,identity);
end

function parameters=read_parameter_report(root,prefix,manifest)
parameters=struct(); path=fullfile(root,[prefix 'parameters.mat']);
if isfield(manifest,'outputs')&&isstruct(manifest.outputs)&&isfield(manifest.outputs,'parameters')
    value=manifest.outputs.parameters;
    if ~((ischar(value)&&isrow(value)&&~isempty(value))|| ...
            (isstring(value)&&isscalar(value)&&strlength(value)>0))
        error('DSC:PathIdentity','The saved parameter path must be scalar text.');
    end
    path=char(value);
    if ~java.io.File(path).isAbsolute(), path=fullfile(root,path); end
end
if ~isfile(path), return; end
loaded=load(path,'parameters');
if ~isfield(loaded,'parameters')||~isstruct(loaded.parameters)||~isscalar(loaded.parameters)
    error('DSC:PathIdentity','The saved parameter report is invalid.');
end
parameters=loaded.parameters;
end

function [settings,summary]=read_retention_settings(manifest,parameters,summary)
settings=struct('burnin',[],'thin',[]);
sources={summary};
if isfield(manifest,'config'), sources{end+1}=manifest.config; end
if isfield(parameters,'retention_config'), sources{end+1}=parameters.retention_config;
elseif isfield(parameters,'estimation_config'), sources{end+1}=parameters.estimation_config; end
names={{'burnin','WarmupIterations'},{'thin','Thin'}}; keys={'burnin','thin'};
for k=1:numel(keys)
    for s=1:numel(sources)
        source=sources{s};
        if ~isstruct(source)||~isscalar(source)
            error('DSC:PathDraws','Saved retention settings must be scalar structures.');
        end
        for f=1:numel(names{k})
            name=names{k}{f};
            if ~isfield(source,name), continue; end
            value=source.(name);
            if ~valid_count(value,k-1)|| ...
                    (~isempty(settings.(keys{k}))&&settings.(keys{k})~=value)
                error('DSC:PathDraws','Invalid or inconsistent saved warm-up/thinning settings.');
            end
            settings.(keys{k})=value;
        end
    end
end
% These are actual saved counts, unlike the requested public NumChains and
% PerDraw options. Keep them in memory; never restore removed JSON fields.
if isfield(parameters,'publication')&&isstruct(parameters.publication)&&isscalar(parameters.publication)
    for key={'saved_draws_per_chain','chain_ids','num_chains'}
        name=key{1};
        if ~isfield(summary,name)&&isfield(parameters.publication,name)
            summary.(name)=parameters.publication.(name);
        end
    end
end
end

function tf=valid_count(value,minimum)
tf=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&& ...
    value>=minimum&&value==floor(value);
end

function [dirs,prefixes]=discover_chain_files(root,prefix)
direct=dir(fullfile(root,[prefix 'posterior_chunk_*.mat']));
files=dir(fullfile(root,[prefix 'chain_*_posterior_chunk_*.mat']));
if ~isempty(direct)&&~isempty(files)
    error('DSC:PathIdentity','Both direct and chain-prefixed chunks match this run.');
end
if isempty(files), dirs={root}; prefixes={prefix}; return; end
prefixes=cell(1,numel(files)); ordinals=zeros(1,numel(files));
pattern=['^' regexptranslate('escape',prefix) 'chain_(\d+)_posterior_chunk_\d+\.mat$'];
for k=1:numel(files)
    token=regexp(files(k).name,pattern,'tokens','once');
    if isempty(token)||~valid_count(str2double(token{1}),1)
        error('DSC:PathIdentity','Invalid chain-prefixed posterior chunk filename.');
    end
    prefixes{k}=[prefix 'chain_' token{1} '_']; ordinals(k)=str2double(token{1});
end
[prefixes,first]=unique(prefixes,'stable'); ordinals=ordinals(first);
if numel(unique(ordinals))~=numel(ordinals)
    error('DSC:PathIdentity','Chain filename indices have ambiguous numeric identities.');
end
[~,order]=sort(ordinals); prefixes=prefixes(order);
dirs=repmat({root},size(prefixes));
end

function validate_parameter_chain_counts(parameters,ids,counts,identity)
if isfield(parameters,'retained_draws')&&parameters.retained_draws~=sum(counts)
    error('DSC:PathDraws','Posterior chunks disagree with the parameter report total count.');
end
if isfield(parameters,'chain_ids')||isfield(parameters,'retained_draws_per_chain')
    if ~all(isfield(parameters,{'chain_ids','retained_draws_per_chain'}))
        error('DSC:PathIdentity','The parameter report needs both chain IDs and retained counts.');
    end
    expected=parameters.chain_ids(:); expected_counts=parameters.retained_draws_per_chain(:);
    if ~isnumeric(expected)||~isreal(expected)||numel(expected)~=numel(ids)|| ...
            any(~isfinite(expected))||any(expected<1)||any(expected~=floor(expected))|| ...
            numel(unique(expected))~=numel(expected)||~isnumeric(expected_counts)|| ...
            ~isreal(expected_counts)||numel(expected_counts)~=numel(ids)|| ...
            any(~isfinite(expected_counts))||any(expected_counts<1)||any(expected_counts~=floor(expected_counts))
        error('DSC:PathIdentity','The parameter report has invalid chain identities or counts.');
    end
    [found,order]=ismember(ids,expected);
    if ~all(found)||~isequal(counts(:),expected_counts(order))
        error('DSC:PathDraws','Saved chunks disagree with parameter-report chain identities/counts.');
    end
end
for field={'tickers','pair_i','pair_j'}
    key=field{1};
    if isfield(parameters,key)&&~isequal(string(parameters.(key)(:)),string(identity.(key)(:)))
        error('DSC:PathIdentity','Posterior chunk series or pair order differs from the parameter report.');
    end
end
if isfield(parameters,'dates')&&~isequal(parameters.dates(:),identity.dates(:))
    error('DSC:PathIdentity','Posterior chunk dates differ from the parameter report.');
end
end

function y_limits=plot_page(path,dates,mean_path,lower,upper,pair_labels,indices,title_text,inference,draw_note,append_page,chain_means,chain_ids)
fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1200 1500]);
cleanup=onCleanup(@()close(fig));
layout = tiledlayout(numel(indices),1,'TileSpacing','compact','Padding','compact');
x = datenum(dates);
cmap = red_blue_map(257);
colormap(fig,cmap);
by_chain=~isempty(chain_ids);
colors=lines(numel(chain_ids));
if numel(chain_ids)>7, colors=turbo(numel(chain_ids)); end
y_limits=zeros(numel(indices),2);
for pos = 1:numel(indices)
    k = indices(pos);
    ax = nexttile; hold(ax,'on');
    if by_chain
        displayed=chain_means(:,k,:);
        handles=gobjects(numel(chain_ids),1);
        for c=1:numel(chain_ids)
            handles(c)=plot(ax,x,chain_means(:,k,c),'Color',colors(c,:),'LineWidth',1.1);
        end
    else
        displayed=[mean_path(:,k);lower(:,k);upper(:,k)];
    end
    extent=[min(displayed(:))-.1 max(displayed(:))+.1];
    if ~all(isfinite(extent))||extent(1)>=extent(2)
        error('DSC:PathLimits','Displayed correlation paths must have finite, nonempty ranges.');
    end
    y_limits(pos,:)=extent;
    if ~by_chain
        draw_bars(ax,x,mean_path(:,k),extent);
        fill(ax,[x; flipud(x)],[lower(:,k); flipud(upper(:,k))],[0.08 0.19 0.31], ...
            'FaceAlpha',0.16,'EdgeColor','none');
        plot(ax,x,mean_path(:,k),'Color',[0.02 0.10 0.20],'LineWidth',1.1);
    end
    yline(ax,0,'Color',[0.55 0.55 0.55],'LineWidth',0.5);
    if ~isempty(inference)
        cutoff=datenum(datetime(inference.parameter_estimation_end,'InputFormat','yyyy-MM-dd'));
        if cutoff>x(1)&&cutoff<x(end)
            xline(ax,cutoff,'--','Color',[.25 .25 .25],'LineWidth',.9,'HandleVisibility','off');
        end
    end
    ylim(ax,extent); xlim(ax,[x(1) x(end)]);
    clim(ax,[-1 1]);
    datetick(ax,'x','yyyy','keeplimits');
    ylabel(ax,'P_{ij,t}');
    title(ax,char(pair_labels(k)),'Interpreter','none','FontWeight','normal');
    box(ax,'off'); grid(ax,'on');
end
if by_chain
    labels=compose('Chain %d',chain_ids);
    lg=legend(ax,handles,labels,'Interpreter','none');
    lg.Layout.Tile='east';
    subtitle='Arithmetic mean across retained draws within each original chain';
else
    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Mean correlation';
    cb.Ticks = -1:0.25:1;
    subtitle='Posterior mean and 68% pointwise band (16th-84th percentiles)';
end
title(layout,{title_text,subtitle,draw_note}, ...
    'FontWeight','normal','Interpreter','none','FontSize',11);
exportgraphics(fig,path,'ContentType','vector','Append',append_page);
end

function inference=read_inference(run_dir,dates,summary,prefix,manifest,parameters)
inference=[];
metadata_path=fullfile(run_dir,[prefix 'inference_metadata.json']);
if isfield(manifest,'inference'), inference=manifest.inference; end
if isfile(metadata_path)
    staged=jsondecode(fileread(metadata_path));
    if isempty(inference), inference=staged;
    else, inference=merge_inference_fields(inference,staged); end
end
if isfield(parameters,'inference')
    if isempty(inference), inference=parameters.inference;
    else, inference=merge_inference_fields(inference,parameters.inference); end
end
if isempty(inference)
    if isfield(summary,'inference_mode')&&strcmp(summary.inference_mode,'fixed_parameter_smoothing')
        error('DSC:InferenceMetadata','Inference metadata is missing; cannot label conditional correlation paths.');
    end
    return
end
if ~isstruct(inference)||~isscalar(inference)|| ...
        ~all(isfield(inference,{'parameter_estimation_start','parameter_estimation_end'}))
    error('DSC:InferenceMetadata','Saved inference metadata must identify the parameter-estimation dates.');
end
% Compact public metadata contains dates. The private MAT report retains the
% inference target and smoothing scope; older complete JSON remains valid.
sources={parameters,summary};
if isfield(parameters,'inference'), sources=[{parameters.inference} sources]; end
for s=1:numel(sources)
    source=sources{s};
    if ~isstruct(source)||~isscalar(source)
        error('DSC:InferenceMetadata','Private inference metadata must be a scalar structure.');
    end
    for key={'inference_mode','parameter_smoothing','smoothing_start','smoothing_end'}
        name=key{1};
        if isfield(source,name)
            if isfield(inference,name)&&~isequaln(inference.(name),source.(name))
                error('DSC:InferenceMetadata','Private and public inference metadata disagree.');
            end
            inference.(name)=source.(name);
        end
    end
end
if ~isfield(inference,'smoothing_end')
    inference.smoothing_end=char(string(dates(end),'yyyy-MM-dd'));
end
if ~isfield(inference,'inference_mode')
    conditional=isfield(parameters,'storage')&&strcmp(parameters.storage,'conditional_parameter_summary_v1');
    if isfield(manifest,'config')&&isfield(manifest.config,'EstimateParameters')
        conditional=conditional||isequal(manifest.config.EstimateParameters,false);
    end
    conditional=conditional||~strcmp(inference.parameter_estimation_end,inference.smoothing_end);
    inference.inference_mode='parameter_estimation';
    if conditional, inference.inference_mode='fixed_parameter_smoothing'; end
end
if ~isfield(inference,'parameter_smoothing')
    if strcmp(inference.inference_mode,'parameter_estimation')
        inference.parameter_smoothing='draws';
    elseif isfield(manifest,'config')&&isfield(manifest.config,'ParameterSmoothing')
        inference.parameter_smoothing=manifest.config.ParameterSmoothing;
    else
        error('DSC:InferenceMetadata','The saved conditional run does not identify its parameter-smoothing scope.');
    end
end
for key={'parameter_estimation_start','parameter_estimation_end','smoothing_end','inference_mode','parameter_smoothing'}
    value=inference.(key{1});
    if ~((ischar(value)&&isrow(value)&&~isempty(value))|| ...
            (isstring(value)&&isscalar(value)&&~ismissing(value)&&strlength(value)>0))
        error('DSC:InferenceMetadata','Inference dates, mode and smoothing scope must be nonempty scalar text.');
    end
    inference.(key{1})=char(value);
end
mode=inference.inference_mode;
parameter_smoothing=inference.parameter_smoothing;
inference.parameter_uncertainty_in_bands= ...
    strcmp(mode,'parameter_estimation')||strcmp(parameter_smoothing,'draws');
if ~ismember(mode,{'parameter_estimation','fixed_parameter_smoothing'})|| ...
        ~ismember(parameter_smoothing,{'draws','mean'})|| ...
        (isfield(inference,'parameter_estimation_draws')&& ...
        ~valid_count(inference.parameter_estimation_draws,2))
    error('DSC:InferenceMetadata','Inference mode, retained parameter count or uncertainty scope is invalid.');
end
try
    start=datetime(inference.parameter_estimation_start,'InputFormat','yyyy-MM-dd');
    cutoff=datetime(inference.parameter_estimation_end,'InputFormat','yyyy-MM-dd');
    last=datetime(inference.smoothing_end,'InputFormat','yyyy-MM-dd');
    if isfield(inference,'parameters_estimated_at')
        estimated=datetime(inference.parameters_estimated_at,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'",'TimeZone','UTC');
        if isnat(estimated), error('DSC:InferenceMetadata','Invalid parameter-estimation timestamp.'); end
    end
catch
    error('DSC:InferenceMetadata','Invalid dates in inference metadata.');
end
if any(isnat([start cutoff last]))||start~=dates(1)||last~=dates(end)|| ...
        cutoff<start||cutoff>last||~any(dates==cutoff)|| ...
        (strcmp(mode,'parameter_estimation')&&cutoff~=last)|| ...
        (isfield(summary,'inference_mode')&&~strcmp(summary.inference_mode,mode))|| ...
        (isfield(summary,'first_date')&&~strcmp(summary.first_date,char(string(dates(1),'yyyy-MM-dd'))))|| ...
        (isfield(summary,'last_date')&&~strcmp(summary.last_date,char(string(dates(end),'yyyy-MM-dd'))))
    error('DSC:InferenceMetadata','Inference dates or mode disagree with the saved correlation sample.');
end
if isfield(inference,'calibration_weeks')
    n=inference.calibration_weeks;
    if ~valid_count(n,1)|| ...
            ~all(isfield(inference,{'calibration_start','calibration_end'}))|| ...
            (isfield(inference,'smoothing_start')&& ...
            ~strcmp(inference.smoothing_start,inference.parameter_estimation_start))|| ...
            (isfield(inference,'full_data_start')&& ...
            ~strcmp(inference.full_data_start,inference.calibration_start))|| ...
            (isfield(summary,'excluded_initial_weeks')&&summary.excluded_initial_weeks~=n)
        error('DSC:InferenceMetadata','Calibration exclusion metadata is inconsistent.');
    end
    initial=datetime(inference.calibration_start,'InputFormat','yyyy-MM-dd');
    final=datetime(inference.calibration_end,'InputFormat','yyyy-MM-dd');
    if isnat(initial)||isnat(final)||start~=final+calweeks(1)||final~=initial+calweeks(n-1)
        error('DSC:InferenceMetadata','Calibration and model dates must be adjacent, non-overlapping weekly blocks.');
    end
end
end

function merged=merge_inference_fields(public,private)
if ~isstruct(public)||~isscalar(public)||~isstruct(private)||~isscalar(private)
    error('DSC:InferenceMetadata','Saved inference metadata must be scalar structures.');
end
merged=private; fields=fieldnames(public);
for k=1:numel(fields)
    name=fields{k};
    if isfield(private,name)&&~isequaln(public.(name),private.(name))
        error('DSC:InferenceMetadata','Saved inference metadata disagree on %s.',name);
    end
    merged.(name)=public.(name);
end
end

function tf=is_final_compact_output(summary)
tf=isfield(summary,'final_output_schema')&& ...
    ismember(summary.final_output_schema,{'compact_outputs_v1','joint_outputs_v1'});
end

function prefix=validate_file_prefix(prefix)
if ~((ischar(prefix)&&(isrow(prefix)||isempty(prefix)))||(isstring(prefix)&&isscalar(prefix)))
    error('DSC:PathIdentity','Filename prefixes must be scalar text.');
end
prefix=char(prefix);
if ~isempty(prefix)&&isempty(regexp(prefix,'^[A-Za-z0-9_-]+$','once'))
    error('DSC:PathIdentity','Filename prefixes may contain only letters, numbers, underscores and hyphens.');
end
end

function draw_bars(ax,x,y,y_limits)
edges = [x(1); (x(1:end-1)+x(2:end))/2; x(end)];
image(ax,'XData',[edges(1) edges(end)],'YData',y_limits, ...
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

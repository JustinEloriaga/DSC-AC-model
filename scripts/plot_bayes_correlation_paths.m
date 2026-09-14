function figure_paths = plot_bayes_correlation_paths(run_dir)
%PLOT_BAYES_CORRELATION_PATHS Plot short-run Bayesian P(i,j,t) paths.
% Uses retained post-warm-up draws from posterior_chunk_*.mat. The colored
% vertical bands use Figure 1's continuous red-white-blue correlation scale.
if nargin<1 || strlength(string(run_dir))==0
    error('DSC:RunDir','Provide a sampler run directory.');
end
run_dir = char(run_dir);
files = dir(fullfile(run_dir,'posterior_chunk_*.mat'));
if isempty(files), error('DSC:NoChunks','No posterior chunks found in %s',run_dir); end
[~,order] = sort({files.name}); files = files(order);

P = []; iterations = [];
for n = 1:numel(files)
    loaded = load(fullfile(run_dir,files(n).name),'chunk');
    chunk = loaded.chunk;
    P = cat(3,P,chunk.P_pairs);
    iterations = [iterations, chunk.iterations]; %#ok<AGROW>
end
if size(P,3)<2, error('DSC:TooFewDraws','At least two retained draws are needed.'); end
dates = chunk.dates;
tickers = string(chunk.tickers);
pair_i = chunk.pair_i;
pair_j = chunk.pair_j;
q = prctile(P,[16 50 84],3);
lower = q(:,:,1); median_path = q(:,:,2); upper = q(:,:,3);
summary_path = fullfile(run_dir,'pilot_summary.json');
if exist(summary_path,'file')
    run_summary = jsondecode(fileread(summary_path));
    warmup = run_summary.burnin;
    convergence_established = logical(run_summary.convergence_established);
else
    warmup = max(0,min(iterations)-1);
    convergence_established = false;
end

fig_dir = fullfile(run_dir,'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
save(fullfile(run_dir,'posterior_correlation_bands.mat'), ...
    'lower','median_path','upper','dates','tickers','pair_i','pair_j','iterations','-v7.3');

pair_labels = tickers(pair_i) + " / " + tickers(pair_j);
selected = select_pairs(median_path,pair_labels);
pages = ceil(numel(pair_labels)/9);
figure_paths = strings(pages+1,1);
figure_paths(1) = plot_page(fullfile(fig_dir,'bayes_correlation_paths_selected.pdf'), ...
    dates, median_path, lower, upper, pair_labels, selected, ...
    'Selected Bayesian correlation paths');

for page = 1:pages
    ix = (page-1)*9 + (1:9);
    ix = ix(ix<=numel(pair_labels));
    figure_paths(page+1) = plot_page(fullfile(fig_dir,sprintf('bayes_correlation_paths_all_%02d.pdf',page)), ...
        dates, median_path, lower, upper, pair_labels, ix, ...
        sprintf('All Bayesian correlation paths, page %d of %d',page,pages));
end
write_json_local(fullfile(fig_dir,'bayes_correlation_paths_manifest.json'), ...
    struct('run_dir',run_dir,'retained_draws',size(P,3),'warmup_removed',warmup, ...
    'paths',{cellstr(figure_paths)},'color_scale',[-1 1], ...
    'color_rule','Continuous RdBu_r scale applied to the posterior median correlation.', ...
    'created_at',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')), ...
    'note',sprintf('%d retained draws after %d warm-up iterations; pointwise 16/50/84 percentiles; convergence established: %s.', ...
    size(P,3),warmup,boolean_text(convergence_established))));
end

function selected = select_pairs(median_path,pair_labels)
range = max(median_path,[],1) - min(median_path,[],1);
score = max(abs(median_path),[],1) + range;
[~,order] = sort(score,'descend');
must = find(pair_labels=="TUKXG / SX5T" | pair_labels=="SX5T / BCOMCOT");
selected = unique([must(:)' order(1:min(6,numel(order)))],'stable');
selected = selected(1:min(6,numel(selected)));
end

function path = plot_page(path,dates,median_path,lower,upper,pair_labels,indices,title_text)
fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1200 1500]);
layout = tiledlayout(numel(indices),1,'TileSpacing','compact','Padding','compact');
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
title(layout,{title_text,'Median and 68% pointwise band | Vertical colors show the median correlation value', ...
    'Short-run estimates from retained draws | Convergence not established'},'FontWeight','bold');
exportgraphics(fig,path,'ContentType','vector');
close(fig);
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

function write_json_local(path,value)
fid = fopen(path,'w');
assert(fid>=0,'DSC:WriteFailed','Cannot open %s',path);
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

function value = boolean_text(value)
if value, value='true'; else, value='false'; end
end

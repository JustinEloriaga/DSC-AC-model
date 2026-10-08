function output_path = plot_chain_mean_h(run_dir)
% Plot chain-specific posterior mean log-volatility paths from saved chunks.

if nargin < 1 || strlength(string(run_dir)) == 0
    run_dir = fullfile(pwd,'outputs','weekly_research_model1', ...
        'runs','20261005-224836-653-chain1-w100-r50000-thin1');
end
run_dir = char(run_dir);

chain_dirs = arrayfun(@(c) fullfile(run_dir,'chains',sprintf('chain_%03d',c)), ...
    1:4,'UniformOutput',false);
colors = lines(4);

for c = 1:4
    files = dir(fullfile(chain_dirs{c},'posterior_chunk_*.mat'));
    if isempty(files)
        error('plot:MissingChunks','No posterior chunks found for chain %d.',c);
    end
    [~,order] = sort({files.name});
    files = files(order);
    first = load(fullfile(files(1).folder,files(1).name),'chunk');
    chunk = first.chunk;
    T = size(chunk.h,1);
    m = size(chunk.h,2);
    sums = zeros(T,m);
    n_draws = 0;
    dates = chunk.dates;
    for f = 1:numel(files)
        loaded = load(fullfile(files(f).folder,files(f).name),'chunk');
        h = loaded.chunk.h;
        sums = sums + sum(h,3);
        n_draws = n_draws + size(h,3);
    end
    means{c} = sums./n_draws;
    draw_counts(c) = n_draws;
end

if any(draw_counts ~= draw_counts(1))
    error('plot:UnequalChains','Chains have unequal retained draw counts.');
end

tickers = string(first.chunk.tickers);
figure('Color','w','Position',[100 100 1400 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for j = 1:m
    nexttile;
    hold on;
    for c = 1:4
        plot(means{c}(:,j),dates,'LineWidth',1.1,'Color',colors(c,:));
    end
    grid on; box on;
    xlabel('Posterior mean h');
    ylabel('t');
    title(sprintf('h for %s',tickers(j)),'Interpreter','none');
    if j == 1
        legend({'Chain 1','Chain 2','Chain 3','Chain 4'}, ...
            'Location','best','Box','off');
    end
end
sgtitle(sprintf('Chain-specific posterior mean log-volatility paths (%d draws per chain)',draw_counts(1)), ...
    'FontWeight','bold');

output_path = fullfile(run_dir,'mean_h_by_chain.png');
exportgraphics(gcf,output_path,'Resolution',180);
close(gcf);
fprintf('Saved %s\n',output_path);
end

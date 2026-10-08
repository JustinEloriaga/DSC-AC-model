function output_path = plot_chain_mean_variances(run_dir)
% Plot chain-specific post-warm-up means of variance parameters.

if nargin < 1 || strlength(string(run_dir)) == 0
    run_dir = fullfile(pwd,'outputs','weekly_research_model1', ...
        'runs','20261005-224836-653-chain1-w100-r50000-thin1');
end
run_dir = char(run_dir);
chain_dirs = arrayfun(@(c) fullfile(run_dir,'chains',sprintf('chain_%03d',c)), ...
    1:4,'UniformOutput',false);
colors = lines(4);

for c = 1:4
    loaded = load(fullfile(chain_dirs{c},'diagnostics.mat'),'diagnostics');
    d = loaded.diagnostics;
    keep = d.iteration > 100;
    values{c} = [mean(d.V_diag(keep,:),1), mean(d.sig2h(keep,:),1), ...
        mean(d.sig2r(keep,:),1)];
    counts(c) = sum(keep);
end

labels = ["V(1,1)" "V(2,2)" "V(3,3)" "V(4,4)" ...
    "sig2h(1)" "sig2h(2)" "sig2h(3)" "sig2h(4)" ...
    "sig2r(1)" "sig2r(2)" "sig2r(3)" "sig2r(4)" "sig2r(5)" "sig2r(6)"];
groups = {1:4,'diag(V)'; 5:8,'sig2h'; 9:14,'sig2r'};

figure('Color','w','Position',[100 100 1500 650]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for g = 1:size(groups,1)
    nexttile;
    idx = groups{g,1};
    matrix = vertcat(values{:});
    bar(matrix(:,idx)','grouped');
    set(gca,'YScale','log','XTick',1:numel(idx),'XTickLabel',labels(idx));
    xtickangle(35); grid on; box on;
    ylabel('Posterior mean variance');
    title(groups{g,2},'Interpreter','none');
    if g == 1
        legend({'Chain 1','Chain 2','Chain 3','Chain 4'}, ...
            'Location','best','Box','off');
    end
end
sgtitle(sprintf('Chain-specific posterior mean variances (%d post-warm-up iterations per chain)',counts(1)), ...
    'FontWeight','bold');

output_path = fullfile(run_dir,'mean_variances_by_chain.png');
exportgraphics(gcf,output_path,'Resolution',180);
close(gcf);
fprintf('Saved %s\n',output_path);
end

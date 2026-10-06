function paths = write_convergence_artifacts(run_dir,result)
%WRITE_CONVERGENCE_ARTIFACTS Save chain-mean paths and a text diagnostic summary.
% Means are averaged over P pairs or h variables at each modeled date, then
% over retained draws within each chain. The convergence metrics are copied
% from dsc_diagnose_chains rather than recomputed here.
outdir = fullfile(char(run_dir),'convergence');
if ~isfolder(outdir), mkdir(outdir); end
chain_dirs = result.chain_dirs;
if ischar(chain_dirs) || isstring(chain_dirs), chain_dirs=cellstr(chain_dirs); end
nChains = numel(chain_dirs);
meanP=[]; meanH=[]; dates=[]; counts=zeros(1,nChains);
for c=1:nChains
    chunks=dir(fullfile(chain_dirs{c},'posterior_chunk_*.mat'));
    [~,order]=sort({chunks.name}); chunks=chunks(order);
    sumP=[]; sumH=[];
    for k=1:numel(chunks)
        loaded=load(fullfile(chunks(k).folder,chunks(k).name),'chunk');
        chunk=loaded.chunk;
        if isempty(dates), dates=chunk.dates(:); end
        h=double(chunk.h); p=double(chunk.P_pairs);
        if isempty(sumP)
            sumP=zeros(size(p,1),1); sumH=zeros(size(h,1),1);
        end
        sumP=sumP+reshape(sum(sum(p,2),3),[],1);
        sumH=sumH+reshape(sum(sum(h,2),3),[],1);
        counts(c)=counts(c)+size(p,3);
    end
    if counts(c)>0
        meanP(:,c)=sumP/counts(c); %#ok<AGROW>
        meanH(:,c)=sumH/counts(c); %#ok<AGROW>
    end
end

pdfPath=fullfile(outdir,'chain_mean_paths.pdf');
fig=figure('Visible','off','Color','w','Position',[100 100 1200 850]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(dates,meanP,'LineWidth',1); grid on;
ylabel('Mean P'); title('Mean correlation across P pairs by chain');
legend(compose('Chain %d',1:nChains),'Location','eastoutside');
nexttile; plot(dates,meanH,'LineWidth',1); grid on;
ylabel('Mean h'); xlabel('Date'); title('Mean h across variables by chain');
legend(compose('Chain %d',1:nChains),'Location','eastoutside');
exportgraphics(fig,pdfPath,'ContentType','vector'); close(fig);

txtPath=fullfile(outdir,'convergence_summary.txt');
fid=fopen(txtPath,'w');
if fid<0, error('dsc:ConvergenceSummary','Cannot write %s.',txtPath); end
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'CONVERGENCE SUMMARY\n====================\n');
fprintf(fid,'Scope: %s\n\n',result.convergence.scope);
fprintf(fid,'STATUS\n------\n%s\n%s\nPassed: %s\n\n', ...
    result.convergence.status,result.convergence.reason,logical_text(result.convergence.passed));
fprintf(fid,'CHAIN DRAW COUNTS (retained draws per chain)\n----------------------------------------------\n');
for c=1:nChains, fprintf(fid,'Chain %d: %d\n',c,counts(c)); end
fprintf(fid,'Total retained draws: %d\n\n',sum(counts));
fprintf(fid,'DIAGNOSTIC EXTREMES\n--------------------\n');
fprintf(fid,'Maximum R-hat: %.8g (threshold %.8g)\n',result.convergence.max_rhat,result.convergence.thresholds.rhat_threshold);
fprintf(fid,'Minimum bulk ESS: %.8g (threshold %.8g)\n',result.convergence.min_ess_bulk,result.convergence.thresholds.min_ess);
fprintf(fid,'Minimum tail ESS: %.8g (threshold %.8g)\n',result.convergence.min_ess_tail,result.convergence.thresholds.min_ess);
fprintf(fid,'Maximum MCSE / posterior SD: %.8g (threshold %.8g)\n', ...
    result.convergence.max_mcse_sd_ratio,result.convergence.thresholds.max_mcse_ratio);
fprintf(fid,'Quantities checked: %d\n\n',result.convergence.quantities_checked);
if isfield(result,'convergence_details') && istable(result.convergence_details)
    d=result.convergence_details;
    fprintf(fid,'QUANTITY-LEVEL SUMMARY\n-----------------------\n');
    fprintf(fid,'Family\tQuantity\tRhat\tBulkESS\tTailESS\tMCSE/SD\tPassed\n');
    for k=1:height(d)
        fprintf(fid,'%s\t%s\t%.8g\t%.8g\t%.8g\t%.8g\t%s\n', ...
            char(d.family(k)),char(d.quantity(k)),d.rhat(k),d.ess_bulk(k), ...
            d.ess_tail(k),d.mcse_sd_ratio(k),logical_text(d.passed(k)));
    end
end
paths=struct('directory',outdir,'pdf',pdfPath,'summary',txtPath);
end

function value=logical_text(flag)
if flag, value='true'; else, value='false'; end
end

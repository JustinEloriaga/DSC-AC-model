function paths = write_convergence_artifacts(run_dir,result)
%WRITE_CONVERGENCE_ARTIFACTS Save entry-by-entry chain-mean paths and a text summary.
% Each panel is one P pair or one h variable, with all chains over time.
% P is a correlation and is displayed on [-1,1]. h is log standard deviation.
outdir = fullfile(char(run_dir),'convergence');
if ~isfolder(outdir), mkdir(outdir); end
chain_dirs = result.chain_dirs;
if ischar(chain_dirs) || isstring(chain_dirs), chain_dirs=cellstr(chain_dirs); end
nChains = numel(chain_dirs);
meanP=[]; meanH=[]; dates=[]; counts=zeros(1,nChains); pair_i=[]; pair_j=[]; tickers=[];
for c=1:nChains
    chunks=dir(fullfile(chain_dirs{c},'posterior_chunk_*.mat'));
    [~,order]=sort({chunks.name}); chunks=chunks(order);
    sumP=[]; sumH=[];
    for k=1:numel(chunks)
        loaded=load(fullfile(chunks(k).folder,chunks(k).name),'chunk');
        chunk=loaded.chunk;
        if isempty(dates)
            dates=chunk.dates(:); pair_i=chunk.pair_i; pair_j=chunk.pair_j; tickers=string(chunk.tickers);
        end
        h=double(chunk.h); p=double(chunk.P_pairs);
        if isempty(sumP)
            sumP=zeros(size(p,1),size(p,2)); sumH=zeros(size(h,1),size(h,2));
        end
        sumP=sumP+sum(p,3);
        sumH=sumH+sum(h,3);
        counts(c)=counts(c)+size(p,3);
    end
    if counts(c)>0
        meanP(:,:,c)=sumP/counts(c); %#ok<AGROW>
        meanH(:,:,c)=sumH/counts(c); %#ok<AGROW>
    end
end

pdfPath=fullfile(outdir,'chain_mean_paths.pdf');
if isfile(pdfPath), delete(pdfPath); end
labels=erase(tickers," Index");
entries=size(meanP,2)+size(meanH,2);
for page=1:ceil(entries/10)
    indices=(page-1)*10+(1:min(10,entries-(page-1)*10));
    fig=figure('Visible','off','Color','w','Position',[50 50 1600 1100]);
    figure_guard=onCleanup(@()close(fig));
    tiledlayout(fig,ceil(numel(indices)/2),2,'TileSpacing','compact','Padding','compact');
    for q=indices
        nexttile;
        if q<=size(meanP,2)
            plot(dates,squeeze(meanP(:,q,:)),'LineWidth',0.8); ylim([-1 1]);
            ylabel('P');
            title(sprintf('Mean P: %s / %s',labels(pair_i(q)),labels(pair_j(q))));
        else
            j=q-size(meanP,2);
            plot(dates,squeeze(meanH(:,j,:)),'LineWidth',0.8);
            ylabel('h'); title(sprintf('Mean h: %s',labels(j)));
        end
        grid on; xlabel('t');
        if q==indices(1), legend(compose('Chain %d',result.chain_ids),'Location','best'); end
    end
    sgtitle('Chain means by entry: P correlations and h log standard deviations');
    exportgraphics(fig,pdfPath,'ContentType','vector','Append',page>1);
    clear figure_guard
end

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

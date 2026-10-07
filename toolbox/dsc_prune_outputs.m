function dsc_prune_outputs(output_root,report_dir,current_run)
% Keep only the newest completed runs and run-linked publication PDFs.
run_root=fullfile(char(output_root),'runs');
if isfolder(run_root)
    entries=dir(run_root);
    entries=entries([entries.isdir]&~ismember({entries.name},{'.','..'}));
    names=string({entries.name});
    [~,order]=sort(names,'descend');
    entries=entries(order);
    protected=false(size(entries));
    for k=1:numel(entries)
        path=fullfile(entries(k).folder,entries(k).name);
        protected(k)=isfile(fullfile(path,'.run_in_progress'))|| ...
            strcmp(path,char(current_run));
    end
    keep_completed=max(0,3-sum(protected));
    completed_kept=0;
    for k=1:numel(entries)
        if protected(k), continue; end
        if completed_kept<keep_completed
            completed_kept=completed_kept+1;
        else
            rmdir(fullfile(entries(k).folder,entries(k).name),'s');
        end
    end
end
prune_pdfs(char(report_dir),'report_*.pdf');
prune_pdfs(char(report_dir),'figure_*.pdf');
end

function prune_pdfs(report_dir,pattern)
if ~isfolder(report_dir), return; end
files=dir(fullfile(report_dir,pattern));
names=string({files.name});
[~,order]=sort(names,'descend');
for k=4:numel(order)
    delete(fullfile(files(order(k)).folder,files(order(k)).name));
end
end

function dsc_prune_outputs(output_root,~,current_run)
% Keep the newest three complete prefixed runs and every artifact they own.
% Active runs and the explicitly selected stopped/current run are protected.
% Incomplete, failed and unclassified runs are eligible for cleanup.
root=char(output_root);
if ~isfolder(root), return; end
manifests=dir(fullfile(root,'*_run_manifest.json'));
if isempty(manifests), return; end
[~,order]=sort({manifests.name},'descend'); manifests=manifests(order);
prefixes=cell(size(manifests));
for k=1:numel(manifests)
    name=manifests(k).name;
    prefixes{k}=name(1:end-numel('run_manifest.json'));
end
current=char(current_run);
if isfolder(current), [~,current]=fileparts(current); end
kept=0;
for k=1:numel(prefixes)
    p=prefixes{k};
    active=isfile(fullfile(root,[p '.run_in_progress']));
    selected=strcmp(p,current)||strcmp([p 'run_manifest.json'],current);
    if active||selected||kept<3
        if ~active&&~selected, kept=kept+1; end
    else
        files=dir(fullfile(root,[p '*']));
        for j=1:numel(files)
            if ~files(j).isdir, delete(fullfile(files(j).folder,files(j).name)); end
        end
    end
end
end

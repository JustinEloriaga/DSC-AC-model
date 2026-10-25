function report=dsc_prune_outputs(output_root)
%DSC_PRUNE_OUTPUTS Keep three complete groups plus active/stopped groups.
% Called only after publication. A run owns every flat file with its fixed
% timestamp prefix, including estimation-stage and temporary artifacts.
root=char(java.io.File(char(output_root)).getCanonicalPath());
assert(isfolder(root),'dsc:OutputDirectory','The output directory must already exist.');
prefixes=run_prefixes(root);
states=cell(size(prefixes)); manifests=cell(size(prefixes));
for k=1:numel(prefixes)
    [states{k},manifests{k}]=group_state(root,prefixes{k});
end
protected=ismember(states,{'active','stopped'});
% A protected smoothing invocation can need an older saved parameter bundle.
% Protect the entire source group, not just its parameters.mat file.
protected=protect_sources(root,prefixes,manifests,protected);
complete=find(strcmp(states,'completed'));
kept=complete(1:min(3,numel(complete)));
keep=protected; keep(kept)=true;
report=struct('kept_complete',{prefixes(kept)}, ...
    'protected',{prefixes(protected)},'deleted',{{}});
for k=find(~keep)
    % Recheck ownership immediately before mutation: another invocation may
    % have started publishing or registered a dependency since the inventory.
    [state,~]=group_state(root,prefixes{k});
    if any(strcmp(state,{'active','stopped'})), continue; end
    live_prefixes=run_prefixes(root);
    current=cell(size(live_prefixes)); current_protected=false(size(live_prefixes));
    for j=1:numel(live_prefixes)
        [current_state,current{j}]=group_state(root,live_prefixes{j});
        current_protected(j)=any(strcmp(current_state,{'active','stopped'}));
    end
    current_protected=protect_sources(root,live_prefixes,current,current_protected);
    if any(strcmp(live_prefixes(current_protected),prefixes{k})), continue; end
    owned=dir(fullfile(root,[prefixes{k} '*']));
    for j=1:numel(owned)
        if ~owned(j).isdir
            delete(fullfile(root,owned(j).name));
        end
    end
    report.deleted{end+1}=prefixes{k};
end
end

function prefixes=run_prefixes(root)
files=dir(root); prefixes={};
for k=1:numel(files)
    if files(k).isdir, continue; end
    prefix=regexp(files(k).name,'^constant_mu_\d{8}T\d{9}_','match','once');
    if ~isempty(prefix), prefixes{end+1}=prefix; end %#ok<AGROW>
end
prefixes=fliplr(sort(unique(prefixes)));
end

function [state,manifest]=group_state(root,prefix)
state='incomplete'; manifest=read_json(fullfile(root,[prefix 'run_manifest.json']));
owner=run_owner(root,prefix,manifest);
if isfield(owner,'status')
    status=char(owner.status);
    if strcmp(status,'active')
        if owner_may_be_running(owner), state='active'; return; end
        % A crash after the final manifest commits can leave an obsolete
        % lease. Only a verified complete publication can classify it below.
    elseif any(strcmp(status,{'stopped','failed'}))
        state=status; return
    end
end
if isfield(manifest,'run_summary'), summary=manifest.run_summary;
else, summary=read_json(fullfile(root,[prefix 'pilot_summary.json'])); end
if ~isfield(summary,'status')
    if complete_publication(manifest,summary,root,prefix), state='completed'; end
    return
end
switch char(summary.status)
    case {'time_limit','stopped'}, state='stopped';
    case 'failed', state='failed';
    case 'iteration_limit'
        if complete_publication(manifest,summary,root,prefix), state='completed'; end
end
end

function yes=owner_may_be_running(owner)
% Unknown ownership/probe errors protect the run. A dead or reused local PID
% does not turn an abandoned active record into permanent protection.
yes=true;
required={'process_id','host','process_started'};
if ~all(isfield(owner,required))||~isnumeric(owner.process_id)|| ...
        ~isscalar(owner.process_id)||~isfinite(owner.process_id)|| ...
        owner.process_id<1||owner.process_id~=floor(owner.process_id)
    return
end
if ~strcmp(owner.host,char(java.net.InetAddress.getLocalHost.getHostName))||~isunix
    return
end
[status,value]=system(sprintf('ps -p %.0f -o lstart= 2>&1',owner.process_id));
value=strtrim(value);
if status==0&&~isempty(value)
    yes=isempty(owner.process_started)||strcmp(value,owner.process_started);
elseif status==1&&isempty(value)
    yes=false;
end
end

function yes=complete_publication(manifest,summary,root,prefix)
yes=false;
if ~isfield(manifest,'config')
    return
end
num_chains=1;
if isfield(summary,'num_chains')&&isnumeric(summary.num_chains), num_chains=summary.num_chains; end
if isfield(manifest,'outputs')
    outputs=manifest.outputs;
    if ~all(isfield(outputs,{'parameters','joint_draw_chunks'})), return; end
    paths=[cellstr(string(outputs.parameters));cellstr(string(outputs.joint_draw_chunks(:)))];
else
    paths={fullfile(root,[prefix 'parameters.mat'])};
    chunks=dir(fullfile(root,[prefix '*posterior_chunk_*.mat']));
    if isempty(chunks), return; end
    paths=[paths; fullfile({chunks.folder},{chunks.name})'];
end
create_pdf=config_option(manifest,'CreatePDF','create_pdf',true);
% Publication metadata records whether this run actually requested figures;
% config accepts both current public option names and historical cfg names.
if isfield(manifest,'pdf')&&isstruct(manifest.pdf)&&isscalar(manifest.pdf)&& ...
        isfield(manifest.pdf,'enabled')
    create_pdf=manifest.pdf.enabled;
end
if create_pdf
    paths{end+1}=fullfile(root,[prefix 'correlation_paths.pdf']);
    if num_chains>1
        paths{end+1}=fullfile(root,[prefix 'convergence_correlation_paths.pdf']);
    end
end
for k=1:numel(paths)
    path=char(java.io.File(paths{k}).getCanonicalPath());
    [directory,name,extension]=fileparts(path);
    if ~strcmp(directory,root)||~startsWith([name extension],prefix)||~isfile(path)
        return
    end
end
yes=true;
end

function protected=protect_sources(root,prefixes,manifests,protected)
changed=true;
while changed
    before=protected;
    for k=find(protected)
        manifest=manifests{k};
        owner=run_owner(root,prefixes{k},manifest);
        paths={config_option(manifest,'ParameterFile','parameter_file','')};
        if isfield(owner,'parameter_file'), paths{end+1}=owner.parameter_file; end
        if isfield(manifest,'inference')&&isfield(manifest.inference,'parameter_file')
            paths{end+1}=manifest.inference.parameter_file;
        end
        % Register before loading parameters. Until automatic selection names
        % its source, preserve candidate bundles during this short setup phase.
        estimate_parameters=config_option(manifest,'EstimateParameters','estimate_parameters',[]);
        if isempty(estimate_parameters)&&isfield(owner,'estimate_parameters')
            estimate_parameters=owner.estimate_parameters;
        end
        loading=isequal(estimate_parameters,false)||isequal(estimate_parameters,0);
        selecting=isfield(owner,'status')&&strcmp(owner.status,'active')&& ...
            loading&&all(cellfun(@(p)strlength(string(p))==0,paths));
        if selecting
            for j=1:numel(prefixes)
                if isfile(fullfile(root,[prefixes{j} 'parameters.mat'])), protected(j)=true; end
            end
        end
        for j=1:numel(paths)
            if strlength(string(paths{j}))==0, continue; end
            path=char(paths{j});
            if ~isfile(path), path=fullfile(root,path); end
            path=char(java.io.File(path).getCanonicalPath());
            [directory,name,extension]=fileparts(path);
            if ~strcmp(directory,root), continue; end
            prefix=regexp([name extension],'^constant_mu_\d{8}T\d{9}_','match','once');
            protected(strcmp(prefixes,prefix))=true;
        end
    end
    changed=~isequal(before,protected);
end
end

function owner=run_owner(root,prefix,manifest)
state_path=fullfile(root,[prefix 'run_state.json']);
if isfile(state_path)
    owner=read_json(state_path);
    if ~isfield(owner,'status'), owner.status='active'; end
elseif isfield(manifest,'lifecycle')
    owner=manifest.lifecycle; % Older unfinished runs remain protected.
else
    owner=struct();
end
end

function value=config_option(manifest,public_name,legacy_name,fallback)
value=fallback;
if ~isfield(manifest,'config')||~isstruct(manifest.config)||~isscalar(manifest.config)
    return
end
if isfield(manifest.config,public_name)
    value=manifest.config.(public_name);
elseif isfield(manifest.config,legacy_name)
    value=manifest.config.(legacy_name);
end
end

function value=read_json(path)
value=struct();
if ~isfile(path), return; end
try
    decoded=jsondecode(fileread(path));
    if isstruct(decoded)&&isscalar(decoded), value=decoded; end
catch
    % Corrupt/unpublished metadata is incomplete, not a completed run.
end
end

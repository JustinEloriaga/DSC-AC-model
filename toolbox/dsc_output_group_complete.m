function [complete,publication]=dsc_output_group_complete(root,prefix,manifest)
%DSC_OUTPUT_GROUP_COMPLETE Verify a published flat run using MAT metadata.
% Public manifests describe the invocation; completion and chain counts live
% in the parameter report and must agree with the requested config.
complete=false; publication=struct();
if ~isstruct(manifest)||~isscalar(manifest)||~isfield(manifest,'config')
    return
end
path=fullfile(root,[prefix 'parameters.mat']);
if ~isfile(path), return; end
try
    loaded=load(path,'parameters'); p=loaded.parameters;
    if ~isstruct(p)||~isscalar(p)||~isfield(p,'storage')|| ...
            ~ismember(string(p.storage),["parameter_report_v1","conditional_parameter_summary_v1"])
        return
    end
    if isfield(p,'publication')
        publication=p.publication;
    elseif isfield(p,'estimation_status')
        publication.status=p.estimation_status;
    end
    if ~isfield(publication,'status')||~strcmp(publication.status,'iteration_limit')
        return
    end
    cfg=manifest.config;
    n=option(cfg,'NumChains','num_chains',[]);
    target=option(cfg,'PerDraw','',[]);
    if isempty(target)&&all(isfield(cfg,{'max_iterations','burnin','thin'}))
        target=floor((cfg.max_iterations-cfg.burnin)/cfg.thin);
    end
    if ~positive_integer(n)||~positive_integer(target)||target<2, return; end
    if isfield(publication,'num_chains')&&~isequal(publication.num_chains,n), return; end
    if isfield(publication,'saved_draws_per_chain')
        counts=publication.saved_draws_per_chain;
    elseif isfield(p,'retained_draws_per_chain')
        counts=p.retained_draws_per_chain;
    else
        return
    end
    if ~isnumeric(counts)||~isreal(counts)||numel(counts)~=n|| ...
            any(~isfinite(counts(:)))||any(counts(:)~=target)
        return
    end
    files=dir(fullfile(root,[prefix 'chain_*_posterior_chunk_*.mat']));
    files=files(~[files.isdir]); ids=zeros(1,numel(files));
    actual=zeros(1,n);
    for k=1:numel(files)
        token=regexp(files(k).name,['^' regexptranslate('escape',prefix) ...
            'chain_(\d+)_posterior_chunk_\d+\.mat$'],'tokens','once');
        if isempty(token)||files(k).bytes==0, return; end
        ids(k)=str2double(token{1});
        if ~ismember(ids(k),1:n), return; end
        % Read only the iteration dataset, never the large state/path cube.
        info=h5info(fullfile(root,files(k).name),'/chunk/iterations');
        actual(ids(k))=actual(ids(k))+prod(info.Dataspace.Size);
    end
    if ~isequal(unique(ids),1:n)||any(actual~=counts(:)'), return; end
    create_pdf=option(cfg,'CreatePDF','create_pdf',true);
    if create_pdf
        names={[prefix 'correlation_paths.pdf']};
        if n>1, names{end+1}=[prefix 'convergence_correlation_paths.pdf']; end
        for k=1:numel(names)
            file=dir(fullfile(root,names{k}));
            if isempty(file)||file.isdir||file.bytes==0, return; end
        end
    end
    complete=true;
catch
    % Invalid or half-published artifacts are never counted as completed.
end
end

function value=option(cfg,name,legacy,fallback)
value=fallback;
if isfield(cfg,name), value=cfg.(name);
elseif ~isempty(legacy)&&isfield(cfg,legacy), value=cfg.(legacy); end
end

function yes=positive_integer(value)
yes=isnumeric(value)&&isreal(value)&&isscalar(value)&& ...
    isfinite(value)&&value>=1&&value==floor(value);
end

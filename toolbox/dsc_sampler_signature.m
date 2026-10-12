function info=dsc_sampler_signature(backend)
% Refuse to resume with a different numerical implementation or native binary.
folder=fileparts(mfilename('fullpath'));
names={'dsc_sample.m','dsc_sampler_signature.m','dsc_corr_matrix.m', ...
    'dsc_corr_matrix_new.m','dsc_random_walk_draw.m','dsc_volatility_likelihood.m'};
mean_backend='matlab';
if strcmp(backend,'mex')||strcmp(backend,'mex_new')
    if strcmp(backend,'mex_new')
        corr_source='dsc_corr_batch_mex_new.cpp';
        corr_binary=['dsc_corr_batch_mex_new.' mexext];
    else
        corr_source='dsc_corr_batch_mex.cpp';
        corr_binary=['dsc_corr_batch_mex.' mexext];
    end
    names=[names,{corr_source,corr_binary}];
    if exist('dsc_mean_mex','file')==3
        mean_backend='mex';
        names=[names,{'dsc_mean_mex.cpp',['dsc_mean_mex.' mexext]}];
    end
end
md=java.security.MessageDigest.getInstance('SHA-256');
for k=1:numel(names)
    md.update(uint8(unicode2native(names{k},'UTF-8'))); md.update(uint8(0));
    fid=fopen(fullfile(folder,names{k}),'rb');
    assert(fid>=0,'dsc:ImplementationFile','Cannot hash sampler file %s',names{k});
    guard=onCleanup(@()fclose(fid));
    md.update(fread(fid,Inf,'*uint8')); md.update(uint8(0)); clear guard
end
info=struct('version',3,'correlation_backend',backend,'mean_backend',mean_backend,'matlab_version',version, ...
    'platform',computer,'source_and_binary_sha256', ...
    lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[])));
end

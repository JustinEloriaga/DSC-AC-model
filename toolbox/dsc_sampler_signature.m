function info=dsc_sampler_signature(backend)
% Refuse to resume with a different numerical implementation or native binary.
folder=fileparts(mfilename('fullpath'));
names={'dsc_sample.m','dsc_sampler_signature.m','dsc_corr_matrix.m','dsc_random_walk_draw.m'};
if strcmp(backend,'mex')
    names=[names,{'dsc_corr_batch_mex.cpp',['dsc_corr_batch_mex.' mexext]}];
end
md=java.security.MessageDigest.getInstance('SHA-256');
for k=1:numel(names)
    md.update(uint8(unicode2native(names{k},'UTF-8'))); md.update(uint8(0));
    fid=fopen(fullfile(folder,names{k}),'rb');
    assert(fid>=0,'dsc:ImplementationFile','Cannot hash sampler file %s',names{k});
    guard=onCleanup(@()fclose(fid));
    md.update(fread(fid,Inf,'*uint8')); md.update(uint8(0)); clear guard
end
info=struct('version',2,'correlation_backend',backend,'matlab_version',version, ...
    'platform',computer,'source_and_binary_sha256', ...
    lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[])));
end

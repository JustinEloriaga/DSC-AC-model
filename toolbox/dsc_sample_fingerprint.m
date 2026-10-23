function fingerprint = dsc_sample_fingerprint(block)
%DSC_SAMPLE_FINGERPRINT Canonical SHA-256 identity of an ordered weekly sample.
% Excludes the source CSV hash so appended observations do not change an
% existing history. Includes missing values, mask, dates, series and pair order.
required={'dates','returns','observation_mask','tickers','pair_i','pair_j'};
if ~isstruct(block)||~isscalar(block)||~all(isfield(block,required))|| ...
        ~isdatetime(block.dates)||any(isnat(block.dates))
    error('dsc:SampleFingerprint','A complete dated sample block is required.');
end
Y=double(block.returns); mask=logical(block.observation_mask);
if ~isreal(Y)||~isequal(size(Y),size(mask))||size(Y,1)~=numel(block.dates)|| ...
        size(Y,2)~=numel(block.tickers)||any(isinf(Y(:)))||any(~isfinite(Y(mask)))
    error('dsc:SampleFingerprint','Invalid sample dimensions, values or observation mask.');
end
metadata=struct('version',1,'dimensions',size(Y),'dates', ...
    {cellstr(string(block.dates(:),'yyyy-MM-dd'))}, ...
    'tickers',{cellstr(string(block.tickers(:)))}, ...
    'pair_i',double(block.pair_i(:))','pair_j',double(block.pair_j(:))');
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native(jsonencode(metadata),'UTF-8'));
Y(isnan(Y))=NaN; % Canonical NaN payload on every supported platform.
[~,~,endian]=computer;
if endian=='B', Y=swapbytes(Y); end
digest.update(typecast(Y(:),'uint8'));
digest.update(uint8(mask(:)));
bytes=typecast(digest.digest(),'uint8');
fingerprint=lower(reshape(dec2hex(bytes,2)',1,[]));
end

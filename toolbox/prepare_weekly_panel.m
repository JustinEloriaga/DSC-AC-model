function panel = prepare_weekly_panel(cfg)
%PREPARE_WEEKLY_PANEL Validate levels and construct Friday as-of log returns.
%   Each market uses only its own last published level at or before Friday.
%   No interpolation, future observation, common-market date, or calendar
%   compression is used. The sole cross-week carry is the documented Nikkei
%   closure ending 2019-05-03. Returns retain the regular weekly time grid.
%
%   cfg.source_file: Bloomberg-style CSV with three header rows.
%   cfg.start_date: earliest allowed raw observation (default 2003-08-04).
%   cfg.closure_policy: 'asof' (default) or 'mask'. With 'mask', values remain
%   available for audit but observation_mask excludes the two NKYTR returns
%   touching the exceptional closure level. The sampler must honor the mask.

arguments
    cfg (1,1) struct
end
if ~isfield(cfg, 'source_file') || ~isfile(cfg.source_file)
    error('weekly:SourceFile', 'cfg.source_file must identify an existing CSV.');
end
if ~isfield(cfg, 'start_date'), cfg.start_date = datetime(2003,8,4); end
if ~isfield(cfg, 'closure_policy'), cfg.closure_policy = 'asof'; end
startDate = datetime(cfg.start_date);
if ~isscalar(startDate) || isnat(startDate)
    error('weekly:StartDate', 'cfg.start_date must be a valid scalar date.');
end
policy = lower(string(cfg.closure_policy));
if ~isscalar(policy) || ~ismember(policy, ["asof", "mask"])
    error('weekly:ClosurePolicy', 'closure_policy must be asof or mask.');
end
canonical = ["AUDUSD Index", "BCOMCOT Index", "BCOMGCTR Index", ...
    "EURUSD Index", "GBPUSD Index", "I02981JP Index", "LEATTREU Index", ...
    "LSG1TRGU Index", "LUATTRUU Index", "NKYTR Index", "SPXT Index", ...
    "SX5T Index", "TUKXG Index", "USDJPY Index"];
sourceFile = char(java.io.File(char(cfg.source_file)).getCanonicalPath());
lines = readlines(sourceFile);
if numel(lines) < 5
    error('weekly:CSVFormat', 'Expected three header rows followed by observations.');
end
header = strtrim(split(lines(2), ','));
rawTickers = header(2:end)';
if numel(unique(rawTickers)) ~= numel(rawTickers)
    error('weekly:Tickers', 'Ticker names must be unique.');
end
[present, columns] = ismember(canonical, rawTickers);
if ~all(present) || any(~ismember(rawTickers, [canonical, "USDCNH Index"]))
    error('weekly:Tickers', 'Expected the 14 supported tickers, optionally with USDCNH.');
end
opts = delimitedTextImportOptions('NumVariables', numel(rawTickers)+1);
opts.DataLines = [4 Inf];
opts.Delimiter = ',';
opts.VariableNames = ["Date", "Series" + (1:numel(rawTickers))];
opts.VariableTypes = ["datetime", repmat("double",1,numel(rawTickers))];
opts = setvaropts(opts, 'Date', 'InputFormat', 'dd/MM/yyyy');
opts.ExtraColumnsRule = 'error';
opts.ImportErrorRule = 'error';
opts.EmptyLineRule = 'skip';
raw = readtable(sourceFile, opts);
rawDates = raw.Date;
rawLevels = raw{:, columns+1};
if isempty(rawDates) || any(isnat(rawDates)) || any(diff(rawDates) <= days(0))
    error('weekly:Dates', 'Dates must be valid, unique, and strictly increasing.');
end
if any(ismember(weekday(rawDates), [1 7]))
    error('weekly:Dates', 'Expected a weekday-only source calendar.');
end
if any(isinf(rawLevels(:))) || any(rawLevels(isfinite(rawLevels)) <= 0)
    error('weekly:Levels', 'Every observed level must be finite and strictly positive.');
end
keep = rawDates >= startDate;
rawDates = rawDates(keep);
rawLevels = rawLevels(keep,:);
if isempty(rawDates), error('weekly:EmptySample', 'No observations after start_date.'); end

% Drop partial initial and final calendar weeks. MATLAB weekday: Sun=1.
firstMonday = rawDates(1) + caldays(mod(2-weekday(rawDates(1)), 7));
firstFriday = firstMonday + caldays(4);
lastFriday = rawDates(end) - caldays(mod(weekday(rawDates(end))-6, 7));
levelDates = (firstFriday:calweeks(1):lastFriday)';
if numel(levelDates) < 2
    error('weekly:ShortSample', 'At least two complete weeks of levels are required.');
end
nLevels = numel(levelDates);
m = numel(canonical);
levels = nan(nLevels,m);
sourceDates = NaT(nLevels,m);
ageDays = nan(nLevels,m);
closure = false(nLevels,m);
for j = 1:m
    observed = find(isfinite(rawLevels(:,j)));
    if isempty(observed)
        error('weekly:MissingEndpoint', 'No observed levels for %s.', canonical(j));
    end
    bins = discretize(datenum(levelDates), [datenum(rawDates(observed)); Inf]);
    if any(isnan(bins))
        bad = find(isnan(bins),1);
        error('weekly:MissingEndpoint', 'No past level for %s at %s.', ...
            canonical(j), string(levelDates(bad),'yyyy-MM-dd'));
    end
    indices = observed(bins);
    levels(:,j) = rawLevels(indices,j);
    sourceDates(:,j) = rawDates(indices);
    ageDays(:,j) = days(levelDates-sourceDates(:,j));
    closure(:,j) = canonical(j) == "NKYTR Index" & ...
        levelDates == datetime(2019,5,3) & ...
        sourceDates(:,j) == datetime(2019,4,26) & ageDays(:,j) == 7;
    bad = find(ageDays(:,j) > 4 & ~closure(:,j),1);
    if ~isempty(bad)
        error('weekly:UnexplainedGap', ...
            '%s endpoint %s is %g days old (source %s); investigate before estimation.', ...
            canonical(j), string(levelDates(bad),'yyyy-MM-dd'), ageDays(bad,j), ...
            string(sourceDates(bad,j),'yyyy-MM-dd'));
    end
end

returns = 100*diff(log(levels),1,1);
closureReturn = closure(1:end-1,:) | closure(2:end,:);
observationMask = true(size(returns));
if policy == "mask", observationMask(closureReturn) = false; end
[pairI,pairJ] = find(tril(true(m),-1));
fid = fopen(sourceFile,'rb');
if fid < 0, error('weekly:SourceFile', 'Cannot read source for SHA-256.'); end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
sourceHash = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
clear cleanup

panel = struct();
panel.dates = levelDates(2:end);
panel.level_dates = levelDates;
panel.tickers = canonical;
panel.labels = ["AUD/USD", "Brent total return", "Gold total return", ...
    "EUR/USD", "GBP/USD", "Japan government bonds", "Euro-area Treasuries", ...
    "UK gilts", "US Treasuries", "Japan equities", "US equities", ...
    "Euro-area equities", "UK equities", "USD/JPY"];
panel.levels = levels;
panel.returns = returns;
panel.source_dates = sourceDates;
panel.age_days = ageDays;
panel.carried = ageDays > 0;
panel.closure_level_mask = closure;
panel.closure_return_mask = closureReturn;
panel.observation_mask = observationMask;
panel.pair_i = pairI;
panel.pair_j = pairJ;
panel.source_file = sourceFile;
panel.source_hash = sourceHash;
panel.closure_policy = char(policy);
panel.currency_basis = 'native quotes; USDJPY retains its original direction';
panel.frequency = 'weekly, Friday as-of, percent log returns';
panel.raw_dates = rawDates;
panel.raw_levels = rawLevels;
panel.start_date_requested = startDate;
panel.excluded_final_rows = sum(rawDates > lastFriday);
panel.closure_reference = 'https://indexes.nikkei.co.jp/nkave/archives/news/20190423E_1.pdf';
end

function limits = correlation_panel_limits(varargin)
%CORRELATION_PANEL_LIMITS Return data bounds padded by one correlation point.
% Finite values define the panel range; NaNs in uncertainty bands are ignored.
% The fixed 0.01 padding is one percentage point, and bounds are not clamped.
values = [];
for k = 1:nargin
    x = varargin{k};
    if ~isnumeric(x) || ~isreal(x)
        error('dsc:PlotLimits', ...
            'Correlation panel limits require numeric real inputs.');
    end
    values = [values; double(x(:))]; %#ok<AGROW>
end
values = values(isfinite(values));
if isempty(values)
    error('dsc:PlotLimits', ...
        'Correlation panel limits require at least one finite value.');
end
limits = [min(values)-0.01, max(values)+0.01];
end

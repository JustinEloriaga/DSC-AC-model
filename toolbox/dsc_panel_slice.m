function sliced=dsc_panel_slice(panel,first,last)
%DSC_PANEL_SLICE Keep only a contiguous interval of weekly model observations.
if first<1||last>numel(panel.dates)||first>last|| ...
        first~=floor(first)||last~=floor(last)
    error('dsc:PanelSlice','Invalid weekly return interval.');
end
sliced=struct('dates',panel.dates(first:last),'returns',panel.returns(first:last,:), ...
    'observation_mask',panel.observation_mask(first:last,:),'tickers',panel.tickers, ...
    'pair_i',panel.pair_i,'pair_j',panel.pair_j,'source_hash',panel.source_hash);
end

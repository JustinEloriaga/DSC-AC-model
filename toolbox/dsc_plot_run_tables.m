function pages = dsc_plot_run_tables(pdf_path,run_dir,file_prefix,render)
%DSC_PLOT_RUN_TABLES Append parameter and return AR(1) tables to a path PDF.
% The optional render=false is for inspecting the prepared page contents
% without authoring a PDF. It is internal to the plotting workflow.
if nargin<4, render=true; end
manifest=struct();
manifest_path=fullfile(run_dir,[char(file_prefix) 'run_manifest.json']);
if isfile(manifest_path), manifest=jsondecode(fileread(manifest_path)); end
parameter_section=read_parameter_section(manifest,run_dir,file_prefix);
ar_section=read_ar_section(manifest);
pages=paginate_sections(parameter_section,ar_section);
if ~render, return; end
for k=1:numel(pages)
    draw_page(pdf_path,pages(k),k,numel(pages));
end
end

function section=read_parameter_section(manifest,root,prefix)
section=empty_section('Parameter priors and posteriors');
path=fullfile(root,[char(prefix) 'parameters.mat']);
if isfield(manifest,'outputs')&&isstruct(manifest.outputs)&& ...
        isfield(manifest.outputs,'parameters')
    value=manifest.outputs.parameters;
    if ~((ischar(value)&&isrow(value)&&~isempty(value))|| ...
            (isstring(value)&&isscalar(value)&&strlength(value)>0))
        error('DSC:ParameterTable','The manifest parameter file must be scalar text.');
    end
    path=char(value);
    if ~java.io.File(path).isAbsolute(), path=fullfile(root,path); end
end
if ~isfile(path)
    section.note={'The parameter prior/posterior table is unavailable in this saved run.'};
    return
end
loaded=load(path,'parameters');
if ~isfield(loaded,'parameters')||~isstruct(loaded.parameters)|| ...
        ~isscalar(loaded.parameters)||~isfield(loaded.parameters,'prior_posterior_summary')
    section.note={'The parameter prior/posterior table is unavailable in this saved run.'};
    return
end
parameters=loaded.parameters;
inference=table_inference(parameters,manifest);
summary=parameters.prior_posterior_summary;
required={'Parameter','Series','PriorDistribution','PriorMean','PriorStd','PosteriorMean','PosteriorStd'};
if ~istable(summary)||~all(ismember(required,summary.Properties.VariableNames))||isempty(summary)
    error('DSC:ParameterTable','The saved parameter comparison must be a nonempty table with prior/posterior moments.');
end
section.headers={'Param.','Series or pair','Prior','Prior mean','Prior SD','Post. mean','Post. SD'};
section.widths=[.075 .225 .10 .11 .11 .19 .19];
section.numeric_columns=4:7;
section.rows=cell(height(summary),numel(required));
for r=1:height(summary)
    section.rows{r,1}=char(string(summary.Parameter(r)));
    section.rows{r,2}=char(erase(string(summary.Series(r))," Index"));
    distribution=string(summary.PriorDistribution(r));
    if startsWith(distribution,'InverseGamma'), distribution="Inv-gamma"; end
    section.rows{r,3}=char(distribution);
    for c=4:7
        values=summary.(required{c});
        if ~isnumeric(values)||~isreal(values)||size(values,2)~=1
            error('DSC:ParameterTable','Each prior/posterior moment must be a real numeric column.');
        end
        if c>=6
            prior_values=summary.(required{c-2});
            if ~isnumeric(prior_values)||~isreal(prior_values)||size(prior_values,2)~=1
                error('DSC:ParameterTable','Each prior moment must be a real numeric column.');
            end
            section.rows{r,c}=relative_posterior_text(values(r),prior_values(r));
        else
            section.rows{r,c}=number_text(values(r),false);
        end
    end
end
section.note={ ...
    'Posterior moments pool retained parameter draws across chains; SD is marginal standard deviation.'};
conditional=isfield(parameters,'storage')&&strcmp(parameters.storage,'conditional_parameter_summary_v1');
conditional=conditional||(isfield(inference,'inference_mode')&& ...
    strcmp(inference.inference_mode,'fixed_parameter_smoothing'));
if isfield(manifest,'config')&&isstruct(manifest.config)&& ...
        isfield(manifest.config,'EstimateParameters')
    conditional=conditional||isequal(manifest.config.EstimateParameters,false);
end
if conditional
    section.title='Conditional parameter summary';
    section.headers(6:7)={'Cond. mean','Cond. SD'};
    smoothing='';
    if isfield(inference,'parameter_smoothing'), smoothing=char(inference.parameter_smoothing); end
    if strcmp(smoothing,'mean')
        section.note={ ...
            'Conditional summary of retained smoothing draws; variance parameters are held fixed.', ...
            'Zero SD for sig2h / sig2r reflects fixed inputs, not estimated posterior uncertainty.'};
    elseif strcmp(smoothing,'draws')
        section.note={ ...
            'Conditional summary of smoothing draws using retained parameter draws from the saved fit.'};
    else
        section.note={'Conditional summary of retained smoothing draws from the saved parameter inputs.'};
    end
end
if all(isfield(parameters,{'estimation_start','estimation_end'}))
    first=date_text(parameters.estimation_start); last=date_text(parameters.estimation_end);
elseif all(isfield(inference,{'parameter_estimation_start','parameter_estimation_end'}))
    first=date_text(inference.parameter_estimation_start);
    last=date_text(inference.parameter_estimation_end);
else
    first=''; last='';
end
if ~isempty(first)&&~isempty(last)
    section.note{end+1}=sprintf('Parameter estimation: %s to %s.',first,last);
end
if conditional&&isfield(inference,'smoothing_end')
    if isfield(inference,'smoothing_start')
        first=date_text(inference.smoothing_start);
    elseif isfield(inference,'parameter_estimation_start')
        first=date_text(inference.parameter_estimation_start);
    else
        first='';
    end
    last=date_text(inference.smoothing_end);
    if ~isempty(first)&&~isempty(last)
        section.note{end+1}=sprintf('Conditional smoothing: %s to %s.',first,last);
    end
end
section.note=[section.note { ...
    'Parentheses: % change = 100*(posterior/prior - 1); N/A when undefined.', ...
    'Prior moments are analytic under empirical-Bayes priors calibrated from the data.', ...
    'mu: constant return mean; sig2h: log-variance-state innovation variance.', ...
    'sig2r: matrix-log correlation-coordinate innovation variance.'}];
end

function inference=table_inference(parameters,manifest)
inference=struct();
% Scope belongs to the saved private report. Full inference JSON is only a
% compatibility fallback for historical reports that predate those fields.
if isfield(manifest,'inference')&&isstruct(manifest.inference)&&isscalar(manifest.inference)
    inference=manifest.inference;
end
if isfield(manifest,'config')&&isstruct(manifest.config)&& ...
        isfield(manifest.config,'ParameterSmoothing')
    inference.parameter_smoothing=manifest.config.ParameterSmoothing;
end
if isfield(parameters,'inference')
    value=parameters.inference;
    if ~isstruct(value)||~isscalar(value)
        error('DSC:ParameterTable','The private inference metadata must be a scalar structure.');
    end
    fields=fieldnames(value);
    for k=1:numel(fields), inference.(fields{k})=value.(fields{k}); end
end
for key={'inference_mode','parameter_smoothing'}
    name=key{1};
    if isfield(parameters,name), inference.(name)=parameters.(name); end
end
if isfield(parameters,'dates')&&isdatetime(parameters.dates)&&~isempty(parameters.dates)
    if ~isfield(inference,'smoothing_start'), inference.smoothing_start=parameters.dates(1); end
    if ~isfield(inference,'smoothing_end'), inference.smoothing_end=parameters.dates(end); end
end
end

function value=date_text(value)
if isdatetime(value)&&isscalar(value)&&~isnat(value)
    value=char(string(value,'yyyy-MM-dd'));
elseif (ischar(value)&&isrow(value))||(isstring(value)&&isscalar(value))
    value=char(value);
    if isempty(regexp(value,'^\d{4}-\d{2}-\d{2}$','once')), value=''; end
else
    value='';
end
end

function section=read_ar_section(manifest)
section=empty_section('Return AR(1): initial prior calibration');
if ~isfield(manifest,'return_ar1')||~isstruct(manifest.return_ar1)|| ...
        ~isscalar(manifest.return_ar1)||~isfield(manifest.return_ar1,'rows')|| ...
        isempty(manifest.return_ar1.rows)
    section.note={'The initial prior-calibration return AR(1) regression is unavailable in this saved run.'};
    return
end
report=manifest.return_ar1; rows=report.rows;
required={'ticker','p_value'};
if ~isstruct(rows)||~all(isfield(rows,required))
    error('DSC:ReturnAR1Table','The saved return AR(1) rows are incomplete.');
end
section.headers={'Series','p-value'};
section.widths=[.7 .3];
section.numeric_columns=2;
section.rows=cell(numel(rows),numel(section.headers));
for r=1:numel(rows)
    section.rows{r,1}=char(erase(string(rows(r).ticker)," Index"));
    value=rows(r).p_value;
    if isempty(value)||(isnumeric(value)&&isscalar(value)&&isnan(value))
        section.rows{r,2}='--';
    elseif ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)||value<0||value>1
        error('DSC:ReturnAR1Table','AR(1) p-values must be finite scalars in [0,1] or unavailable.');
    else
        section.rows{r,2}=number_text(value,true);
    end
end
if isfield(report,'sample_start')&&isfield(report,'sample_end')
    section.title=sprintf('Return AR(1): initial prior calibration | %s to %s',report.sample_start,report.sample_end);
elseif isfield(manifest,'inference')&&isstruct(manifest.inference)&& ...
        all(isfield(manifest.inference,{'calibration_start','calibration_end'}))
    first=date_text(manifest.inference.calibration_start);
    last=date_text(manifest.inference.calibration_end);
    if ~isempty(first)&&~isempty(last)
        section.title=sprintf('Return AR(1): initial prior calibration | %s to %s',first,last);
    end
end
section.note={ ...
    'x_t = a + rho*x_(t-1) + e_t; H0: rho = 0, two-sided standard-normal p-values.', ...
    'Weekly percent log returns; OLS with Newey-West HAC (Bartlett), N/(N-2) adjustment.', ...
    'Initial prior-calibration sample; 103 consecutive-week pairs with 104 complete returns.', ...
    'Failure to reject rho = 0 does not establish i.i.d. returns.'};
end

function section=empty_section(title_text)
section=struct('title',title_text,'headers',{{}},'widths',[], ...
    'numeric_columns',[],'rows',{cell(0,0)},'note',{{}});
end

function pages=paginate_sections(parameters,ar)
% A four-series run fits its 14 parameter rows and four AR rows on one page.
% Larger tables retain the same type size and row spacing across pages.
pages=struct('sections',{}); budget=.80; row_height=.023;
sources={parameters,ar}; used=0; page=struct('sections',{{}});
for s=1:numel(sources)
    source=sources{s}; offset=0; total=size(source.rows,1);
    if total==0
        cost=.045+.017*numel(source.note);
        if used+cost>budget&&~isempty(page.sections)
            pages(end+1)=page; page=struct('sections',{{}}); used=0; %#ok<AGROW>
        end
        page.sections{end+1}=source; used=used+cost;
        continue
    end
    while offset<total
        overhead=.075+.017*numel(source.note);
        capacity=floor((budget-used-overhead)/row_height);
        remaining=total-offset;
        if capacity<1||(s==2&&remaining>capacity&&~isempty(page.sections))
            pages(end+1)=page; page=struct('sections',{{}}); used=0; %#ok<AGROW>
            continue
        end
        count=min(remaining,capacity); part=source;
        part.rows=source.rows(offset+(1:count),:);
        if count<total
            part.title=sprintf('%s | Rows %d-%d of %d',source.title,offset+1,offset+count,total);
        end
        page.sections{end+1}=part; used=used+overhead+count*row_height;
        offset=offset+count;
        if offset<total
            pages(end+1)=page; page=struct('sections',{{}}); used=0; %#ok<AGROW>
        end
    end
end
if ~isempty(page.sections), pages(end+1)=page; end
end

function text_value=number_text(value,is_p)
if isempty(value)||~isnumeric(value)||~isreal(value)||~isscalar(value)||isnan(value)
    text_value='--';
elseif isinf(value)
    if value>0, text_value='Inf'; else, text_value='-Inf'; end
elseif is_p&&value==0
    text_value='<1e-300';
elseif value==0
    text_value='0';
elseif abs(value)<1e-3||abs(value)>=1e4
    text_value=sprintf('%.3e',value);
else
    text_value=sprintf('%.4g',value);
end
end

function text_value=relative_posterior_text(posterior,prior)
base=number_text(posterior,false);
if ~isnumeric(posterior)||~isreal(posterior)||~isscalar(posterior)||~isfinite(posterior)|| ...
        ~isnumeric(prior)||~isreal(prior)||~isscalar(prior)||~isfinite(prior)||prior==0
    text_value=[base ' (N/A)'];
    return
end
change=100*(posterior/prior-1);
if ~isfinite(change)
    text_value=[base ' (N/A)'];
else
    text_value=sprintf('%s (%+.1f%%)',base,change);
end
end

function draw_page(path,page,index,count)
fig=figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1200 1500]);
cleanup=onCleanup(@()close(fig));
ax=axes(fig,'Position',[.045 .035 .91 .93],'XLim',[0 1],'YLim',[0 1],'Visible','off');
hold(ax,'on'); navy=[.09 .19 .30]; muted=[.34 .39 .44];
% Relative font sizes preserve text-to-page proportions when the PDF
% exporter uses the display's device-pixel scaling (including Retina).
font_unit=1/700;
text(ax,0,.985,sprintf('Convergence: parameters and return AR(1) | Supplement %d of %d',index,count), ...
    'FontUnits','normalized','FontSize',14*font_unit,'FontWeight','bold', ...
    'Color',navy,'Interpreter','none','VerticalAlignment','top');
y=.91; row_height=.023;
for s=1:numel(page.sections)
    section=page.sections{s};
    text(ax,0,y,section.title,'FontUnits','normalized','FontSize',11*font_unit, ...
        'FontWeight','bold','Color',navy, ...
        'Interpreter','none','VerticalAlignment','top');
    y=y-.035;
    if ~isempty(section.rows)
        edges=[0 cumsum(section.widths)];
        rectangle(ax,'Position',[0 y-.031 1 .031],'FaceColor',navy,'EdgeColor','none');
        for c=1:numel(section.headers)
            [x,align]=cell_position(edges,c,section.numeric_columns);
            text(ax,x,y-.0155,section.headers{c},'FontUnits','normalized', ...
                'FontSize',9*font_unit,'FontWeight','bold', ...
                'Color','w','Interpreter','none','HorizontalAlignment',align,'VerticalAlignment','middle');
        end
        y=y-.037;
        for r=1:size(section.rows,1)
            if mod(r,2)==1
                rectangle(ax,'Position',[0 y-row_height 1 row_height], ...
                    'FaceColor',[.945 .960 .975],'EdgeColor','none');
            end
            for c=1:numel(section.headers)
                [x,align]=cell_position(edges,c,section.numeric_columns);
                text(ax,x,y-row_height/2,section.rows{r,c},'FontUnits','normalized', ...
                    'FontSize',9.5*font_unit,'Color',navy, ...
                    'Interpreter','none','HorizontalAlignment',align,'VerticalAlignment','middle');
            end
            y=y-row_height;
        end
        y=y-.008;
    end
    for n=1:numel(section.note)
        text(ax,0,y,section.note{n},'FontUnits','normalized', ...
            'FontSize',8.5*font_unit,'Color',muted, ...
            'Interpreter','none','VerticalAlignment','top');
        y=y-.017;
    end
    y=y-.012;
end
exportgraphics(fig,path,'ContentType','vector','Append',true);
end

function [x,align]=cell_position(edges,column,numeric_columns)
if ismember(column,numeric_columns)
    x=edges(column+1)-.009; align='right';
else
    x=edges(column)+.009; align='left';
end
end

function result = run_weekly_model1_report(varargin)
%RUN_WEEKLY_MODEL1_REPORT Run the weekly report on variables 10-13 only.
% Model 1 uses NKYTR, SPXT, SX5T, and TUKXG in the source file order.
root = fileparts(mfilename('fullpath'));
model1Variables = ["NKYTR Index", "SPXT Index", "SX5T Index", "TUKXG Index"];
args = varargin;
if ~has_name_argument(args, "OutputRoot")
    args = [args, {'OutputRoot', fullfile(root,'outputs','weekly_research_model1')}];
end
if ~has_name_argument(args, "BuildReport")
    args = [args, {'BuildReport', false}];
end
if ~has_name_argument(args, "VerifyReport")
    args = [args, {'VerifyReport', false}];
end
if ~has_name_argument(args, "ParameterSmoothing")
    args = [args, {'ParameterSmoothing', "draws"}];
end
if has_name_argument(args, "VariableNames")
    error('DSC:Model1Options', 'run_weekly_model1_report fixes VariableNames to the model1 ticker set.');
end
result = run_weekly_model_report(args{:}, 'VariableNames', model1Variables);
end

function tf = has_name_argument(args, name)
tf = false;
for k = 1:2:numel(args)
    if (ischar(args{k}) || isstring(args{k})) && strcmpi(string(args{k}), name)
        tf = true;
        return
    end
end
end

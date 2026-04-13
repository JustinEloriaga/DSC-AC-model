function set_T = seq_date(eval_T0, eval_T1)
% generate list of date sequences (in cell str) in our format
% eval_T0: start
% eval_T1: end

% matlab format
T0_vec = datevec(eval_T0, 'mm/dd/yyyy');
T1_vec = datevec(eval_T1, 'mm/dd/yyyy');

% cell string
set_T = cell(1,1);


% at t=1
% set_T{1} = eval_T0;

% rest
cur_T = T0_vec;

while ( any(cur_T ~= T1_vec) )
       
    % Make it in our format
    temp_str = datestr(cur_T,'mm/dd/yyyy');
    temp_str(strfind(temp_str(1:5), '0')) = [];
    temp_str(strfind(temp_str,'/')) = '-';
    
    % Stack in our vector
    set_T = [set_T; temp_str];
    
    % Move forward by one quarter (+3m)
    if (cur_T(2) + 3) == 15 %year
        cur_T(1) = cur_T(1) + 1; %move year
        cur_T(2) = 3; %reset month
    else
        cur_T(2) = cur_T(2) + 3; %move month
    end
    
end

% ---
% Store last period
% Make it in our format
temp_str = datestr(cur_T,'mm/dd/yyyy');
temp_str(strfind(temp_str(1:5), '0')) = [];
temp_str(strfind(temp_str,'/')) = '-';

% Stack in our vector
set_T = [set_T; temp_str];


% ---
% kill the first empty cell
set_T(1) = [];


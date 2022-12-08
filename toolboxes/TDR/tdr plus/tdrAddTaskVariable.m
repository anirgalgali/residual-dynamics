function data_add = tdrAddTaskVariable(data,addpars)
% tdrAddTaskVariable Add task variable to data structure
%
% Inputs:
%  data
%  addpars.task_variable
%
% Output:
%  data_add
%
% data_add = tdrAddTaskVariable(data,addpars);


% Initialize
data_add = data;

if isfield(data,'task_variable') && ~isempty(data.task_variable)
    % New task variables
    variable_name = fieldnames(addpars.task_variable);
    nvr = length(variable_name);
    
    % Loop over new variables
    for ivr = 1:nvr
        % Add task variable
        data_add.task_variable.(variable_name{ivr}) = ...
            addpars.task_variable.(variable_name{ivr});
        
        % Initialize task index
        ntr = length(addpars.task_variable.(variable_name{ivr}));
        data_add.task_index.(variable_name{ivr}) = zeros(ntr,1);
        
        % Determine task indeces
        variable_unique = unique(addpars.task_variable.(variable_name{ivr}));
        nvl = length(variable_unique);
        for ivl = 1:nvl
            data_add.task_index.(variable_name{ivr})(addpars.task_variable.(variable_name{ivr}) == variable_unique(ivl)) = ivl;
        end
        
        % The variable-index pairs
        data_add.task_variable_index.(variable_name{ivr}) = [variable_unique (1:nvl)'];
    end
end


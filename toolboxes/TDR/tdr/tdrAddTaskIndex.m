function data = tdrAddTaskIndex(data)
% tdrAddTaskIndex add task indeces based on task variables
%
% Input: data.task_variable
%
% data = tdrAddTaskIndex(data)

variable_name = fieldnames(data.task_variable);
nvr = length(variable_name);
for ivr = 1:nvr
    [~,~,data.task_index.(variable_name{ivr})(:,1)] = ...
        unique(data.task_variable.(variable_name{ivr}));
end

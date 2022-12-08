function data_boot = tdrBootstrapResample(data,bootpars,jboot)

% Initialize
data_boot = data;

% Check if data is sequentially or simultaneously recorded
if isfield(data,'unit') && ~isempty(data.unit)
    
    %--- Sequential recordings ---
    
    % Number of units
    nun = length(data.unit);
    
    % Loop over units
    for iun = 1:nun
        
        % Trials to use
        itrial = bootpars.iboot{iun}(jboot,:);
        
        % Resample responses
        data_boot.unit(iun).response = ...
            data.unit(iun).response(itrial,:);
        
        % Resample task variable
        variable_name = fieldnames(data.unit(iun).task_variable);
        for ivr = 1:length(variable_name)
            data_boot.unit(iun).task_variable.(variable_name{ivr}) = ...
                data.unit(iun).task_variable.(variable_name{ivr})(itrial);
        end
        
    end
    
else
    
    
    %--- Simultaneous recordings ---
    
    % Trials to use
    itrial = bootpars.iboot(jboot,:);
    
    % Resample responses
    data_boot.response = data.response(:,:,itrial);
    
    % Resample task variables
    variable_name = fieldnames(data.task_variable);
    for ivr = 1:length(variable_name)
        data_boot.task_variable.(variable_name{ivr}) = ...
            data.task_variable.(variable_name{ivr})(itrial);
    end
    
    % Resample task indeces
    index_name = fieldnames(data.task_index);
    for ivr = 1:length(index_name)
        data_boot.task_index.(index_name{ivr}) = ...
            data.task_index.(index_name{ivr})(itrial);
    end
    
    % Resample ntrial
    data_boot.n_trial = data.n_trial(itrial);
    
end

function [data_avg,data_sub] = tdrAverageCondition(data,task_index)
% tdrAverageCondition compute condition averaged responses
%
% Inputs:
%  data: population response (sequential or simultaneous, see 'help tdr')
%  task_index: indeces of task variables to average.
%
%     e.g. task_index.stim_dir = [1 2 3]
%          task_index.targ_dir = [2 2 2]
%
%     would result in three averages, over the following sets of trials:
%      (1) stim_dir==unique(stim_dir)(1) & targ_dir==unique(stim_dir)(2) 
%      (2) stim_dir==unique(stim_dir)(2) & targ_dir==unique(stim_dir)(2) 
%      (3) stim_dir==unique(stim_dir)(3) & targ_dir==unique(stim_dir)(2) 
%
%     Setting an index to zero corresponds to averaging out the
%     corresponding task variable.
%
% Outputs:
%  data_avg: condition averages. Same format as simultaneous trial-by-trial
%     responses (see 'help tdr').
%  data_sub: condition-average-subtracted trial-by-trial responses. Here the 
%     responses in data_avg are subtracted from corresponding trials in data.
%     Only meaningful if each trial contributes to at most 1 condition. 
%
% [data_avg,data_sub] = tdrAverageCondition(data,task_index)



% All the task indeces
variable_name = fieldnames(task_index);
nvr = length(variable_name);

% Number of conditions
ncd = length(task_index.(variable_name{1}));

% Number of time samples
npt = length(data.time);
    
% Check if data is sequentially or simultaneously recorded
if isfield(data,'unit') && ~isempty(data.unit)
    
    
    %--- Sequential recordings ---    
    % Initialize
    data_sub = data;
    
    % Number of units
    nun = length(data.unit);
    
    % Initialize
    response_avg = NaN(nun,npt,ncd); 
    n_trial = zeros(nun,ncd); 
    variable_condition_value = NaN(nvr,ncd,nun);
    
    for iun = 1:nun
        
        % Number of trials
        ntr = size(data.unit(iun).response,1);
        
        % Initialize trials to average
        jj_cnd = false(ncd,ntr);
        
        % Unique variable names
        variable_unique = cell(1,nvr);
        for ivr = 1:nvr
            variable_unique{ivr} = unique(data.unit(iun).task_variable.(variable_name{ivr}));
        end           
        
        % Loop over conditions
        for icd = 1:ncd
            
            % Trials to average
            jj_var = zeros(nvr,ntr);
            
            % Loop over variables
            for ivr = 1:nvr
                
                if task_index.(variable_name{ivr})(icd)==0
                    % Take all trials
                    jj_var(ivr,:) = true(ntr,1);
                else
                    % Find matching trials
                    jj_var(ivr,:) = ...
                        data.unit(iun).task_variable.(variable_name{ivr}) == ...
                        variable_unique{ivr}(task_index.(variable_name{ivr})(icd));
                    
                    % Keep variable value
                    variable_condition_value(ivr,icd,iun) = ...
                        variable_unique{ivr}(task_index.(variable_name{ivr})(icd));
                end
            end
            
            % Fullfill constraints on all indeces
            jj_cnd(icd,:) = prod(jj_var,1);
            
            % Number of trials for this condition
            n_trial(iun,icd) = sum(jj_cnd(icd,:));
            
            % Condition average
            if n_trial(iun,icd) > 0
                response_avg(iun,:,icd) = mean(data.unit(iun).response(jj_cnd(icd,:),:));
            end
            
            % Average-subtracted responses
            if n_trial(iun,icd) > 0
                data_sub.unit(iun).response(jj_cnd(icd,:),:) = ...
                    data.unit(iun).response(jj_cnd(icd,:),:) - ...
                    repmat(response_avg(iun,:,icd),[n_trial(iun,icd) 1]);
            end
            
        end
    end
    
    % Keep average values of the task variables
    clear task_variable
    for ivr = 1:nvr
       task_variable.(variable_name{ivr}) = squeeze(nanmean(variable_condition_value(ivr,:,:),3))'; 
    end
    
    % The task events
    task_event = [];
    
    % Keep dimension names
    dimension = {data.unit(:).dimension}';
    
else
    
    
    %--- Simultaneous recordings ---    
    % Dimensions
    [nun npt ntr] = size(data.response);
    
    % Initialize
    data_sub = data; 
    response_avg = NaN(nun,npt,ncd); 
    n_trial = zeros(1,ncd); 
    clear task_variable
    for ivr = 1:nvr
        task_variable.(variable_name{ivr}) = NaN(ncd,1);
    end
    
    % Initialize trials to average
    jj_cnd = false(ncd,ntr);
    
    % Unique variable names
    variable_unique = cell(1,nvr);
    for ivr = 1:nvr
        variable_unique{ivr} = unique(data.task_variable.(variable_name{ivr}));
    end
    
    % Initialize task events
    if isfield(data,'task_event') && ~isempty(data.task_event)
        event_flag = 1;
        event_names = fieldnames(data.task_event);
        nev = length(event_names);
        for iev = 1:nev
            task_event.(event_names{iev}) = NaN(ncd,1);
        end
    else
        event_flag = 0;
    end
        
    % Loop over conditions
    for icd = 1:ncd
        
        % Trials to average
        jj_var = zeros(nvr,ntr);
        
        % Loop over variables
        for ivr = 1:nvr
            
            if task_index.(variable_name{ivr})(icd)==0
                % Take all trials
                jj_var(ivr,:) = true(ntr,1);
            else
                % Find matching trials
                jj_var(ivr,:) = ...
                    data.task_variable.(variable_name{ivr}) == ...
                    variable_unique{ivr}(task_index.(variable_name{ivr})(icd));
                
                % Keep variable value
                task_variable.(variable_name{ivr})(icd) = ...
                    variable_unique{ivr}(task_index.(variable_name{ivr})(icd));
            end
        end
        
        % Fullfill constraints on all indeces
        jj_cnd(icd,:) = prod(jj_var,1);
        
        % Number of trials for this condition
        n_trial(1,icd) = sum(jj_cnd(icd,:));
        
        % Condition average
        if n_trial(1,icd) > 0
            response_avg(:,:,icd) = mean(data.response(:,:,jj_cnd(icd,:)),3);
        end
        
        % Average-subtracted response
        if n_trial(1,icd) > 0
            data_sub.response(:,:,jj_cnd(icd,:)) = ...
                data.response(:,:,jj_cnd(icd,:)) - repmat(response_avg(:,:,icd),[1 1 n_trial(1,icd)]);
        end
        
        % Average the task events
        if event_flag
            for iev = 1:nev
                if n_trial(1,icd) > 0
                    task_event.(event_names{iev})(icd) = mean(data.task_event.(event_names{iev})(jj_cnd(icd,:)));
                end
            end
        else
            task_event = [];
        end
    end
    
    % Keep dimension names
    if(isfield(data,'dimension'))
    
        dimension = data.dimension;
        
    end
end

% The variable-index pairs
for ivr = 1:nvr
    
    % The unique indeces
    [~,unique_ind] = unique(task_index.(variable_name{ivr}));
    
    % The pairs
    task_variable_index.(variable_name{ivr}) = [...
        task_variable.(variable_name{ivr})(unique_ind) task_index.(variable_name{ivr})(unique_ind)];
    
end

% Keep what you need
data_avg.response = response_avg;
data_avg.task_variable = task_variable;
data_avg.task_index = task_index;
data_avg.task_variable_index = task_variable_index;
data_avg.n_trial = n_trial;
data_avg.time = data.time;

if(isfield(data,'time_rel'))
    data_avg.time_rel = data.time_rel;
end

if(isfield(data,'time_iev'))
    data_avg.time_iev = data.time_iev;
end

if(isfield(data,task_event))
   data_avg.task_event = taks_event; 
end
if(isfield(data,'dimension'))
    data_avg.dimension = dimension;
end


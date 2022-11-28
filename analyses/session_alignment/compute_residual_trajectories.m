function [data_out, excluded_trials_idxs] = compute_residual_trajectories(data, condition_vars, exclude_trials)
%{ Computes residuals based on defined conditions
%  Input
%  - data (cell array of structs) - of size n_experiments x 1. Each element 
%    of the cell array is a struct that contains the single trial data
%    belonging to a single experiment
%  - condition_vars (cell array) - specifies the task variables based on 
%    which trial are sorted prior to computing residuals
%  - exclude_trials (cell array) - specifies the task variables and their
%    corresponding settings which can be used to filter out trials from the
%    dataset
%  Output
%  - data_out (cell array of structs) - similar to data, but contains the
%    aligned single-trial responses.
%  - excluded_trials_idxs (cell array of logicals) - contains the indices
%    of the excluded trials in each experiment
%
%  Author: Aniruddh Galgali
%
assert(iscell(data),'variable: data - must be a cell array');
assert(iscell(condition_vars),'variable: task_conditions - must be a cell array');

num_sessions = length(data);

data_out = cell(1,num_sessions);
excluded_trials_idxs = cell(1,num_sessions);

for i_session = 1:num_sessions
    
    unique_time_labels = unique(data{i_session}.time_iev);
    num_alignments = length(unique_time_labels);
    assert(length(condition_vars) == num_alignments, 'task-conditions have to be specified separately for each alignment');
    
    if(~isempty(exclude_trials))
       num_exclusion_conditions = length(exclude_trials);
       trial_to_exclude = true(size(data{i_session}.response,3),1);
       
       for iexc = 1: num_exclusion_conditions
           cond_name_of_exclusion = exclude_trials{iexc}{1};
           cond_idx_of_exclusion = exclude_trials{iexc}{2};
           trial_to_exclude = trial_to_exclude & (data{i_session}.task_index.(cond_name_of_exclusion) == cond_idx_of_exclusion);
       end
            
      getpars.trial_get = ~trial_to_exclude;
      data{i_session} = tdrGetTrials(data{i_session}, getpars);  
      excluded_trials_idxs{i_session} =  trial_to_exclude;
    
    end
    
    data_out{i_session} = data{i_session};
    
    for iev = 1:num_alignments
        
        time_mask = data{i_session}.time_iev == unique_time_labels(iev);
        time_ev = data{i_session}.time(time_mask);
        time_window.tStart = time_ev(1);
        time_window.tEnd = time_ev(end);
        D_ev = extractAnalysisTimeWindow(data{i_session},time_window);

        if(~isempty(condition_vars{iev}))
            
            [D_cond,T] = sort_trials_by_condition(D_ev,condition_vars{iev});
            cond_var_names = T.Properties.VariableNames;
            
            for i_cond = 1: length(D_cond)
                
                D_cond_mean = mean(D_cond{i_cond}.response,3);
                
                valid_inds = true(size(D_ev.response,3),1);
                for icond_vars = 1:length(cond_var_names)
                    
                    if(iscellstr(D_ev.task_variable.(cond_var_names{icond_vars})))
                        valid_inds = valid_inds & cell2mat(cellfun(@(x) strcmp(x, T.(cond_var_names{icond_vars}){i_cond}), ...
                            D_ev.task_variable.(cond_var_names{icond_vars}),'uni',false));
                    else
                        valid_inds = valid_inds & D_ev.task_variable.(cond_var_names{icond_vars}) == T.(cond_var_names{icond_vars}){i_cond};
                    end
                    
                end
                
                data_out{i_session}.response(:,time_mask,valid_inds) = D_ev.response(:,:,valid_inds) - D_cond_mean;
                
            end
            
        end
        
        data_out{i_session}.time(time_mask) = D_ev.time;
        data_out{i_session}.time_rel(time_mask) = D_ev.time_rel;
        data_out{i_session}.time_iev(time_mask) = D_ev.time_iev;
        
    end
end

end

    
    
    
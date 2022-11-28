function [D_out,T] = sort_trials_by_condition(D,sort_vars)

assert(isfield(D,'response'),'no response field in data');
assert(isfield(D,'task_variable'),'no task_variable field in data');

if(~isempty(sort_vars))
    num_cond_labels = length(sort_vars);
    all_task_field_names = fieldnames(D.task_variable);
    assert(sum(ismember(all_task_field_names, sort_vars)) ~= 0,'field does not exist')
    [~, ~, num_trials] = size(D.response);

    % Creating a table for the different task configurations.
    unique_cond_var_values = cell(1,num_cond_labels);

    for icond_label = 1:num_cond_labels
        unique_cond_var_values{icond_label} = unique(D.task_variable.(sort_vars{icond_label}))';
        if(~iscell(unique_cond_var_values{icond_label}))
            unique_cond_var_values{icond_label} = num2cell(unique_cond_var_values{icond_label},1);
        end
    end

    cond_var_values = allcomb(unique_cond_var_values{:});
    T = table();
    for icond_label = 1:num_cond_labels
       T.(sort_vars{icond_label}) = cond_var_values(:,icond_label);
    end

    num_conds = size(cond_var_values,1);
    D_out = cell(num_conds,1);
    
    for icond = 1:num_conds
        D_out{icond} = D;
        T.Properties.RowNames{icond} = sprintf('%s%d','cond-',icond);
        cond_var_names = T.Properties.VariableNames;
        valid_inds = true(num_trials,1);

        for icond_vars = 1:length(cond_var_names)
            if(iscellstr(D.task_variable.(cond_var_names{icond_vars})))
%                 assert(iscellstr(T.(cond_var_names{icond_vars})(icond)),'sort variable should be a cell array of strings')
                valid_inds = valid_inds & cell2mat(cellfun(@(x) strcmp(x, T.(cond_var_names{icond_vars}){icond}), ...
                    D.task_variable.(cond_var_names{icond_vars}),'uni',false));
            else              
                valid_inds = valid_inds & D.task_variable.(cond_var_names{icond_vars}) == T.(cond_var_names{icond_vars}){icond};
            end

        end

        D_out{icond}.response = D.response(:,:,valid_inds);
        for i_tf_name = 1: length(all_task_field_names)
           D_out{icond}.task_variable.(all_task_field_names{i_tf_name}) =  D.task_variable.(all_task_field_names{i_tf_name})(valid_inds); 
        end

    end
           
else
    D_out{1} = D;
    T = table();
end

invalid_conds = cell2mat(cellfun(@(x) isempty(x.response),D_out,'uni',false));
D_out = D_out(~invalid_conds);
if(~isempty(T))
    T = T(~invalid_conds,:);
end

end
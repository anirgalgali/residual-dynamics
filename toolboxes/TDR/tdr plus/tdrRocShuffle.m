function [roc0,rocp] = tdrRocShuffle(data,pars)
% tdrRocShuffle ROC analysis along TDR dimensions
%
% Inputs
%  - data: trial-by-trial simultaneous population data
%  - pars.time_win.(dim_name): [tmin tmax] window to use along dimension dim_name
%  - pars.nshuffle: number of shufflings for p-value computation (def 0)
%  - pars.slide_flag: ROC at all times (1) or across average response (0)
%  - pars.task_index.(index_name): indeces specifying conditions to use [nindx ncd]
%  - pars.roc_index.(index_name): indeces specifying ROC classes [nindx 2]
%
% Output:
%  - roc0: roc values. {ndim 1}[npt ncd] (sliding window) or [ndim 1 ncd] (average)
%  - rocp: pvalue. Same format as roc0, only computed if nshuffle>0
%
% [roc0,rocp] = tdrRocShuffle(data,pars)


% The vector names
vec_names = fieldnames(pars.time_win);
nvec = length(vec_names);

% Shuffle or not
if isfield(pars,'n_shuffle') && ~isempty(pars.n_shuffle)
    n_shuffle = pars.n_shuffle;
else
    n_shuffle = 0;
end

% Sliding window
if isfield(pars,'slide_flag') && ~isempty(pars.slide_flag)
    slide_flag = pars.slide_flag;
else
    slide_flag = 0;
end

% Check if serial or simultaneous recording
if isfield(data,'unit') && ~isempty(data.unit)

    % Talk to me
    disp('not fit for sequential recordings');
    
    % Output
    roc0 = [];
    rocp = [];
    return
    
else
    
    % Dimensions
    [ndm, npt, ntr] = size(data.response);
    
    
    %--- Split trials into two groups (for ROC) ---
    if isfield(pars,'roc_index') && ~isempty(pars.roc_index)
        
        % The indeces
        roc_index = pars.roc_index;
        
        % All the task indeces
        variable_name = fieldnames(roc_index);
        nvr = length(variable_name);
        
        % Number of conditions
        ncd = length(roc_index.(variable_name{1}));
        
        % Initialize trials to average
        jj_roc = false(ncd,ntr);
        
        % Unique variable names
        variable_unique = cell(1,nvr);
        for ivr = 1:nvr
            % All values
            variable_all = data.task_variable.(variable_name{ivr});
            
            % Unique values (without NaNs)
            variable_unique{ivr} = unique(variable_all(~isnan(variable_all)));
        end
        
        % Loop over conditions
        for icd = 1:ncd
            
            % Trials to average
            jj_var = zeros(nvr,ntr);
            
            % Loop over variables
            for ivr = 1:nvr
                
                if roc_index.(variable_name{ivr})(icd)==0
                    % Take all trials
                    jj_var(ivr,:) = true(ntr,1);
                else
                    % Find matching trials
                    jj_var(ivr,:) = ...
                        data.task_variable.(variable_name{ivr}) == ...
                        variable_unique{ivr}(roc_index.(variable_name{ivr})(icd));
                end
            end
            
            % Fullfill constraints on all indeces
            jj_roc(icd,:) = prod(jj_var,1);
                        
        end
                
    else
        
        % Initialize
        jj_roc = NaN;
        
    end
    
    
    %--- Split trials into conditions ---
    if isfield(pars,'task_index') && ~isempty(pars.task_index)
        
        % The indeces
        task_index = pars.task_index;
        
        % All the task indeces
        variable_name = fieldnames(task_index);
        nvr = length(variable_name);
        
        % Number of conditions
        ncd = length(task_index.(variable_name{1}));
        
        % Initialize trials to average
        jj_cnd = false(ncd,ntr);
        
        % Unique variable names
        variable_unique = cell(1,nvr);
        for ivr = 1:nvr
            % All values
            variable_all = data.task_variable.(variable_name{ivr});
            
            % Unique values (without NaNs)
            variable_unique{ivr} = unique(variable_all(~isnan(variable_all)));
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
                end
            end
            
            % Fullfill constraints on all indeces
            jj_cnd(icd,:) = prod(jj_var,1);
                        
        end
                
    else
        
        % Initialize
        ncd = 1;
        jj_cnd = ones(ncd,ntr);
        
    end
    
    
    %--- ROC analysis ---
    % Initialize output
    if slide_flag
        roc0 = cell(nvec,1);
        rocp = cell(nvec,1);
        for ivec = 1:nvec
            roc0{ivec} = zeros(npt,ncd);
            rocp{ivec} = NaN(npt,ncd);
        end
    else
        roc0 = zeros(nvec,1,ncd);
        rocp = NaN(nvec,1,ncd);
    end
    
    % Loop over dimensions
    for ivec = 1:nvec
        
        % Vector name
        vec_name = vec_names{ivec};
        
        % The corresponding regressor
        if isnan(jj_roc)
            isep = strfind(vec_name,'_');
            reg_name = vec_name(1:isep(end)-1);
            % vec_twin = str2double(vec_name(isep(end)+1:end));
        end
        
        % Time window to use
        time_win = pars.time_win.(vec_name);
        vectt = data.time>=time_win(1) & data.time<=time_win(2);
        
        % Loop over conditions
        for icd = 1:ncd
            
            % Sliding window of not
            if slide_flag
                % Sliding window
                vecti = find(vectt);
                npt = length(vecti);
                
                % Loop over time points
                for ipt = 1:npt
                    
                    % Response at this time
                    ravg = squeeze(data.response(ivec,vecti(ipt),:));
                    
                    % Trials to compare
                    if isnan(jj_roc)
                        jtr1 = data.task_variable.(reg_name) == max(data.task_variable.(reg_name));
                        jtr2 = data.task_variable.(reg_name) == min(data.task_variable.(reg_name));
                    else
                        jtr1 = jj_roc(2,:) & jj_cnd(icd,:);
                        jtr2 = jj_roc(1,:) & jj_cnd(icd,:);
                    end
                    
                    % ROC value
                    if n_shuffle > 0
                        [roc0{ivec}(ipt,icd),rocp{ivec}(ipt)] = rocshuf(ravg(jtr1),ravg(jtr2),n_shuffle);
                    else
                        roc0{ivec}(ipt,icd) = rocN(ravg(jtr1),ravg(jtr2),[],0);
                    end
                    
                end
                
            else
                % Average responses within time window
                ravg = squeeze(mean(data.response(ivec,vectt,:),2));
                
                % Trials to compare
                if isnan(jj_roc)
                    jtr1 = data.task_variable.(reg_name) == max(data.task_variable.(reg_name));
                    jtr2 = data.task_variable.(reg_name) == min(data.task_variable.(reg_name));
                else
                    jtr1 = jj_roc(2,:) & jj_cnd(icd,:);
                    jtr2 = jj_roc(1,:) & jj_cnd(icd,:);
                end
                
                % ROC value
                if n_shuffle > 0
                    [roc0(ivec,1,icd),rocp(ivec,1,icd)] = rocshuf(ravg(jtr1),ravg(jtr2),n_shuffle);
                else
                    roc0(ivec,1,icd) = rocN(ravg(jtr1),ravg(jtr2),[],0);
                end
                
            end
        end
    end
end


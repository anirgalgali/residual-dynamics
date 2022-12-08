function [roc0,rocp] = tdrRocShuffle(data,pars)


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



% Initialize
if slide_flag
    roc0 = cell(nvec,1);
    rocp = cell(nvec,1);
else
    roc0 = zeros(nvec,1);
    rocp = NaN(nvec,1);
end

% Loop over vectors
for ivec = 1:nvec
    
    % Vector name
    vec_name = vec_names{ivec};
    
    % The corresponding regressor
    isep = strfind(vec_name,'_');
    reg_name = vec_name(1:isep(end)-1);
    vec_twin = str2double(vec_name(isep(end)+1:end));
    
    % Time window to use
    time_win = pars.time_win.(vec_name);
    vectt = data.time>=time_win(1) & data.time<=time_win(2);
        
    
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
            jtr1 = data.task_variable.(reg_name) == max(data.task_variable.(reg_name));
            jtr2 = data.task_variable.(reg_name) == min(data.task_variable.(reg_name));
            
            % ROC value
            if n_shuffle > 0
                [roc0{ivec}(ipt),rocp{ivec}(ipt)] = rocshuf(ravg(jtr1),ravg(jtr2),n_shuffle);
            else
                roc0{ivec}(ipt) = rocN(ravg(jtr1),ravg(jtr2),[],0);
            end
            
        end
    
    else
        % Average responses within time window
        ravg = squeeze(mean(data.response(ivec,vectt,:),2));
        
        % Trials to compare
        jtr1 = data.task_variable.(reg_name) == max(data.task_variable.(reg_name));
        jtr2 = data.task_variable.(reg_name) == min(data.task_variable.(reg_name));
        
        % ROC value
        if n_shuffle > 0
            [roc0(ivec),rocp(ivec)] = rocshuf(ravg(jtr1),ravg(jtr2),n_shuffle);
        else
            roc0(ivec) = rocN(ravg(jtr1),ravg(jtr2),[],0);
        end
        
    end
end



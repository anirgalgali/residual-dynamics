function [dataM] = mergeDataAcrossSessions(data)
%{ Used to merge a set of data corresponding to different sessions
% Input
% - data(cell array of structs) - of size n_sessions x 1, contianing 
%   the neural data and associated variables for each session
% Output
% - dataM (struct) - merged data structure, where trials and their 
%   associated task variables are pooled across sessions
% Author : Aniruddh Galgali (Jan 2018). Adapted from the TDR toolbox
% developed by Valerio Mante .
%}
n_sessions = length(data);

if(iscell(data))
    dataM = data{1};
else
    dataM = data(1);
end

if(isfield(dataM,'U'))
    dataM.U_align = cell(n_sessions,1);
end

if(isfield(dataM,'X_mu'))
    dataM.X_mu = cell(n_sessions,1);
end

if(isfield(dataM,'dimension'))
    dataM = rmfield(dataM,'dimension');
end

if(isfield(dataM,'task_event_align'))    
    dataM = rmfield(dataM,'task_event_align');   
end

if(isfield(dataM,'task_event_reference'))
    dataM = rmfield(dataM,'task_event_reference');
end

if(isfield(dataM,'task_variable_index'))
    dataM = rmfield(dataM,'task_variable_index');
end

if(isfield(dataM,'dspInfo'))
    dataM = rmfield(dataM,'dspInfo');
end

if(isfield(dataM,'U'))
    dataM = rmfield(dataM,'U');
end


% All task variables
varnames = fieldnames(dataM.task_variable);
nvar = length(varnames);
session_idx = cell(n_sessions,1);
% Loop over task variables
for ivar = 1:nvar
    
    dataM.task_variable.(varnames{ivar}) = [];
    
    if(isfield(dataM,'task_index'))
        dataM.task_index.(varnames{ivar}) = [];
    end
    
    for i_session = 1: n_sessions
        
%         if(~isfield(data(iSessions).task_variable,varnames{ivar}))
%            data(iSessions).task_variable.(varnames{ivar}) = NaN(size(data(iSessions).response,3),1);
%            data(iSessions).task_index.(varnames{ivar}) = NaN(size(data(iSessions).response,3),1);
%         end
        
        if(iscell(data))
            dataM.task_variable.(varnames{ivar}) = ...
                [dataM.task_variable.(varnames{ivar}); data{i_session}.task_variable.(varnames{ivar})];
        else
            dataM.task_variable.(varnames{ivar}) = ...
                [dataM.task_variable.(varnames{ivar}); data(i_session).task_variable.(varnames{ivar})];
        end
        
        if(isfield(dataM,'task_index'))
            if(iscell(data))
                dataM.task_index.(varnames{ivar}) = ...
                    [dataM.task_index.(varnames{ivar}); data{i_session}.task_index.(varnames{ivar})];
            else
                dataM.task_index.(varnames{ivar}) = ...
                    [dataM.task_index.(varnames{ivar}); data(i_session).task_index.(varnames{ivar})];
            end
        end
        if(iscell(data))
            session_idx{i_session} = i_session.*ones(size(data{i_session}.response,3),1);
        else
            session_idx{i_session} = i_session.*ones(size(data(i_session).response,3),1);
        end
        
    end
end
dataM.task_variable.session_id = cell2mat(session_idx);
dataM.response = [];

for i_session = 1: n_sessions
    if(iscell(data))
        dataM.response = cat(3, dataM.response, data{i_session}.response);
        
        if(isfield(dataM,'U_align'))
            dataM.U_align{i_session} = data{i_session}.U;
        end
        if(isfield(dataM,'X_mu'))
            dataM.X_mu{i_session} = data{i_session}.X_mu;
        end
        
    else
        dataM.response = cat(3, dataM.response, data(i_session).response);
        if(isfield(dataM,'U_align'))
            dataM.U_align{i_session} = data(i_session).U;
        end
        if(isfield(dataM,'X_mu'))
            dataM.X_mu{i_session} = data(i_session).X_mu;
        end
    end
end

end
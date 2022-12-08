function dataW = tdrGetTimeWindow(data,timeWin)
% tdrGetTimeWindow Extract data in a time window
%
% Inputs:
%  data: data structure
%  timeWin: [t1 t2]
% 
% Output:
%  dataW: data structure
% 
% dataW = tdrGetTimeWindow(data,timeWin)

% Initialize
dataW = data;

% Find time window
jtt = data.time>=timeWin(1) & data.time<=timeWin(2);

% Extract times
dataW.time = data.time(jtt);

if(isfield(dataW,'time_rel'))
   dataW.time = dataW.time_rel(jtt); 
   dataW = rmfield(dataW,'time_rel');
end

if(isfield(dataW,'time_iev'))
    dataW = rmfield(dataW,'time_iev');
end

% Check if data is sequentially or simultaneously recorded
if isfield(data,'unit') && ~isempty(data.unit)
    
    % Number of units
    nun = length(data.unit);
    
    % Loop over units
    for iun = 1:nun
        % Extract response
        dataW.unit(iun).response = data.unit(iun).response(:,jtt);
    end
    
    % Update alignment info
    if isfield(data,'task_event_align') && ~isempty(data.task_event_align)
        dataW.task_event_align.time_window = timeWin;
    end
    
else
    % Extract response
    dataW.response = data.response(:,jtt,:);
    
    % Update alignment info
    if isfield(data,'task_event_align') && ~isempty(data.task_event_align)
        dataW.task_event_align.time_window = timeWin;
    end
end


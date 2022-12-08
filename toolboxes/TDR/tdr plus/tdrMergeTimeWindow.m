function dataM = tdrMergeTimeWindow(data,pars)



% Merge events
nev = length(data);

% Times
time_rel = data(1).time;
time_iev = ones(size(time_rel));
for iev = 2:nev
    time_rel = cat(2,time_rel,data(iev).time);
    time_iev = cat(2,time_iev,iev*ones(size(data(iev).time)));
end

% New absolute time
time_abs = (0:length(time_rel)-1)/pars.Fs;

% Check if data is sequentially or simultaneously recorded
if isfield(data(1),'unit') && ~isempty(data(1).unit)
    
    % Number of units
    nun = length(data(1).unit);
    
    % Merged data
    dataM = data(1);

    % Loop over units
    for iun = 1:nun
        % Response
        response = data(1).unit(iun).response;
        for iev = 2:nev
            response = cat(2,response,data(iev).unit(iun).response);
        end
        
        % Keep response
        dataM.unit(iun).response = response;
    end
    
    % Keep times
    dataM.time = time_abs;
    dataM.time_rel = time_rel;
    dataM.time_iev = time_iev;
    
else

    % Response and relative time
    response = data(1).response;
    for iev = 2:nev
        response = cat(2,response,data(iev).response);
    end

    % Merged data
    dataM = data(1);
    dataM.response = response;
    dataM.time = time_abs;
    dataM.time_rel = time_rel;
    dataM.time_iev = time_iev;

end

% Remove relative event timing
if isfield(data(1),'task_event_align')
    dataM = rmfield(dataM,'task_event');
end

% Update event alignment
if isfield(data(1),'task_event_align') && ~isempty(data(1).task_event_align)
    for iev = 1:nev
        task_event_align.name{iev} = data(iev).task_event_align.name;
        task_event_align.time_window{iev} = data(iev).task_event_align.time_window;
    end
    dataM.task_event_align = task_event_align;
end

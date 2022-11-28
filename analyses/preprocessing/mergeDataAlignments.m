function [dataM] = mergeDataAlignments(data,varargin)
%{ Used to merge a dataset across different task epochs,
% Input
% - data(cell array of structs) - of size n_alignments x 1, each element is
% data corresponding to a specific temporal alignment of the neural data
% Output
% - dataM (struct) - marged data structure, where the time dimensions is
% combined across different alignment
% Author : Aniruddh Galgali (Jan 2018). Adapted from the TDR toolbox
% developed by Valerio Mante .
%}


nev = length(data);
if(nev <= 1)
    error('Only one alignment provided. Need to provide atleast two')
end

% Times
if(~isfield(data,'time_rel'))

    time_rel = data(1).time;
    time_iev = ones(size(time_rel));
    dt = data(1).time(2) - data(1).time(1) ;
    
    for iev = 2:nev
        time_rel = cat(2,time_rel,data(iev).time);
        time_iev = cat(2,time_iev,iev*ones(size(data(iev).time)));
    end
    
else
    dt = data(1).time(2) - data(1).time(1); 
    time_rel = data.time_rel;
    time_iev = data.time_iev;
    
end


% New absolute time
time_abs = (0:length(time_rel)-1).*dt;
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


if(nargin > 1)
    eventLabels = varargin{1};
    for jj = 1:length(unique(time_iev))
        
        dataM.time_iev_label{jj} = eventLabels{jj};
        
    end
    
end


if isfield(data(1),'task_event_align')
    if(isfield(data(1),'task_event'))
        dataM = rmfield(dataM,'task_event');
    end
end

% Update event alignment
if isfield(data(1),'task_event_align') && ~isempty(data(1).task_event_align)
    for iev = 1:nev
        task_event_align.name{iev} = data(iev).task_event_align.name;
        task_event_align.time_window{iev} = data(iev).task_event_align.time_window;
    end
    dataM.task_event_align = task_event_align;
end

end
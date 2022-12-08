function [data,metadata] = tdrLoadPfcData(datadir,animal,event,unitList,Fd)
% tdrLoadPfcData load PFC units and downsample the responses
%
% Inputs:
%  datadir: the data directory
%  animal: the animal to load ('ar' or 'fe')
%  event: temporal alignment ('Vstim','Tsac',Targ') 
%  unitList: cell array of unit names
%  Fd: sampling frequency after downsampling (def Fd=20)
%
% Outputs:
%  data.unit(i).response: binned firing rates for unit i [ntr npt]
%  data.unit(i).task_variable: all the the task variables [ntr 1]
%     stim_dir: motion coherence of the random dots (>1 = towards pref)
%     stim_col: color coherence of the random dots (>1 = towards green)
%     targ_dir: location of the chosen target (+1=pref, -1=antipref)
%     targ_col: color of the chosen target (+1=green, -1=red)
%     context: motion (+1) or color context (-1)
%     correct: correct (1) or wrong (0)
%     congruent: congruent (1) or incongruent (0)
%     stim_dir2col: motion coherence towards green (>1) or red (<1) target
%     stim_col2dir: color coherence towards pref (>1) or antipref (<1) target
%     stim_trial: trial number
%  data.unit(i).dimension: i-th state space dimension
%  data.time: time axis (s) [1 npt]
%  metadata: All fields not added to data are kept here.
%
% [data,pars] = tdrLoadPfcData(datadir,animal,event,unitList,Fd)


% Inputs
if nargin < 2
    animal = [];
end
if nargin < 3
    event = [];
end
if nargin < 4
    unitList = [];
end
if nargin < 5 || isempty(Fd)
    Fd = 20;
end

% All candidate units
unitDir = dir(datadir);
unitName = setdiff({unitDir(:).name},{'.','..','.DS_Store'})';

% Find units to load
goodUnit = ones(length(unitName),3);
for iun = 1:length(unitName)
    
    % Compare animal
    if ~isempty(animal) && length(unitName) > 2
        goodUnit(iun,1) = strncmp(animal,unitName{iun},2);
    end
    
    % Compare event
    if ~isempty(event) && length(unitName) > length(event)
        goodUnit(iun,2) = any(strfind(unitName{iun},event));
    end
    
    % Compare unit name
    iunder = strfind(unitName{iun},'_');
    if ~isempty(unitList) && length(iunder) > 2
        goodUnit(iun,3) = any(strcmp(unitName{iun}(1:iunder(3)-1),unitList));
    end
    
end

% The units to load
loadUnit = unitName(logical(prod(goodUnit,2)));
nun = length(loadUnit);

% The data fields
datafield = {'response';'task_variable';'time'};

% Loop over units
data = [];
metadata = [];
for iun = 1:nun
    
    % Load unit
    load(fullfile(datadir,loadUnit{iun}));
    
    % The original sampling frequency
    dt = unit.time(2)-unit.time(1);
    Fs = round(1/dt);
    
    % Box size
    R = round(Fs/Fd);
    
    % Resample
    [response,time] = boxresample(Fs*unit.response,R,R,Fs,unit.time(1));
    
    % Keep data
    data.unit(iun,1).response = response;
    data.unit(iun,1).task_variable = unit.task_variable;
    data.unit(iun,1).dimension = sprintf('unit_%03i',iun);
    
    % Metadata fields
    metafield = setdiff(fieldnames(unit),datafield);
    for ifd = 1:length(metafield)
        metadata.unit(iun,1).(metafield{ifd}) = unit.(metafield{ifd});
    end
    
end
data.time = time;


return;








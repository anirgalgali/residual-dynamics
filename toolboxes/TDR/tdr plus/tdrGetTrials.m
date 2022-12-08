function data_get = tdrGetTrials(data,getpars)
% tdrGetTrials extract trials from data
%
% Inputs:
%  data
%  getpars.trial_get (logical or indeces)
%
% Ouputs:
%  data_get: same format as data
%
% data_get = tdrGetTrials(data,getpars)

% Initialize
data_get = [];

% Trials to get
if islogical(getpars.trial_get)
    itrial = find(getpars.trial_get);
else
    itrial = getpars.trial_get;
end

% All fields
fnames = fieldnames(data);
fnames_get = {};

% Get responses
if isfield(data,'response') && ~isempty(data.response)
    data_get.response = data.response(:,:,itrial);
    fnames_get{end+1,1} = 'response';
end

% Get task variables
if isfield(data,'task_variable') && ~isempty(data.task_variable)
    % All task variables
    varnames = fieldnames(data.task_variable);
    nvar = length(varnames);
    
    % Loop over task variables
    for ivar = 1:nvar
        data_get.task_variable.(varnames{ivar}) = ...
            data.task_variable.(varnames{ivar})(itrial,1);
    end
    fnames_get{end+1,1} = 'task_variable';
end

% Get task indeces
if isfield(data,'task_index') && ~isempty(data.task_index)
    % All task indeces
    varnames = fieldnames(data.task_variable);
    nvar = length(varnames);
    
    % Loop over task variables
    for ivar = 1:nvar
        data_get.task_index.(varnames{ivar}) = ...
            data.task_index.(varnames{ivar})(itrial,1);
    end
    fnames_get{end+1,1} = 'task_index';
end

% Add variable-index pairs
if isfield(data_get,'task_variable') && isfield(data_get,'task_index')
    % All task variables
    varnames = fieldnames(data.task_variable);
    nvar = length(varnames);
    
    % Loop over task variables
    for ivar = 1:nvar
        % The unique variable values
        [uval,iuval] = unique(data_get.task_variable.(varnames{ivar}));
        % The pairs
        data_get.task_variable_index.(varnames{ivar}) = [...
            data_get.task_variable.(varnames{ivar})(iuval) ...
            data_get.task_index.(varnames{ivar})(iuval) ];
    end
    fnames_get{end+1,1} = 'task_variable_index';
end

% Get task events
if isfield(data,'task_event') && ~isempty(data.task_event)
    % All task events
    varnames = fieldnames(data.task_event);
    nvar = length(varnames);
    
    % Loop over task variables
    for ivar = 1:nvar
        data_get.task_event.(varnames{ivar}) = ...
            data.task_event.(varnames{ivar})(itrial,1);
    end
    fnames_get{end+1,1} = 'task_event';
end

% Get event reference
if isfield(data,'task_event_reference') && ~isempty(data.task_event_reference)
    data_get.task_event_reference.name = data.task_event_reference.name;
    data_get.task_event_reference.time_absolute = data.task_event_reference.time_absolute(itrial);
    fnames_get{end+1,1} = 'task_event_reference';
end

% Time averaged responses (Added by Aniruddh)

if isfield(data,'time_averaged') && ~isempty(data.time_averaged)
    data_get.time_averaged.response = data.time_averaged.response(:,itrial);
    data_get.time_averaged.population_mean = data.time_averaged.population_mean(itrial);
    data_get.time_averaged.analysis_window = data.time_averaged.analysis_window;
    fnames_get{end+1,1} = 'time_averaged';
end


% Alignment
if isfield(data,'task_event_align') && ~isempty(data.task_event_align)
    data_get.task_event_align = data.task_event_align;
    fnames_get{end+1,1} = 'task_event_align';
end

% Time
if isfield(data,'time') && ~isempty(data.time)
    data_get.time = data.time;
    fnames_get{end+1,1} = 'time';
end

% Relative time
if isfield(data,'time_rel') && ~isempty(data.time_rel)
    data_get.time_rel = data.time_rel;
    fnames_get{end+1,1} = 'time_rel';
end

% Events
if isfield(data,'time_iev') && ~isempty(data.time_iev)
    data_get.time_iev = data.time_iev;
    fnames_get{end+1,1} = 'time_iev';
end

% Dimensions
if isfield(data,'dimension') && ~isempty(data.dimension)
    data_get.dimension = data.dimension;
    fnames_get{end+1,1} = 'dimension';
end

if isfield(data,'dspInfo') && ~isempty(data.dimension)
    data_get.dspInfo = data.dspInfo;
    fnames_get{end+1,1} = 'dimension';
end

if isfield(data,'TaskCondsInfo') && ~isempty(data.TaskCondsInfo)
    data_get.TaskCondsInfo = data.TaskCondsInfo;
    fnames_get{end+1,1} = 'TaskCondsInfo';
end

if isfield(data,'n_trial') && ~isempty(data.n_trial)
    data_get.n_trial = data.n_trial;
    fnames_get{end+1,1} = 'n_trial';
end

% Fields that were missed
missed_fields = setdiff(fnames,fnames_get);
if ~isempty(missed_fields)
    disp('WARNING: trdGetTrials missed fields');
end

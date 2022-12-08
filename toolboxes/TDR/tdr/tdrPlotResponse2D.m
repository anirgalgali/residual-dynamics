function [handles] = tdrPlotResponse2D(data,plotpars)
% tdrPlotResponse2D population trajectories in 2D subspace
%
% Inputs:
%  data: population response
%  plotpars: plotting parameters
%     .dimension: labels of dimensions to plot {2 1}
%     .average_trials: average trials, if more than one found (1 or 0)
%     .colormap: the colormap to draw from
%     .handle.figure: handle to a figure to use (optional)
%     .handle.axes: handle to axes to use (optional)
%     .handle.plot: handle to plots in the axes to use (optional)
%     .type: label type (e.g. 'o-')
%     .markersize: markersize
%     .linewidth: linewidth
%     .markerface: markerface ('full' or 'empty')
%     .dataaspectratio: axis dataspectration
%     .plotboxaspectratio: axis plotboxaspectratio
%     .task_index: indeces to the trials/conditions to plot. Fields of task
%         index refer to task variables in data. Only one field can contain
%         more than one value (which determines the color of the trajectory
%         based on colormap).
%
%         e.g. task_index.context  = 1;
%              task_index.correct  = 2;
%              task_index.stim_dir = [6 5 4]
%         
%         results in 3 trajectories, each corresponding to one the the
%         values in stim_dir.
%
% Outputs:
%  handles: object handles, each an array.
%     .figure
%     .axes
%     .plot
%
% [handles] = tdrPlotResponse2D(data,plotpars)




% Dimensions
[ndm npt ntr] = size(data.response);

% Default markersize
if isempty(plotpars.markersize)
    plotpars.markersize = 6;
end

% Default linewidth
if isempty(plotpars.linewidth)
    plotpars.linewidth = 1;
end

% Time window to plot
if isfield(plotpars,'jtime_plot') && ~isempty(plotpars.jtime_plot)
    jtime = plotpars.jtime_plot;
    it1 = find(plotpars.jtime_plot,1,'first');
    it2 = find(plotpars.jtime_plot,1,'last');
    itt = it1:it2;
else
    jtime = true(1,npt);
    itt = 1:npt;
end

% Indeces that are constrained
index_plot = fieldnames(plotpars.task_index);
nip = length(index_plot);

% All indeces
index_data = fieldnames(data.task_index);
nid = length(index_data);

% Find number of values to plot for each index
index_data_nval = zeros(1,nid);
for iid = 1:nid
    if isfield(plotpars.task_index,index_data{iid}) && ~isempty(plotpars.task_index.(index_data{iid}))
        index_data_nval(iid) = length(plotpars.task_index.(index_data{iid}));
    end
end

% Find task indeces that are variable, fixed, or over which to average
jvariable = find(index_data_nval>1);
jfixed    = find(index_data_nval==1);
% javerage  = find(index_data_nval==0);
% index_name.variable = index_data(jvariable);
% index_name.fixed    = index_data(jfixed);
% index_name.average  = index_data(javerage);

% Colors to use
nval = index_data_nval(jvariable);
lc = linecolors(nval,plotpars.colormap);

% Find trials based on fixed indeces
nfix = length(jfixed);
jtr_fixed = zeros(nfix,ntr);    
for ifix = 1:nfix
    jtr_fixed(ifix,:) = data.task_index.(index_data{jfixed(ifix)}) == ...
        plotpars.task_index.(index_data{jfixed(ifix)});
end

% Find trials based on variable indeces
nvar = index_data_nval(jvariable);
jtr_variable = zeros(nvar,ntr);
for ivar = 1:nvar
    jtr_variable(ivar,:) = data.task_index.(index_data{jvariable}) == ...
        plotpars.task_index.(index_data{jvariable})(ivar);
end

% Dimensions to plot
idimx = find(strcmp(data.dimension,plotpars.dimension{1}));
idimy = find(strcmp(data.dimension,plotpars.dimension{2}));


% PLOT

% The figure
if ~isfield(plotpars,'handle') || ~isfield(plotpars.handle,'figure') || ...
        isempty(plotpars.handle.figure)
    handles.figure = figure;
else
    figure(plotpars.handle.figure);
    handles.figure = plotpars.handle.figure;
end

% The axes
if ~isfield(plotpars,'handle') || ~isfield(plotpars.handle,'figure') || ...
        isempty(plotpars.handle.axes)
    handles.axes = axes;
else
    axes(plotpars.handle.axes);
    handles.axes = plotpars.handle.axes;
end
hold on;

% The plots
if ~isfield(plotpars,'handle') || ~isfield(plotpars.handle,'figure') || ...
        isempty(plotpars.handle.plot);
    hp = [];
else
    hp = plotpars.handle.plot;
end

% Plot responses
for ivar = 1:nvar
    
    % The trial/condition to plot
    jtr = prod(jtr_fixed,1) & jtr_variable(ivar,:);
    
    % The responses
    Rx = data.response(idimx,:,jtr);
    Ry = data.response(idimy,:,jtr);
    
    % Erase response that are not to be plotted
    Rx(:,~jtime,:) = NaN;
    Ry(:,~jtime,:) = NaN;
    
    % Average
    if plotpars.average_trials && sum(jtr)>1
        rx = nanmean(Rx,3);
        ry = nanmean(Ry,3);
    else
        rx = Rx;
        ry = Ry;
    end
    
    % Plot
    np = size(rx,3);
    for ip = 1:np
        hp(end+1) = plot(rx(:,itt,ip),ry(:,itt,ip),plotpars.type);
        
        % Appearance
        switch plotpars.markerface
            case 'full'
                set(hp(end),'color',lc(ivar,:),'markerfacecolor',lc(ivar,:));
            case 'empty'
                set(hp(end),'color',lc(ivar,:),'markerfacecolor',[1 1 1]);
        end
        set(hp(end),'markersize',plotpars.markersize,'linewidth',plotpars.linewidth);
    end
    
end

% Keep handles
handles.plot = hp;

% Axes properties
set(gca,'plotbox',plotpars.plotboxaspectratio,'dataaspectratio',plotpars.dataaspectratio);
hx=xlabel(plotpars.dimension{1});
hy=ylabel(plotpars.dimension{2});
set([hx hy],'interpreter','none');

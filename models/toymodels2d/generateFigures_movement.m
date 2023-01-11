%{ This script generates the main plots in Fig 1 and Fig 2 that concern the 
% models of movement. Plots the recurrent dynamics flow-field along with
% trajectories, time-course of inputs, and estimates of residual dynamcis.
%
%}
%% Load Data
clearvars -except DIRS
close all 
clc
data_load_path = fullfile(DIRS.analysis,'/simulations/toymodels/');
file_name = 'movement_models_newAnalysisPipeline_22-Nov-2018.mat';
load(fullfile(data_load_path,file_name));

%% Extract relevant variables

phase_data = model_analyses.data.phase;
state_data = model_analyses.data.states;
input_data = model_analyses.data.inputs;
dynamics_analysis = model_analyses.analysis;
model_names = fieldnames(phase_data);
rseed = 543612;

%% Setting top-level plot parameters

flow_line_lengths = [60 12 3]; % for movement models only
scale_ = 3.5;
for ii = 1 : length(model_names)

    plotpars.(model_names{ii}).phase.lineWidth  = 1;
    plotpars.(model_names{ii}).phase.LineLength = flow_line_lengths(ii);
    plotpars.(model_names{ii}).phase.arrowShape = [1 0.6 scale_*0.004 scale_*0.06 scale_*0.01 scale_*0.025 scale_*0.000 scale_*0.005];
    plotpars.(model_names{ii}).phase.plot_locs = false;
    % For movement
    plotpars.(model_names{ii}).phase.start_X = 2;
    plotpars.(model_names{ii}).phase.start_Y = 2;
    plotpars.(model_names{ii}).phase.end_X = length(phase_data.(model_names{ii}).grid_points);
    plotpars.(model_names{ii}).phase.end_Y = length(phase_data.(model_names{ii}).grid_points);
    
    plotpars.(model_names{ii}).phase.step_X = 3;
    plotpars.(model_names{ii}).phase.step_Y = 3;


    
  % For movement    
    plotpars.(model_names{ii}).phase.xLim =  [-1.4 1.4];
    plotpars.(model_names{ii}).phase.yLim =  [-1.4 1.4];
    plotpars.(model_names{ii}).phase.xTick =  [-1 -0.5 0 0.5 1 1.5];
    plotpars.(model_names{ii}).phase.yTick=  [-1 -0.5 0 0.5 1];
    plotpars.(model_names{ii}).phase.xTickLabel =  [-1 -0.5 0 0.5 1];
    plotpars.(model_names{ii}).phase.yTickLabel =  [-1 -0.5 0 0.5 1];



    plotpars.(model_names{ii}).phase.figPosition = [130 69 635 636];
    plotpars.(model_names{ii}).phase.Color = [0.75 0.75 0.75];
    plotpars.(model_names{ii}).phase.dotColor = [0 0 0];
    plotpars.(model_names{ii}).phase.dotSize = 4;
    
    plotpars.(model_names{ii}).phase.arrow_start_idx = 5;
    
    plotpars.(model_names{ii}).trajs.single_trials.lineWidth = 1;
    plotpars.(model_names{ii}).trajs.single_trials.figPosition = [130 69 635 636];

    plotpars.(model_names{ii}).trajs.single_trials.doMarkers = false;
    plotpars.(model_names{ii}).trajs.single_trials.Color = [0.3 0.3 0.3];
    plotpars.(model_names{ii}).trajs.single_trials.markersize = 4;
    plotpars.(model_names{ii}).trajs.single_trials.markerFaceCol = [0 0 0];
    plotpars.(model_names{ii}).trajs.single_trials.num_trials_to_plot = 50;
    plotpars.(model_names{ii}).trajs.single_trials.rseed = rseed;
    
    plotpars.(model_names{ii}).trajs.averages.lineWidth = 2;
    plotpars.(model_names{ii}).trajs.averages.figPosition = [130 69 635 636];

    plotpars.(model_names{ii}).trajs.averages.doMarkers = true;
    plotpars.(model_names{ii}).trajs.averages.Color = [0 0 1;1 0 0];
    plotpars.(model_names{ii}).trajs.averages.markersize = 10;
    plotpars.(model_names{ii}).trajs.averages.markerFaceCol = [1 1 1];
    plotpars.(model_names{ii}).trajs.averages.markerEdgeCol = ...
                    plotpars.(model_names{ii}).trajs.averages.Color;
    plotpars.(model_names{ii}).trajs.averages.step = 65;
    
end


%% Plotting phase portrait extended in time (i.e trajectories starting from grid points)

for ii = 1: length(model_names)
    [X,Y] = meshgrid(phase_data.(model_names{ii}).grid_points,phase_data.(model_names{ii}).grid_points);
    V = reshape(phase_data.(model_names{ii}).trajectories,[2 size(phase_data.(model_names{ii}).trajectories,2) size(X,1) size(X,1)]);
    [plot_out.(model_names{ii})] = plotExtendedFlowField_v2(squeeze(V(1,:,:,:)),squeeze(V(2,:,:,:)),plotpars.(model_names{ii}).phase); 
    
end

%% Plotting phase portrait arrows

for ii = 1: length(model_names)
    V = reshape(phase_data.(model_names{ii}).trajectories,[2 size(phase_data.(model_names{ii}).trajectories,2) size(X,1) size(X,1)]);
    locs = squeeze(V(:,1,:,:));
    flow = squeeze(V(:,2,:,:)) - squeeze(V(:,1,:,:));
    [h,ah] = plotFlowField_v2(locs,flow,plotpars.(model_names{ii}).phase,...
         plot_out.(model_names{ii}).fig_handle, plot_out.(model_names{ii}).axes_handle);
     plot_out.(model_names{ii}).fig_handle = h; 
     plot_out.(model_names{ii}).axes_handle = ah;
end

%% Plotting single trials 

for ii = 1:length(model_names)
   
    X = state_data.(model_names{ii}).response;  
    labels = state_data.(model_names{ii}).choice_labels; 
    [plot_out.(model_names{ii})] = plotStateTrajectories(X,plotpars.(model_names{ii}).trajs.single_trials,labels,...
        plot_out.(model_names{ii}).fig_handle, plot_out.(model_names{ii}).axes_handle);
    
end


%% Plotting Condition Averages
        
for ii = 1:length(model_names)
 
    X = state_data.(model_names{ii}).condition_average;  
    labels = unique(state_data.(model_names{ii}).choice_labels); 
    [plot_out.(model_names{ii})] = plotStateTrajectories(X,plotpars.(model_names{ii}).trajs.averages,labels,...
        plot_out.(model_names{ii}).fig_handle, plot_out.(model_names{ii}).axes_handle);
    
    set(plot_out.(model_names{ii}).axes_handle, 'FontSize', 16)
    set(plot_out.(model_names{ii}).axes_handle,'plotbox',[1 1 1]);
    set(plot_out.(model_names{ii}).axes_handle,'plotbox',[1 1 1]);
    set(get(plot_out.(model_names{ii}).axes_handle,'YLabel'),'String','Dim - 2','FontName','Myriad-Pro','FontSize', 16)
    set(get(plot_out.(model_names{ii}).axes_handle,'XLabel'),'String','Dim - 1','FontName','Myriad-Pro','FontSize', 16)

end

%% Plotting inputs time-series

for ii = 1 : length(model_names)

    plotpars.(model_names{ii}).input.linewidth  = 2;
    plotpars.(model_names{ii}).input.lineStyle = {'-';'--'};
    plotpars.(model_names{ii}).input.nCols = 1;
    plotpars.(model_names{ii}).input.nRows = 1;
    plotpars.(model_names{ii}).input.figPosition = [130 69 635 636];
    plotpars.(model_names{ii}).input.Color = [0 0 1;1 0 0];
    
    
    X = permute(input_data.(model_names{ii}).mean,[3 2 1]);
    
    timeAxes = input_data.(model_names{ii}).time;
    [plot_input_handles{ii}] = plotXvTime(X,timeAxes,plotpars.(model_names{ii}).input);
    set(get(plot_input_handles{ii}.axes_handle,'YLabel'),'String','input (a.u)','FontName','Myriad-Pro','FontSize', 16)
    yLims_All(:,ii) = plot_input_handles{ii}.axes_handle.YLim;
end

for ii = 1: length(model_names)
   
%     set(plot_input_handles{ii}.axes_handle,'YLim',[min(yLims_All(1,:)) max(yLims_All(2,:))]);

    % These yLims have been derived from the integration models plot
    
    set(plot_input_handles{ii}.axes_handle,'YLim',[-47.0268 46.9890]);
    set(plot_input_handles{ii}.axes_handle, 'FontSize', 16)
    
end
%% Plotting singular and eigen values

dt_in = 0.001;
dt_out = 0.015;
for ii = 1 : length(model_names)
    X_ev{ii} = NaN(2,size(dynamics_analysis.(model_names{ii}).dynamics{1}.eigVal,2),2);
    X_sv{ii} = NaN(2,size(dynamics_analysis.(model_names{ii}).dynamics{1}.eigVal,2),2);
    X_ang{ii} = NaN(2,size(dynamics_analysis.(model_names{ii}).dynamics{1}.eigVal,2),2);
    X_rf{ii} = NaN(2,size(dynamics_analysis.(model_names{ii}).dynamics{1}.eigVal,2),2);
    for jj = 1:2
        
        [dd_out,ss_out] = rebin_eigenvalue_dynamics(dynamics_analysis.(model_names{ii}).dynamics{jj}.A,...
            dt_in,dt_out,'shuffle');

        X_ev{ii}(:,:,jj) = abs(dd_out); 
        X_sv{ii}(:,:,jj) = ss_out;
        ang_ev = abs(angle(dd_out));
        X_ang{ii}(1,:,jj) = ang_ev(1,:) ;
        X_ang{ii}(2,:,jj) = -1.*ang_ev(1,:) ;

        rf_ev = abs(angle(dd_out))./(2*pi*dt_out);
        X_rf{ii}(1,:,jj) = rf_ev(1,:) ;
        X_rf{ii}(2,:,jj) = -1.*rf_ev(1,:) ;
    
    end
end
cols_lines{1} = [0 0 1;0 0.5 1];
cols_lines{2} = [1 0 0;1 0.5 0];

clear plotpars
for ii = 1 : length(model_names)

    plotpars.(model_names{ii}).dyn.linewidth  = 2;
    plotpars.(model_names{ii}).dyn.lineStyle = {'-','-'};
    plotpars.(model_names{ii}).dyn.nCols = 1;
    plotpars.(model_names{ii}).dyn.nRows = 1;
    plotpars.(model_names{ii}).dyn.figPosition = [1200 275 400 350];
    plotpars.(model_names{ii}).dyn.Color = cols_lines;
    plotpars.(model_names{ii}).dyn.xLim = [0 1];
    pL = dynamics_analysis.(model_names{ii}).dynamics{1}.fitPars.pastLag;  
    timeAxes = state_data.(model_names{ii}).time(pL+1:end-1);

    [plot_eVal_handles{ii}] = plotXvTime(X_ev{ii},timeAxes,plotpars.(model_names{ii}).dyn);
    set(get(plot_eVal_handles{ii}.axes_handle,'YLabel'),'String','|e.v|','FontName','Helvetica','FontSize', 8)
    yLims_All_ev(:,ii) = plot_eVal_handles{ii}.axes_handle.YLim;
    
   
end

clear plotpars
for ii = 1 : length(model_names)

    plotpars.(model_names{ii}).dyn.linewidth  = 2;
    plotpars.(model_names{ii}).dyn.lineStyle = {'-','-'};
    plotpars.(model_names{ii}).dyn.nCols = 1;
    plotpars.(model_names{ii}).dyn.nRows = 1;
    plotpars.(model_names{ii}).dyn.figPosition = [1200 275 400 350];
    plotpars.(model_names{ii}).dyn.Color = cols_lines;
    plotpars.(model_names{ii}).dyn.xLim = [0 1];
    pL = dynamics_analysis.(model_names{ii}).dynamics{1}.fitPars.pastLag;  
    timeAxes = state_data.(model_names{ii}).time(pL+1:end-1);

    [plot_sVal_handles{ii}] = plotXvTime(X_sv{ii},timeAxes,plotpars.(model_names{ii}).dyn);
    set(get(plot_sVal_handles{ii}.axes_handle,'YLabel'),'String','|s.v|','FontName','Helvetica','FontSize', 8)
    yLims_All_sv(:,ii) = plot_sVal_handles{ii}.axes_handle.YLim;
end

clear plotpars
for ii = 1: length(model_names)
    
    plotpars.(model_names{ii}).dyn.linewidth  = 2;
    plotpars.(model_names{ii}).dyn.lineStyle = {'-','-'};
    plotpars.(model_names{ii}).dyn.nCols = 1;
    plotpars.(model_names{ii}).dyn.nRows = 1;
    plotpars.(model_names{ii}).dyn.figPosition = [1200 275 400 350];
    plotpars.(model_names{ii}).dyn.Color = cols_lines;
    plotpars.(model_names{ii}).dyn.xLim = [0 1];
    pL = dynamics_analysis.(model_names{ii}).dynamics{1}.fitPars.pastLag;  
    timeAxes = state_data.(model_names{ii}).time(pL+1:end-1);
    plotpars.(model_names{ii}).dyn.yLim = [-0.1 0.1];
    
    [plot_ang_handles{ii}] = plotXvTime(X_ang{ii},timeAxes,plotpars.(model_names{ii}).dyn);
    set(get(plot_ang_handles{ii}.axes_handle,'YLabel'),'String','angle','FontName','Helvetica','FontSize', 8)
     
end


clear plotpars
for ii = 1: length(model_names)
    
    plotpars.(model_names{ii}).dyn.linewidth  = 2;
    plotpars.(model_names{ii}).dyn.lineStyle = {'-','-'};
    plotpars.(model_names{ii}).dyn.nCols = 1;
    plotpars.(model_names{ii}).dyn.nRows = 1;
    plotpars.(model_names{ii}).dyn.figPosition = [1200 275 400 350];
    plotpars.(model_names{ii}).dyn.Color = cols_lines;
    plotpars.(model_names{ii}).dyn.xLim = [0 1];
    pL = dynamics_analysis.(model_names{ii}).dynamics{1}.fitPars.pastLag;  
    timeAxes = state_data.(model_names{ii}).time(pL+1:end-1);
    plotpars.(model_names{ii}).dyn.yLim = [-0.5 0.5];
    
    [plot_rf_handles{ii}] = plotXvTime(X_rf{ii},timeAxes,plotpars.(model_names{ii}).dyn);
    set(get(plot_rf_handles{ii}.axes_handle,'YLabel'),'String','rot freq (Hz)','FontName','Helvetica','FontSize', 8)
     
end

for ii = 1: length(model_names)
   
    if(ii == 1)
        
        set(plot_eVal_handles{ii}.axes_handle,'YLim',[0.85 1.1]);
        set(plot_sVal_handles{ii}.axes_handle,'YLim',[0.85 1.1]);
        set(plot_ang_handles{ii}.axes_handle,'YLim',[-0.1 0.1]);
        set(plot_rf_handles{ii}.axes_handle,'YLim',[-0.5 0.5]);
    
    elseif(ii == 2)
        set(plot_eVal_handles{ii}.axes_handle,'YLim',[0.7 1.1]);
        set(plot_sVal_handles{ii}.axes_handle,'YLim',[0.7 1.1]);
        set(plot_ang_handles{ii}.axes_handle,'YLim',[-0.1 0.1]);
        set(plot_rf_handles{ii}.axes_handle,'YLim',[-0.5 0.5]);
    
    elseif(ii == 3)
        set(plot_eVal_handles{ii}.axes_handle,'YLim',[0.4 1.1]);
        set(plot_sVal_handles{ii}.axes_handle,'YLim',[0.4 1.1]);
        set(plot_ang_handles{ii}.axes_handle,'YLim',[-0.1 0.1]);
        set(plot_rf_handles{ii}.axes_handle,'YLim',[-0.5 0.5]);

        
    end
    
   
    set(plot_eVal_handles{ii}.axes_handle, 'FontSize', 8)
    set(plot_sVal_handles{ii}.axes_handle, 'FontSize', 8)
    set(plot_ang_handles{ii}.axes_handle, 'FontSize', 8)
    set(plot_rf_handles{ii}.axes_handle, 'FontSize', 8)
end
%%
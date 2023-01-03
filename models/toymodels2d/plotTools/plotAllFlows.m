function plotAllFlows(analyses, model_type, model_idx, c_idx, t_idx, pars)
%{ This function produces the phase portraits of the recurrent, effective and
% residual dynamics (Extended Data Fig. 1). 

% INPUTS
%
% - analyses (struct) - contains the simulated data of the different
%   simulated models
% - model_type (string) - specifies whether the models in 'analyses' are
%   models of decisions OR models of movement
% - model_idx (int) - specfies the index of the model. 
%                      For model_type == 'integration' ('movement'):
%                       '1' = saddle_point (rotations)
%                       '2' = line attractor (channel)
%                       '3' = point attractor
% - c_idx (int) - choice index to plot                      
% - t_idx(int) - time index within trial for which you want to make pltos
%
% Author: Aniruddh Galgali (Oct 2018)
%}
model_names = fieldnames(analyses.data.states);
model_to_plot = model_names{model_idx};

avg_response_at_t = analyses.data.states.(model_to_plot).condition_average(:,t_idx,c_idx);
avg_response_at_tplus1 = analyses.data.states.(model_to_plot).condition_average(:,t_idx + 1, c_idx);


boxX_start = avg_response_at_t(1) - pars.box.halfWidthX;
boxX_end = avg_response_at_t(1) + pars.box.halfWidthX;
boxY_end = avg_response_at_t(2) + pars.box.halfWidthY;
boxY_start = avg_response_at_t(2) - pars.box.halfWidthY;

boxX = [boxX_start boxX_end boxX_end boxX_start];
boxY = [boxY_end boxY_end boxY_start boxY_start];
box_f_patch = [boxX' boxY'];

grid_points_X = boxX_start:pars.box.step:boxX_end;
grid_points_Y = boxY_start:pars.box.step:boxY_end;


switch model_type
    
    case 'movement'
        
        if(model_idx == 3)
            [trajectories_boxed] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,grid_points_X,grid_points_Y,pars.time_steps,false,[0;0]);
        elseif(model_idx == 1)
            [trajectories_boxed] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,grid_points_X,grid_points_Y,pars.time_steps,true,[0;0]);
        elseif(model_idx == 2)
            forcing_term = [analyses.simPars.(model_to_plot).pars.beta;0];
            [trajectories_boxed] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,grid_points_X,grid_points_Y,pars.time_steps,true,forcing_term);
        end

    case 'integration'
      
        [trajectories_boxed] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,grid_points_X,grid_points_Y,pars.time_steps,false,[0;0]);
        
end


[X_b,Y_b] = meshgrid(grid_points_X,grid_points_Y);
V_rec = reshape(trajectories_boxed,[2 size(trajectories_boxed,2) size(X_b,1) size(Y_b,1)]);

locs = squeeze(V_rec(:,1,:,:));
flow_rec = (1./analyses.simPars.(model_to_plot).model_obj.TIME_STEP).*(squeeze(V_rec(:,2,:,:)) - squeeze(V_rec(:,1,:,:)));

if(strcmp(model_type,'movement') && model_idx == 2)

    flow_inp = zeros(size(locs));
    
else

    flow_inp =  repmat(analyses.data.inputs.(model_to_plot).mean(:,t_idx,c_idx),[1 size(locs,2) size(locs,3)]);
end

flow_eff = flow_rec + flow_inp;
switch model_type
    
    case 'movement'
        
        if(model_idx == 3)
            [change_at_avg] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,avg_response_at_t(1),avg_response_at_t(2),2,false,[0;0]);
            flow_at_avg = (1./analyses.simPars.(model_to_plot).model_obj.TIME_STEP).*(change_at_avg(:,2) - change_at_avg(:,1)); 
            flow_at_avg = flow_at_avg + analyses.data.inputs.(model_to_plot).mean(:,t_idx,c_idx);
        elseif(model_idx == 1)
            [change_at_avg] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,avg_response_at_t(1),avg_response_at_t(2),2,true,[0;0]);
            flow_at_avg = (1./analyses.simPars.(model_to_plot).model_obj.TIME_STEP).*(change_at_avg(:,2) - change_at_avg(:,1)); 
        elseif(model_idx == 2)
            forcing_term = [analyses.simPars.(model_to_plot).pars.beta;0];
            [change_at_avg] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,avg_response_at_t(1),avg_response_at_t(2),2,true,forcing_term);
            flow_at_avg = (1./analyses.simPars.(model_to_plot).model_obj.TIME_STEP).*(change_at_avg(:,2) - change_at_avg(:,1)); 
        end
  
    case 'integration'

        [change_at_avg] = generateGridEvolution(analyses.simPars.(model_to_plot).model_obj,avg_response_at_t(1),avg_response_at_t(2),2,false,[0;0]);
        flow_at_avg = (1./analyses.simPars.(model_to_plot).model_obj.TIME_STEP).*(change_at_avg(:,2) - change_at_avg(:,1)); 
        flow_at_avg = flow_at_avg + analyses.data.inputs.(model_to_plot).mean(:,t_idx,c_idx);

end


flow_res_true = flow_eff - flow_at_avg ;


%%

opt_lag = analyses.analysis.(model_to_plot).dynamics{c_idx}.fitPars.pastLag;
time_axes = analyses.simPars.simulation_time;
opt_lambda =  analyses.analysis.(model_to_plot).dynamics{c_idx}.fitPars.regConst;

[A,time_idxs_dyn,X_fs,time_idxs_fs] = estimateDynamicsTwoStageLS(analyses.data.states.(model_to_plot).residuals{c_idx},...
    opt_lag,opt_lambda,ones(length(time_axes),1));

avg_trajectory = analyses.data.states.(model_to_plot).condition_average(:,time_idxs_fs,:);


[analysis_ff] = computeTimeDependentFlowField_v2({A}, [1 0; 0 1], {[1 0; 0 1]},...
    grid_points_X, grid_points_Y, avg_trajectory(:,:,c_idx), X_fs,...
    ones(length(time_idxs_fs),1), time_axes(time_idxs_fs), time_axes(time_idxs_fs));
%%

flow_res_est = (1/(time_axes(2) - time_axes(1))).* analysis_ff.x_dot_dyn{1,1}(:,:,t_idx - opt_lag);
flow_res_est = reshape(flow_res_est,[2 length(grid_points_X) length(grid_points_Y)]);

flow_mean_est = (1/(time_axes(2) - time_axes(1))).* analysis_ff.x_dot_avg{1,1}(:,:,t_idx - opt_lag);
flow_mean_est = reshape(flow_mean_est,[2 length(grid_points_X) length(grid_points_Y)]);

flow_combined_est = (1/(time_axes(2) - time_axes(1))).* analysis_ff.x_dot_combined{1,1}(:,:,t_idx - opt_lag);
flow_combined_est = reshape(flow_combined_est,[2 length(grid_points_X) length(grid_points_Y)]);

%%

pars_flow.start_X = 3;
pars_flow.start_Y = 3;
pars_flow.step_X = 5;
pars_flow.step_Y = 5;
pars_flow.end_X = length(grid_points_X);
pars_flow.end_Y = length(grid_points_Y);
pars_flow.arrowShape = pars.arrowShape;
pars_flow.lineWidth = 1;
pars_flow.LineLength = 0.001 * pars.scale_fac;
pars_flow.xLim =  [boxX_start boxX_end];
pars_flow.yLim =  [boxY_start boxY_end];
pars_flow.xTick =  [boxX_start boxX_end];
pars_flow.yTick=  [boxY_start boxY_end];
pars_flow.xTickLabel =  '';
pars_flow.yTickLabel =  '';
pars_flow.plot_locs = pars.plot_locs;
pars_flow.dotColor = [0 0 0];
pars_flow.dotSize = 0.7;

plotpars_single_trials.lineWidth = 0.5;
plotpars_single_trials.doMarkers = false;
plotpars_single_trials.Color = [0.6 0.6 0.6];
plotpars_single_trials.markersize = 4;
plotpars_single_trials.markerFaceCol = [0 0 0];
plotpars_single_trials.num_trials_to_plot = 50;
plotpars_single_trials.rseed = pars.rseed;

plotpars_averages.lineWidth = 1;

plotpars_averages.doMarkers = true;
plotpars_averages.Color = [0 0 1;1 0 0];
plotpars_averages.markersize = 5;
plotpars_averages.markerFaceCol = [1 1 1];
plotpars_averages.markerEdgeCol = ...
    plotpars_averages.Color;
plotpars_averages.step = 20;

h = figure;
%%
ah1 = subplot(1,6,1);hold(ah1);set(ah1,'plotbox',[1 1 1]);
X = analyses.data.states.(model_to_plot).response;
labels = analyses.data.states.(model_to_plot).choice_labels;
[plot_box] = plotStateTrajectories(X,plotpars_single_trials,labels,...
     h, ah1);

X = analyses.data.states.(model_to_plot).condition_average;
labels = unique(analyses.data.states.(model_to_plot).choice_labels);

[plot_box] = plotStateTrajectories(X,plotpars_averages,labels,...
    plot_box.fig_handle, plot_box.axes_handle);    

plotFlowField_v2(locs,flow_rec,pars_flow,plot_box.fig_handle,plot_box.axes_handle); 
plot(gca,avg_response_at_t(1),avg_response_at_t(2),'o','markersize',8,'markeredgecolor',pars.col_marker,'markerfacecolor',pars.col_marker); 

%%

ah2 = subplot(1,6,2);hold(ah2);set(ah2,'plotbox',[1 1 1]);
X = analyses.data.states.(model_to_plot).response;
labels = analyses.data.states.(model_to_plot).choice_labels;
[plot_box] = plotStateTrajectories(X,plotpars_single_trials,labels,...
     h, ah2);

X = analyses.data.states.(model_to_plot).condition_average;
labels = unique(analyses.data.states.(model_to_plot).choice_labels);

[plot_box] = plotStateTrajectories(X,plotpars_averages,labels,...
    plot_box.fig_handle, plot_box.axes_handle);    
plotFlowField_v2(locs,flow_inp,pars_flow,plot_box.fig_handle,plot_box.axes_handle); 
plot(gca,avg_response_at_t(1),avg_response_at_t(2),'o','markersize',8,'markeredgecolor',pars.col_marker,'markerfacecolor',pars.col_marker); 


%%

ah3 = subplot(1,6,3);hold(ah3);set(ah3,'plotbox',[1 1 1]);
X = analyses.data.states.(model_to_plot).response;
labels = analyses.data.states.(model_to_plot).choice_labels;
[plot_box] = plotStateTrajectories(X,plotpars_single_trials,labels,...
     h, ah3);

X = analyses.data.states.(model_to_plot).condition_average;
labels = unique(analyses.data.states.(model_to_plot).choice_labels);

[plot_box] = plotStateTrajectories(X,plotpars_averages,labels,...
    plot_box.fig_handle, plot_box.axes_handle);    
plotFlowField_v2(locs,flow_eff,pars_flow,plot_box.fig_handle,plot_box.axes_handle); 
plot(gca,avg_response_at_t(1),avg_response_at_t(2),'o','markersize',8,'markeredgecolor',pars.col_marker,'markerfacecolor',pars.col_marker); 

%%
ah4 = subplot(1,6,4);hold(ah4);set(ah4,'plotbox',[1 1 1]);
X = analyses.data.states.(model_to_plot).response;
labels = analyses.data.states.(model_to_plot).choice_labels;
[plot_box] = plotStateTrajectories(X,plotpars_single_trials,labels,...
     h, ah4);

X = analyses.data.states.(model_to_plot).condition_average;
labels = unique(analyses.data.states.(model_to_plot).choice_labels);

[plot_box] = plotStateTrajectories(X,plotpars_averages,labels,...
    plot_box.fig_handle, plot_box.axes_handle);    
plotFlowField_v2(locs,flow_res_true,pars_flow,plot_box.fig_handle,plot_box.axes_handle); 
plot(gca,avg_response_at_t(1),avg_response_at_t(2),'o','markersize',8,'markeredgecolor',pars.col_marker,'markerfacecolor',pars.col_marker); 

%%
ah5 = subplot(1,6,5);hold(ah5);set(ah5,'plotbox',[1 1 1]);
X = analyses.data.states.(model_to_plot).response;
labels = analyses.data.states.(model_to_plot).choice_labels;
[plot_box] = plotStateTrajectories(X,plotpars_single_trials,labels,...
     h, ah5);

X = analyses.data.states.(model_to_plot).condition_average;
labels = unique(analyses.data.states.(model_to_plot).choice_labels);

[plot_box] = plotStateTrajectories(X,plotpars_averages,labels,...
    plot_box.fig_handle, plot_box.axes_handle);    
plotFlowField_v2(locs,flow_combined_est,pars_flow,plot_box.fig_handle,plot_box.axes_handle); 
plot(gca,avg_response_at_t(1),avg_response_at_t(2),'o','markersize',8,'markeredgecolor',pars.col_marker,'markerfacecolor',pars.col_marker); 

%%
ah6 = subplot(1,6,6);hold(ah6);set(ah6,'plotbox',[1 1 1]);
X = analyses.data.states.(model_to_plot).response;
labels = analyses.data.states.(model_to_plot).choice_labels;
[plot_box] = plotStateTrajectories(X,plotpars_single_trials,labels,...
     h, ah6);

X = analyses.data.states.(model_to_plot).condition_average;
labels = unique(analyses.data.states.(model_to_plot).choice_labels);

[plot_box] = plotStateTrajectories(X,plotpars_averages,labels,...
    plot_box.fig_handle, plot_box.axes_handle);    
plotFlowField_v2(locs,flow_res_est,pars_flow,plot_box.fig_handle,plot_box.axes_handle); 
plot(gca,avg_response_at_t(1),avg_response_at_t(2),'o','markersize',8,'markeredgecolor',pars.col_marker,'markerfacecolor',pars.col_marker); 



end



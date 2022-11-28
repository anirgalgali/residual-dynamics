function [varargout] = plot_evsvpanel(analysis, cond_idx_to_plot, base_plot_pars)
%{ This function plots the ev and sv as a function of time (used mostly for the lds models)
% Input
%
% -analysis (Struct) - containing the results of the residual dynamcis
% pipeline
% cond_idx_to_plot - condition (index) for which to plot residual dynamics
% base_plot_pars (struct) - contains high level plotting parameters.

%}
% Generate figure and high-level plotting params

h = figure; set(gcf,'Position',[1742 419 1215 310]);
ax = cell(3,1);
ah = subplot(1, 3, 1);set(ah,'plotbox',[1 1 1]);hold(ah)
col_lines = flipud(cbrewer(base_plot_pars.col_map_type, base_plot_pars.col_map_cols, 20));

%% PLOT EV magnitude and time constants

clear data_to_plot
plotpars_ev.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis.result.final_model{cond_idx_to_plot}.A,1))),:);
plotpars_ev.do_plot_markers = base_plot_pars.do_plot_markers;
plotpars_ev.markersize = base_plot_pars.markersize;
plotpars_ev.linewidth = base_plot_pars.linewidth;
plotpars_ev.ylim = [0.5 1.1];
plotpars_ev.ytick = [0.5 : 0.2 : 1.1];
plotpars_ev.yticklabel = [0.5 : 0.2 : 1.1];
plotpars_ev.yLabel = '| eigen-value | (a.u)';
plotpars_ev.xLabel = 'time (s)';
plotpars_ev.do_plot_stability = true;

data_to_plot = num2cell(abs(analysis.result.final_model{cond_idx_to_plot}.eigVal),2);
plotpars_ev.line_style = repmat(base_plot_pars.line_style,[numel(data_to_plot),1]);
plotpars_ev.marker_style = repmat(base_plot_pars.marker_style,[length(data_to_plot) 1]);
plotpars_ev.markerfacecol = repmat(base_plot_pars.marker_face_color,[length(data_to_plot) 1]);

time_abs = analysis.result.final_model{cond_idx_to_plot}.time;
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_ev);
[ah] = formatXvsTimeplot(ah,[],plotpars_ev.do_plot_stability);


if(isfield(analysis.result, 'final_model_ols') & base_plot_pars.do_plot_ols)
   col_lines_ols = flipud(cbrewer(base_plot_pars.col_map_type,'Greys',20));
   plotpars_ev.col_lines = col_lines_ols(round(linspace(1, size(col_lines_ols,1) - 10, size(analysis.result.final_model_ols{cond_idx_to_plot}.A,1))),:); 
   data_to_plot = num2cell(abs(analysis.result.final_model_ols{cond_idx_to_plot}.eigVal),2);
   [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_ev);
   [ah] = formatXvsTimeplot(ah,[],plotpars_ev.do_plot_stability); 
end

if(isfield(analysis,'sim_opts'))
    
    assert(isfield(analysis.sim_opts, 'binned'),'need binned version of ground truth')
    [~,eigenvalue_groundtruth] = my_eigenshuffle(analysis.sim_opts.binned.dynamics.A,'abs');
    plotpars_ev.col_lines = zeros(size(analysis.sim_opts.binned.dynamics.A,1), 3);
    plotpars_ev.do_plot_markers = false;
    plotpars_ev.linewidth = 1;
    data_to_plot = num2cell(abs(eigenvalue_groundtruth(:,analysis.result.fit_pars.lag{1} + 1 :end-1)),2);
    plotpars_ev.marker_style = repmat({''},[length(data_to_plot) 1]);
    plotpars_ev.line_style = repmat({'-'},[numel(data_to_plot),1]);
    plotpars_ev.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
    [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_ev);
end

yyaxis right
taus_to_show = [-0.1 -0.15 -0.2 -0.35 -0.5 -1 Inf 0.5];
yticks_tau = exp((analysis.sim_opts.bin_size./1000)./taus_to_show);
set(ah,'YColor',[0 0 0],'ylim',plotpars_ev.ylim,'ytick',yticks_tau)
y_ticks_tau = arrayfun( @(x) sprintf('%.2f',x), taus_to_show,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_tau);
set(ah,'Parent',gcf);
ax{1} = ah;

%% PLOT EV phase and rotation frequencies

ah = subplot(1,3,2);set(ah,'plotbox',[1 1 1]);hold(ah)
clear data_to_plot
plotpars_rf.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis.result.final_model{cond_idx_to_plot}.A,1))),:);
plotpars_rf.do_plot_markers = base_plot_pars.do_plot_markers;
plotpars_rf.markersize = base_plot_pars.markersize;
plotpars_rf.linewidth = base_plot_pars.linewidth;
plotpars_rf.ylim = [-0.5 0.5];
plotpars_rf.ytick = [-0.5 : 0.25 : 0.5];
plotpars_rf.yticklabel = [-0.5 : 0.25 : 0.5];
plotpars_rf.yLabel = 'Angle (eigen-value) (a.u)';
plotpars_rf.do_plot_stability = false;
plotpars_rf.xLabel = 'time (s)';

data_to_plot = num2cell(analysis.result.final_model{cond_idx_to_plot}.eigAngle,2);
plotpars_rf.line_style = repmat(base_plot_pars.line_style,[numel(data_to_plot),1]);
plotpars_rf.marker_style = repmat(base_plot_pars.marker_style,[length(data_to_plot) 1]);
plotpars_rf.markerfacecol = repmat(base_plot_pars.marker_face_color,[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_rf);
[ah] = formatXvsTimeplot(ah,[],plotpars_rf.do_plot_stability); 

if(isfield(analysis.result, 'final_model_ols') & base_plot_pars.do_plot_ols)
   col_lines_ols = flipud(cbrewer(base_plot_pars.col_map_type,'Greys',20));
   plotpars_rf.col_lines = col_lines_ols(round(linspace(1, size(col_lines_ols,1) - 10, size(analysis.result.final_model_ols{cond_idx_to_plot}.A,1))),:); 
   data_to_plot = num2cell(abs(analysis.result.final_model_ols{cond_idx_to_plot}.eigAngle),2);
   [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_rf);
   [ah] = formatXvsTimeplot(ah,[],plotpars_rf.do_plot_stability); 
end

if(isfield(analysis,'sim_opts'))
    
    plotpars_rf.col_lines = zeros(size(analysis.sim_opts.binned.dynamics.A,1), 3);
    plotpars_rf.do_plot_markers = false;
    plotpars_rf.linewidth = 1;
    data_to_plot = num2cell(angle(eigenvalue_groundtruth(:,analysis.result.fit_pars.lag{1} + 1 :end-1)),2);
    plotpars_rf.marker_style = repmat({''},[length(data_to_plot) 1]);
    plotpars_rf.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
    plotpars_rf.line_style = repmat({'-'},[numel(data_to_plot),1]);
    [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_rf);
    
end


yyaxis right
set(ah,'YColor',[0 0 0],'ylim',plotpars_rf.ylim,'ytick',plotpars_rf.ytick)
f_to_show = [-2 -1 -0.5 -0.25 0 0.25 0.5 1 2];
yticks_f = 2*pi*(analysis.sim_opts.bin_size/1000)*f_to_show;
set(ah,'YColor',[0 0 0],'ylim',plotpars_rf.ylim,'ytick',yticks_f)
y_ticks_f = arrayfun( @(x) sprintf('%.2f',x), f_to_show ,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_f);
set(ah,'Parent',gcf);
ax{2} = ah;

%% Plot SV magnitudes
 
ah = subplot(1,3,3);set(ah,'plotbox',[1 1 1]);hold(ah)
clear data_to_plot

plotpars_sv.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis.result.final_model{cond_idx_to_plot}.A,1))),:);
plotpars_sv.do_plot_markers = base_plot_pars.do_plot_markers;
plotpars_sv.markersize = base_plot_pars.markersize;
plotpars_sv.linewidth = base_plot_pars.linewidth;
plotpars_sv.ylim = [0.5 1.1];
plotpars_sv.ytick = [0.5 : 0.2 : 1.1];
plotpars_sv.yticklabel = [0.5 : 0.2 : 1.1];
plotpars_sv.yLabel = '|singular-value | (a.u)';
plotpars_sv.xLabel = 'time (s)';
plotpars_sv.do_plot_stability = true;

data_to_plot = num2cell(abs(analysis.result.final_model{cond_idx_to_plot}.sVal),2);
plotpars_sv.line_style = repmat(base_plot_pars.line_style,[numel(data_to_plot),1]);
plotpars_sv.marker_style = repmat(base_plot_pars.marker_style,[length(data_to_plot) 1]);
plotpars_sv.markerfacecol = repmat(base_plot_pars.marker_face_color,[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_sv);
[ah] = formatXvsTimeplot(ah,[],plotpars_sv.do_plot_stability);

if(isfield(analysis.result, 'final_model_ols') & base_plot_pars.do_plot_ols)
   col_lines_ols = flipud(cbrewer(base_plot_pars.col_map_type,'Greys',20));
   plotpars_sv.col_lines = col_lines_ols(round(linspace(1, size(col_lines_ols,1) - 10, size(analysis.result.final_model_ols{cond_idx_to_plot}.A,1))),:); 
   data_to_plot = num2cell(abs(analysis.result.final_model_ols{cond_idx_to_plot}.sVal),2);
   [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_sv);
   [ah] = formatXvsTimeplot(ah,[],plotpars_sv.do_plot_stability); 
end

if(isfield(analysis,'sim_opts'))
    
    [~,singularvalue_groundtruth,~] = my_singularshuffle(analysis.sim_opts.binned.dynamics.A,'right');
    plotpars_sv.col_lines = zeros(size(analysis.sim_opts.binned.dynamics.A,1), 3);
    plotpars_sv.do_plot_markers = false;
    plotpars_sv.linewidth = 1;
    data_to_plot = num2cell(abs(singularvalue_groundtruth(:,analysis.result.fit_pars.lag{1} + 1 :end-1)),2);
    plotpars_sv.marker_style = repmat({''},[length(data_to_plot) 1]);
    plotpars_sv.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
    plotpars_sv.line_style = repmat({'-'},[numel(data_to_plot),1]);
    [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars_sv);
end
ax{3} = ah;
%% Retunr figure and axes handles
varargout{1} = h;
varargout{2} = ax;
end
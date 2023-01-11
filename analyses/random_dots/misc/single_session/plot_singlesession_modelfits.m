%{ This script computes and plots the final residual dynamics fit for 
%  non-alinged residuals from a single session (2008_08_05T). Crossval of
%  the hyperparams needs to be done in a separate step before running this
%  script. See Extended Data Fig. 7
%}

%% Clear Workspace
clearvars -except DIRS
clc
rseed = rng('default');
data_path = './data/analyses/';

%% Specifying the session you want to fit. 

% Make sure that the results of the CV are stored somewhere before running
% this.
dir_contents = dir(data_path);
dir_contents = {dir_contents.name};
session_name_short = '2008_08_05T';
idx_file = cell2mat(cellfun(@(x) ~isempty(strfind(x,session_name_short)),dir_contents,'uni',false));
file_name = dir_contents{idx_file};
load(fullfile(data_path,file_name));

% Properties
save_fig_format = '-pdf';
save_fig_path = '';
do_save_fig = false;

%% Fitting the residual dynamics with optimal hyperparams

% using either the smallest or largest possible alpha obtained through cross-validation
smooth_type = 'max'; 

% Determining the final properties for the residual dynamics fits
opt_final_fit.lag = num2cell(analysis.cv_smooth_pars.lag,2);
opt_final_fit.dim = num2cell(analysis.cv_smooth_pars.dim,2);
opt_final_fit.hankel_order = 5;
opt_final_fit.hankel_rank = analysis.opt_hankel_rank;
opt_final_fit.do_cross_validate_noise_model = false;
nconds = length(analysis.cv_lagdim_opt);
switch smooth_type
    case 'max'
        alpha_all = cell2mat(cellfun(@(x) x.alpha(1), analysis.cv_smooth_opt,'uni',false));
    case 'min'
        alpha_all = cell2mat(cellfun(@(x) x.alpha(end), analysis.cv_smooth_opt,'uni',false));
end

for icond = 1:nconds
    opt_final_fit.alpha{icond} = alpha_all(:,icond);
end

% Fitting the model
[analysis.result] = fit_residual_dynamics(analysis.data_residuals{1}.response, analysis.data_residuals{1}.time, ...
    analysis.data_residuals{1}.time_rel, analysis.data_residuals{1}.time_iev,...
    analysis.data_residuals{1}.task_index.targ_dir, opt_final_fit);

%% Plotting the results
% Choice-index to plot

cidx = 1;
figure;set(gcf,'Position',[449  576 1215 310]);

ah = subplot(1,3,1);set(ah,'plotbox',[1 1 1]);hold(ah)

time_abs = analysis.result.final_model{cidx}.time;
time_iev = analysis.result.final_model{cidx}.time_iev;
time_rel = analysis.result.final_model{cidx}.time_rel;
[timeOfEvents] = findEventMarkers(time_abs,time_rel,time_iev);
col_lines = flipud(cbrewer('seq','GnBu',20));
plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis.result.final_model{cidx}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = 4;
plotpars.linewidth = 2;
plotpars.ylim = [0 1.2];
plotpars.ytick = [0 : 0.2 : 1.2];
plotpars.yticklabel = [0 : 0.2 : 1.2];
plotpars.yLabel = '| eigen-value | (a.u)';
plotpars.xLabel = 'time (s)';
plotpars.do_plot_stability = true;
data_to_plot = num2cell(abs(analysis.result.final_model{cidx}.eigVal),2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_iev,'time_rel',time_rel);
[ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);
yyaxis right
taus_to_show = [-0.1 -0.15 -0.2 -0.35 -0.5 -1 Inf 0.5];
bin_size = analysis.data_residuals{1}.time(2) - analysis.data_residuals{1}.time(1);
yticks_tau = exp(bin_size./taus_to_show);
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_tau)
y_ticks_tau = arrayfun( @(x) sprintf('%.2f',x), taus_to_show,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_tau);
set(ah,'Parent',gcf);

ah = subplot(1,3,2);set(ah,'plotbox',[1 1 1]);hold(ah)
clear plotpars
clear data_to_plot
plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis.result.final_model{cidx}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = 4;
plotpars.linewidth = 2;
plotpars.ylim = [-0.5 0.5];
plotpars.ytick = [-0.5 : 0.25 : 0.5];
plotpars.yticklabel = [-0.5 : 0.25 : 0.5];
plotpars.yLabel = 'Angle (eigen-value) (a.u)';
plotpars.do_plot_stability = false;
plotpars.xLabel = 'time (s)';
plotpars.do_plot_stability = false;
data_to_plot = num2cell(analysis.result.final_model{cidx}.eigAngle,2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_iev,'time_rel',time_rel);
[ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);

yyaxis right
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',plotpars.ytick)
f_to_show = [-2 -1 -0.5 -0.25 0 0.25 0.5 1 2];
yticks_f = 2*pi*(bin_size)*f_to_show;
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_f)
y_ticks_f = arrayfun( @(x) sprintf('%.2f',x), f_to_show ,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_f);
set(ah,'Parent',gcf);

ah = subplot(1,3,3);set(ah,'plotbox',[1 1 1]);hold(ah)
clear plotpars
clear data_to_plot

plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis.result.final_model{cidx}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = 4;
plotpars.linewidth = 2;
plotpars.ylim = [0 1.2];
plotpars.ytick = [0 : 0.2 : 1.2];
plotpars.yticklabel = [0 : 0.2 : 1.2];
plotpars.yLabel = '|singular-value | (a.u)';
plotpars.xLabel = 'time (s)';
plotpars.do_plot_stability = true;
data_to_plot = num2cell(abs(analysis.result.final_model{cidx}.sVal),2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_iev,'time_rel',time_rel);
[ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);
set(ah,'Parent',gcf);

if(do_save_fig)
    fig_name = sprintf('%s%s%s%d%s%s%s',exp_name,'_ModelFitsSingleSession','Choice',cidx,'_',smooth_type,'smoothness');
    export_fig(fullfile(save_fig_path, fig_name),'-painters','-transparent',save_fig_format);
end
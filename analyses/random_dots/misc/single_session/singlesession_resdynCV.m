%{ This script uses the non-aligned residuals from single experiments and
% runs the residual dynamics pipeline with cross-validaiton of
% hyperparameters. 
%
%}
%% Clear Workspace
clearvars -except DIRS
clc
rseed = rng('default');

%% Choosing animal and dataset to fit
animal = 'Tex';
bin_size = 45;
cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics')
data_path = './data/processedneuraldata/';
file_name = sprintf('%s%s%d%s',animal,'_nonaligned_reppasdotsTask_binsize=',bin_size,'ms.mat');
load(fullfile(data_path,file_name)); 

%% Choosing the single session you want to fit

% Choose among the following sessions for Tex -
% { 2008_08_05T, 2008_08_06T, 2009_01_23T, 2008_08_22T, 2008_08_04T,2008_09_11T}
% The above session list is arranged in descending order of the number of 
% trials avaiable in a single contiguous experiment.

n_configs = length(data_aligned.aligned_single_trials);
session_name = 'Tex-05-Aug-2008';

for icfg = 1:n_configs
    
    valid_idxs = cellfun(@(x) strfind(x,session_name),data_aligned.expt_names{icfg},'uni',false);
    if(length(cell2mat(valid_idxs)) > 0)
        config_idx_to_analyze = icfg;
        exp_idx_to_analyze =  find(cell2mat(cellfun(@(x) ~isempty(x), valid_idxs, 'uni',false)));
    end
    
end
data_to_analyze = data_aligned.residuals{config_idx_to_analyze}{exp_idx_to_analyze};
clear data_aligned % Doing this to save memory.

%% Options for the cross-validation of the pipeline

opts_analysis.hankel_order = 5;
opts_analysis.hankel_cv.hankel_ranks = [1 2 4 6 8 10];
opts_analysis.hankel_cv.n_cv = 5;
opts_analysis.hankel_cv.criterion = 'sem';
opts_analysis.hankel_cv.type= 'align';
opts_analysis.lag_cv.n_cv = 1;
opts_analysis.lag_cv.n_folds = 5;
opts_analysis.lag_cv.grid_type = 'uniform';
opts_analysis.lag_cv.do_display = false;
opts_analysis.doLagDimCV = true;
opts_analysis.doFinalFit = false;
opts_analysis.doSmoothCV = false; 
opts_analysis.final_fit.lag_cv_metric = 'mse';
opts_analysis.final_fit.use_common_lagdim_across_conditions = true;

switch opts_analysis.lag_cv.grid_type
    
    case 'uniform'
        opts_analysis.lag_cv.iv_lag = [1 2 3 4];
        opts_analysis.lag_cv.sub_dim = [1 2 4 6 8 10];
        
    case 'random'
        opts_analysis.lag_cv.numPars = 20;
        opts_analysis.lag_cv.max_lag = 10;
        opts_analysis.lag_cv.max_dim = 6;
             
end

if(~(opts_analysis.doSmoothCV))
   fprintf('CV for second-stage hyperparameters are done in a separate step after inspecting first stage CV results\n'); 
end

%% Run cross-validation for lag (l) and dimensionality (d) and store results to disk
   
time_labels = data_to_analyze.time_iev;
trial_labels = data_to_analyze.task_index.targ_dir;
[result_lagdimCV] = run_pipeline(data_to_analyze.response,data_to_analyze.time,...
    data_to_analyze.time_rel,time_labels, trial_labels, opts_analysis, rseed);


file_name = sprintf('%s%s%d%s',session_name,'_nonaligned_hankelagdimCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
save(fullfile(data_path,file_name),'result','-v7.3');

%% Plotting cross-validation results for optimal hankel rank

clear plotpars
plotpars.col_test = [0 0 1; 1 0 0];
plotpars.linewidth = 1.0;
plotpars.bar_length = 0.5;
plotpars.bar_width = 1.0;
plotpars.markersize = 4;
alignments_to_plot = {'stim','sacc'};

    
figure;
nconds = length(result_lagdimCV.cv_hankel);
n_cv = result_lagdimCV.fit_pars.hankel_cv.n_cv;
hankel_ranks = result_lagdimCV.fit_pars.hankel_cv.hankel_ranks;
n_align = length(result_lagdimCV.cv_hankel{1}.stats.test_mean.align);

for ialign = 1:n_align
    
    ah = subplot(1,n_align,ialign); hold(ah);set(ah,'plotbox',[1 1 1]);
    
    for icond = 1:nconds
        
        err_test_mean = result_lagdimCV.cv_hankel{icond}.stats.test_mean.align{ialign};
        err_test_std = result_lagdimCV.cv_hankel{icond}.stats.test_std.align{ialign};
        err_test_sem = (1/sqrt(n_cv)).* err_test_std;
        
        hh_test = ploterr(hankel_ranks, err_test_mean,[],err_test_sem,'-','abshhy', 1.0);
        set(hh_test(1),'Color',plotpars.col_test(icond,:),'linewidth',plotpars.linewidth);
        set(hh_test(2),'Color',plotpars.col_test(icond,:),'linewidth',plotpars.bar_width);
        set(hh_test,'Parent',ah);
        
        hm_test = plot(ah,hankel_ranks, err_test_mean,'o');
        set(hm_test,'markeredgecolor',plotpars.col_test(icond,:) ,'markersize',plotpars.markersize,'markerfacecolor',[1 1 1]);
        set(hm_test,'Parent',ah);
    end
    xlabel(ah,'hankel ranks')
    if(ialign == 1)
        ylabel(ah, 'reconstruction error (mean +/- s.e.m)')
    end
    title(ah,alignments_to_plot{ialign})
end




%% Plotting cross-validation results for lag (l) and dimensionality (d)

do_save_fig = false;
clear plotpars
col_trains = flipud(cbrewer('seq','OrRd',10));
col_tests = flipud(cbrewer('seq','Greys',10));
plotpars.err_type = 'mse';
plotpars.x_type = 'dim';
plotpars.err_scale = 1.0;
plotpars.err_metric = 'sem';
plotpars.linewidth = 1.0;
plotpars.bar_length = 0.5;
plotpars.bar_width = 1.0;
plotpars.do_plot_train = false;
plotpars.markersize = 4;

switch plotpars.x_type
    case 'dim'
        plotpars.col_test = col_tests(1:length(result_lagdimCV.fit_pars.lag_cv.iv_lag),:);
        plotpars.col_train = col_trains(1:length(result_lagdimCV.fit_pars.lag_cv.iv_lag),:);
    case 'lag'
        plotpars.col_test = col_tests(1:length(result_lagdimCV.fit_pars.lag_cv.sub_dim),:);
        plotpars.col_train = col_trains(1:length(result_lagdimCV.fit_pars.lag_cv.sub_dim),:);
end

h = plotCV_lagDim(result_lagdimCV.cv_lagdim,result_lagdimCV.cv_lagdim_hyp,plotpars);
ax = get(h,'Children');
ax = flipud(ax);
for iax = 1: length(ax)
    title(ax(iax),sprintf('%s%s%d%s%d',session_name,', optLag = ',result_lagdimCV.cv_lagdim_opt{iax}.mse_opt.lags(1),...
        ', optDim =',result_lagdimCV.cv_lagdim_opt{iax}.mse_opt.dims(1)),'FontSize',8);
end
set(gcf,'Position',[566 -318 973 1131]);

if(do_save_fig)
    save_fig_path = ''; % enter a valid path indicating where to save the figures
    fig_name = sprintf('%s%s%s%s%.1f%s',session_name,'_LagDimCV',plotpars.err_type,'_meanand',plotpars.err_scale,plotpars.err_metric);
    export_fig(fullfile(save_fig_path, fig_name),'-painters','-transparent',save_fig_format);
end


%% Run cross-validation for alpha(smoothness) and store results to disk

if(~opts_analysis.doSmoothCV)
    
    opts_smooth_cv.alpha_all = [1e-5 1e-3 1e-1 1e1 5e1 1e2 5e2 1e3 1e4];
    opts_smooth_cv.n_cv = 5;
    opts_smooth_cv.n_folds = 5;
    opts_smooth_cv.combine_align = false;
    n_conds = length(result{1}.cv_lagdim_opt);
    opt_lag = 2.*ones(n_conds,1); % These parameters have been obtained through cross-validation
    opt_dim = 4.*ones(n_conds,1); % These parameters have been obtained through cross-validation
    

    result_smoothness = cell(nExps,1);
        
    time_labels = data_to_analyze.time_iev;
    trial_labels = data_to_analyze.task_index.targ_dir;
    
    
    [U_dyn] = compute_dynamics_subspace(data_to_analyze.response,...
        trial_labels,data_to_analyze.time,time_labels,...
        result_lagdimCV.fit_pars.hankel_order, result_lagdimCV.opt_hankel_rank);
    
    [cv_smoothness, cv_smooth_opt] = hyperparamsearch_alpha(data_to_analyze.response,...
        U_dyn, opt_lag, opt_dim, trial_labels, time_labels, opts_smooth_cv, rseed);
    
    
    result_smoothness.cv_smoothness = cv_smoothness;
    result_smoothness.cv_smooth_opt = cv_smooth_opt;
    result_smoothness.cv_pars = opts_smooth_cv;
    result_smoothness.rseed = rseed;
    result_smoothness.U_dyn = U_dyn;
    result_smoothness.opt_lag = opt_lag;
    result_smoothness.opt_dim = opt_dim;

    
    
    file_name = sprintf('%s%s%d%s',session_name,'_aligned_smoothnessCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
    save(fullfile(data_path,file_name),'result_smoothness','-v7.3');
    
end
%% Plotting cross-validation results for alpha.

do_save_fig = false;
clear plotpars
plotpars.col_train = [1 0 0];
plotpars.col_test = [0 0 0];
plotpars.err_scale = 2.0;
plotpars.err_metric = 'std';
plotpars.linewidth = 1.0;
plotpars.bar_length = 0.5;
plotpars.bar_width = 1.0;
plotpars.do_plot_train = true;
plotpars.markersize = 4;
plotpars.marker_type = 'o';
plotpars.err_type = 'mse';
alignments = {'stim';'sacc'};

nconds = size(result_smoothness.cv_smooth_opt,2);
for icond = 1: nconds
    h = plotCV_smoothness(result_smoothness.cv_smoothness{icond},...
        result_smoothness.cv_pars.alpha_all, plotpars);
    ax = get(h,'Children');
    ax = flipud(ax);
    for iax = 1: length(ax)
        title(ax(iax),sprintf('%s%s%s%s%d%s%d',session_name,', ', alignments{iax},': optLag = ',opt_lag(icond),...
            ', optDim =',opt_dim(icond)),'FontSize',8);
    end
    set(gcf,'Position',[566 -318 973 1131]);
    
    if(do_save_fig)
        save_fig_path = '';
        fig_name = sprintf('%s%s%s%d%s%.1f%s',session_name,'_SmoothCV','Choice',icond,'_repeatedCVmeanand',plotpars.err_scale,plotpars.err_metric);
        export_fig(fullfile(save_fig_path, fig_name),'-painters','-transparent',save_fig_format);
    end
    
end

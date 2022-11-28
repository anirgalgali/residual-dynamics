%{ This script runs the cross-validation of the hyper-parameters 
% (d,l and alpha) of the residual dynamics using aligned residuals from a
% single  animal for a specified bin-size. The 'aligned' residuals are
% extracted in a separate previous step and stored to disk prior to running
% this script. 
% 
% We run the pipeline in two separate steps instead of end-to-end, to 
% inspect the intermediate results. First, we run a hyper-param search for
% l (lag), d (dimensionality) and hankel rank (r). We then set an 'optimal' 
% value of l and d by inspecting the cross-validation error (as a function of l and d).
% Using these optimal values, we then do a hyper-param search for alpha
% (smoothness). The results of the CV of l and d aare stored
% in 'result' and the CV for alpha is stored in result_smoothness.
%
% The pipeline is run separately for residuals from each task
% configuration.
% Author: Aniruddh Galgali (Modified: May 2021)
%
%}

%% Clear wokspace
clearvars -except DIRS
clc
rseed = rng('default');

%% Choose animal and bin-size
animal = 'Vito';
bin_size = 45;

%% Load aligned residuals previously stored to disk

% this loads the cell array 'data_aligned', which is a n_configs x 1 cell
% array of struct. Each struct contains the aligned trial and residuals of
% a single task configuration
cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics')
data_path = './data/processedneuraldata/';
file_name = sprintf('%s%s%d%s',animal,'_aligned_reppasdotsTask_binsize=',bin_size,'ms.mat');
load(fullfile(data_path,file_name)); 

%% Set pipeline parameters for cross-validation

opts_analysis.hankel_order = 5;
opts_analysis.hankel_cv.hankel_ranks = [1:12];
opts_analysis.hankel_cv.n_cv = 20;
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
        opts_analysis.lag_cv.iv_lag = [1 2 3 4 5];
        opts_analysis.lag_cv.sub_dim = [1 2 3 4 5 6 8 10 12];
        
    case 'random'
        opts_analysis.lag_cv.numPars = 20;
        opts_analysis.lag_cv.max_lag = 10;
        opts_analysis.lag_cv.max_dim = 6;
             
end

if(~(opts_analysis.doSmoothCV))
   fprintf('CV for second-stage hyperparameters are done in a separate step after inspecting first stage CV results\n'); 
end

%% Run cross-validation for lag (l) and dimensionality (d) and store results to disk

nExps = length(data_aligned.residuals);
result = cell(nExps,1);
for iexp = 1:nExps
   
    time_labels = data_aligned.residuals{iexp}.time_iev;  
    trial_labels = data_aligned.residuals{iexp}.task_index.targ_dir;    
    [result{iexp}] = run_pipeline(data_aligned.residuals{iexp}.response, data_aligned.residuals{iexp}.time,...
        data_aligned.residuals{iexp}.time_rel, time_labels, trial_labels, opts_analysis, rseed);
    fprintf('Done with Config-%d\n',iexp);
end
file_name = sprintf('%s%s%d%s',animal,'_aligned_hankelagdimCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
save(fullfile(data_path,file_name),'result','-v7.3');

%% Plotting cross-validation results for optimal hankel rank

clear plotpars
plotpars.col_test = [0 0 1; 1 0 0];
plotpars.linewidth = 1.0;
plotpars.bar_length = 0.5;
plotpars.bar_width = 1.0;
plotpars.markersize = 4;
alignments_to_plot = {'stim','sacc'};

for iexp = 1:nExps
    
    figure;
    nconds = length(result{iexp}.cv_hankel);
    n_cv = result{iexp}.fit_pars.hankel_cv.n_cv;
    hankel_ranks = result{iexp}.fit_pars.hankel_cv.hankel_ranks;
    n_align = length(result{iexp}.cv_hankel{1}.stats.test_mean.align);
    
    for ialign = 1:n_align
        
        ah = subplot(1,n_align,ialign); hold(ah);set(ah,'plotbox',[1 1 1]);
        
        for icond = 1:nconds
            
            err_test_mean = result{iexp}.cv_hankel{icond}.stats.test_mean.align{ialign};
            err_test_std = result{iexp}.cv_hankel{icond}.stats.test_std.align{ialign};
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
    suptitle(sprintf('%s%s%s%d','HankelCV: ', animal,'Config-',iexp));
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

for iexp = 1:nExps
    switch plotpars.x_type
        case 'dim'
            plotpars.col_test = col_tests(1:length(result{iexp}.fit_pars.lag_cv.iv_lag),:);
            plotpars.col_train = col_trains(1:length(result{iexp}.fit_pars.lag_cv.iv_lag),:);
        case 'lag'
            plotpars.col_test = col_tests(1:length(result{iexp}.fit_pars.lag_cv.sub_dim),:);
            plotpars.col_train = col_trains(1:length(result{iexp}.fit_pars.lag_cv.sub_dim),:);
    end
    
    h = plotCV_lagDim(result{iexp}.cv_lagdim,result{iexp}.cv_lagdim_hyp,plotpars);
    ax = get(h,'Children');
    ax = flipud(ax);
    for iax = 1: length(ax)
        title(ax(iax),sprintf('%s%s%d%s%d%s%d',animal,':Config-',iexp,', optLag = ',result{iexp}.cv_lagdim_opt{iax}.mse_opt.lags(1),...
            ', optDim =',result{iexp}.cv_lagdim_opt{iax}.mse_opt.dims(1)),'FontSize',8);
    end
    set(gcf,'Position',[566 -318 973 1131]);

    if(do_save_fig)
        save_fig_path = ''; % enter a valid path indicating where to save the figures
        fig_name = sprintf('%s%s%d%s%s%s%.1f%s',animal,'_Config=',iexp,'_LagDimCV',plotpars.err_type,'_meanand',plotpars.err_scale,plotpars.err_metric);
        export_fig(fullfile(save_fig_path, fig_name),'-painters','-transparent',save_fig_format);
    end
end

%% Run cross-validation for alpha(smoothness) and store results to disk

if(~opts_analysis.doSmoothCV)
    
    opts_smooth_cv.alpha_all = [1e-5 1e-3 1e-1 1e1 5e1 1e2 2e2 5e2 1e3 1e4 1e5];
    opts_smooth_cv.n_cv = 5;
    opts_smooth_cv.n_folds = 5;
    opts_smooth_cv.combine_align = false;
    n_conds = length(result{1}.cv_lagdim_opt);

    if(strcmp(animal,'Tex'))
        opt_lag = 3.*ones(n_conds,1); % These parameters have been obtained through coss-validation
        opt_dim = 8.*ones(n_conds,1);
    elseif(strcmp(animal,'Vito'))
        opt_lag = 3.*ones(n_conds,1);
        opt_dim = 8.*ones(n_conds,1);
    end

    result_smoothness = cell(nExps,1);
    for iexp = 1:nExps
        
        time_labels = data_aligned.residuals{iexp}.time_iev;
        trial_labels = data_aligned.residuals{iexp}.task_index.targ_dir;
        
        
        [U_dyn] = compute_dynamics_subspace(data_aligned.residuals{iexp}.response,...
            trial_labels,data_aligned.residuals{iexp}.time,time_labels,...
            result{iexp}.fit_pars.hankel_order, result{iexp}.opt_hankel_rank);
        
        [cv_smoothness, cv_smooth_opt] = hyperparamsearch_alpha(data_aligned.residuals{iexp}.response,...
            U_dyn, opt_lag, opt_dim, trial_labels, time_labels, opts_smooth_cv, rseed);
        
     
        result_smoothness{iexp}.cv_smoothness = cv_smoothness;
        result_smoothness{iexp}.cv_smooth_opt = cv_smooth_opt;
        result_smoothness{iexp}.cv_pars = opts_smooth_cv;
        result_smoothness{iexp}.rseed = rseed;
        result_smoothness{iexp}.U_dyn = U_dyn;
        result_smoothness{iexp}.opt_lag = opt_lag;
        result_smoothness{iexp}.opt_dim = opt_dim;
        fprintf('Done with Config-%d\n',iexp);

    end
    
    file_name = sprintf('%s%s%d%s',animal,'_aligned_smoothnessCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
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

for iexp = 1:nExps
    nconds = size(result_smoothness{iexp}.cv_smooth_opt,2);
    for icond = 1: nconds
        h = plotCV_smoothness(result_smoothness{iexp}.cv_smoothness{icond},...
            result_smoothness{iexp}.cv_pars.alpha_all, plotpars);
        ax = get(h,'Children');
        ax = flipud(ax);
        for iax = 1: length(ax)
            title(ax(iax),sprintf('%s%s%d%s%s%s%d%s%d',animal,':Config-',iexp,', ', alignments{iax},': optLag = ',opt_lag(icond),...
                ', optDim =',opt_dim(icond)),'FontSize',8);
        end
        set(gcf,'Position',[566 -318 973 1131]);
        
        if(do_save_fig)
            save_fig_path = '';
            fig_name = sprintf('%s%s%d%s%s%d%s%.1f%s',animal,'_Config=',iexp,'_SmoothCV','Choice',icond,'_repeatedCVmeanand',plotpars.err_scale,plotpars.err_metric);
            export_fig(fullfile(save_fig_path, fig_name),'-painters','-transparent',save_fig_format);
        end
        
    end
end
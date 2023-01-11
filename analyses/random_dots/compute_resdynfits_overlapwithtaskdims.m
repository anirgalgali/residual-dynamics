%{This script compute the residual dynamics for all task configurations,
% epochs, and conditions of a specified dataset (monkey/bin-size). These
% fits are the final fits that take into account the optimal value of all
% hyperparameters. This script also compute the aligned "task-planes", and
% the overlap between the task plaes and the residual dynamics.
%
% Running this script requires storing results of the cross-valiation of 
% the hyperparams to disk.
%
% Author - Aniruddh Galgali (Nov 2018, Modified : May 2020)
%}
clearvars -except DIRS
clc
rseed = rng('default');
close all

%% Choosing dataset

animal = 'Tex';
bin_size = 45;
do_save_vars = false;

%% Load aligned data

data_path = './data/analyses/';
file_name = sprintf('%s%s%d%s',animal,'_aligned_reppasdotsTask_binsize=',bin_size,'ms.mat');
load(fullfile(data_path,file_name)); 

%% Load sross-validation results (lagDimCV and smoothCV)

file_name = sprintf('%s%s%d%s',animal,'_aligned_hankelagdimCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
R = load(fullfile(data_path,file_name));
result_lagdimCV = R.result;

file_name = sprintf('%s%s%d%s',animal,'_aligned_smoothnessCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
load(fullfile(data_path,file_name));
clear R
%% Estimating condition-averaged trajectories in selected 2D "task" planes


alignments = {'stim';'sacc'};

% jPC params
jPC_pars.n_jpc_planes = 2;
jPC_pars.do_pre_smooth = true;
jPC_pars.time_win_analyze = {[0.49 1], [-0.25 0.25]}; 

%choice-time pars
choice_time_pars.condition_vars = {'targ_dir'};
pref_idx = [2 1 1 1]; % config-1 has flipped indices for the choice (can be made cleaner)
anti_idx = [1 2 2 2];

n_configs = length(data_aligned.aligned_single_trials);
analysis_condavg = cell(n_configs,1);

for iexp = 1:n_configs
    
    choice_time_pars.pref_idx = pref_idx(iexp);
    choice_time_pars.anti_idx = anti_idx(iexp);
    trial_labels = data_aligned.aligned_single_trials{iexp}.task_index.targ_dir;
    [analysis_condavg{iexp}] = extract_alignedtaskplanes...
        (data_aligned.aligned_single_trials{iexp}, trial_labels, jPC_pars, choice_time_pars, alignments);
end

choice_time_pars.pref_idx = pref_idx;
choice_time_pars.anti_idx = anti_idx;

jPC_pars.alignments = alignments;
choice_time_pars.alignments = alignments;

%% Estimating dynamics for all configurations in the dataset


n_dim = 8;
lag = 3;
alpha = [200;50];
opts_final.dim = {n_dim;n_dim};
opts_final.lag = {lag;lag};
opts_final.alpha = {[alpha];[alpha]};
opts_final.order_across_alignment = true;
opts_final.do_cross_validate_noise_model = false;


do_bootstrap_fits = false; % !!! Set this to false by default !!!
% WARNING: Bootstrap take a loooooong time even for a single task
% configuration. Better to run bootstrap separately for each task
% configuration 
if(do_bootstrap_fits)
   
    opts_final.do_boot = true;
    opts_final.boot_opts.n_boot = 1000;
    opts_final.boot_opts.type = 'ci';
    opts_final.boot_opts.use_parallel = true;
    
end


result_final = cell(1,n_configs);
opts_final_ = opts_final;

for iexp = 1:n_configs
    
    time_labels = data_aligned.residuals{iexp}.time_iev;
    trial_labels = data_aligned.residuals{iexp}.task_index.targ_dir;
    opts_final_.hankel_order = result_lagdimCV{iexp}.fit_pars.hankel_order;
    opts_final_.hankel_rank = result_lagdimCV{iexp}.opt_hankel_rank;
    
    [result_final{iexp}] = fit_residual_dynamics(data_aligned.residuals{iexp}.response, ...
        data_aligned.residuals{iexp}.time, data_aligned.residuals{iexp}.time_rel,...
        time_labels, trial_labels, opts_final_);
 
end


%% Computing alignment with condition-averages and aggregating results across all configurations in the dataset

[data_agg] =  compute_ovelap_resdyn_avgtaskplanes(analysis_condavg, result_final);

%% Computing pairwise angles betwen the different aligned task planes
pairwise_angle_task_planes = table();
c_count = 1;

for iexp = 1: length(analysis_condavg)
   for ialign = 1: length(alignments)
        
       dir_labels = analysis_condavg{iexp}.(alignments{ialign}).direction_labels;
       
       for idir1 = 1:length(dir_labels)-1
           for idir2 = idir1+1:length(dir_labels)
          
               pairwise_angle = rad2deg(subspace(analysis_condavg{iexp}.(alignments{ialign}).directions{idir1},...
                   analysis_condavg{iexp}.(alignments{ialign}).directions{idir2}));
                   
               pairwise_angle_task_planes.(['angle_' dir_labels{idir1} '_' dir_labels{idir2}])(c_count) = ...
                   pairwise_angle;              
               
           end       
       end 
       pairwise_angle_task_planes.alignment(c_count) = alignments(ialign);
       pairwise_angle_task_planes.config_idx(c_count) = iexp;
       c_count = c_count + 1;
   end 
end


%% Sotring results
if(do_save_vars)
    save_file_path = './data/analyses/';
    result.data_identifier = sprintf('%s%s%d%s%d%s%d%s%d',animal,'_ndim=',n_dim,'_lag=',lag,'_alpha=',opts_final.alpha{1}(1),'_',opts_final.alpha{1}(2));
    result.resdynfit = data_agg;
    result.resdynfit_final_raw = result_final;
    result.resdynfit_cv = result_lagdim_cv_all;
    result.cond_avgs = analysis_condavg;
    result.resdynfit_pars = opts_final;
    result.cond_avg_pars.jPC = jPC_pars;
    result.cond_avg_pars.choice_time = choice_time_pars;
    result.task_subspace_overlap = pairwise_angle_task_planes;
    file_name = strcat(result.data_identifier,'_allconfigsresdynresults','.mat');
    save(fullfile(save_file_path,file_name),'result','-v7.3')
end
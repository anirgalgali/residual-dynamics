%{ This script generates simulated neural data from a two-area RNN, consisting 
% of two areas or "modules" - PPC and PFC. The within-area recurrence strength
% and the across-area recurrence strength areswept across a range of values. 
% The other parameters of the network are set to values specified in the 
% original paper by Murray et al, 2017).

%}
%% Clear workspace

close all
clc
rseed = rng('default');
rng(rseed);
do_save_vars = false;
%% Simulation parameters

sim_pars.dt = 0.0001;
sim_pars.T = 12000; 
sim_pars.T_on = 4000;
sim_pars.K = 6000;
sim_pars.lat_dim = 2; % number of latent dims (populations) per area moduel (always = 2)
sim_pars.obs_dim = 10; % number of observed neurons per area
sim_pars.doFbk = false;
sim_pars.area_to_decode = 'pfc';
sim_pars.rseed = rseed;

%% Input parameters that determine input to PPC module of the network
input_pars.mu0 = 25;
input_pars.I0 = 0.334; %nA
input_pars.sigma_noise = 0.07; %nA 
input_pars.coh_set = [0];
input_pars.dir_stim_set = [-1 1];
input_pars.x0 = 0.01*ones(2*sim_pars.lat_dim,1);

%% Parameters determining the observation model

obs_pars.bin_size = [45];
C = RandOrthMat(sim_pars.obs_dim);
obs_pars.C_ppc = C(:,1:sim_pars.lat_dim);
obs_pars.C_pfc = C(:,sim_pars.lat_dim + 1: 2*sim_pars.lat_dim);
obs_pars.C_obs = blkdiag(obs_pars.C_ppc, obs_pars.C_pfc);

obs_pars.sigma_obs_noise = 0.0006; % original value
obs_pars.R = diag(obs_pars.sigma_obs_noise.*ones(size(obs_pars.C_obs,1),1));
obs_pars.residual_condition_vars{1} = {'targ_dir','dots_coh'};

%% Residual dynamics pipeline parameters

opts_analysis.hankel_order = 5;
opts_analysis.hankel_cv.hankel_ranks = [2 3 4 5 6 8 10];
opts_analysis.hankel_cv.n_cv = 20;
opts_analysis.hankel_cv.criterion = 'min';
opts_analysis.hankel_cv.type= 'combined';

opts_analysis.lag_cv.n_cv = 1;
opts_analysis.lag_cv.n_folds = 5;
opts_analysis.lag_cv.grid_type = 'uniform';
opts_analysis.lag_cv.do_display = false;

opts_analysis.doSmoothCV = true;
opts_analysis.doFinalFit = false;
opts_analysis.doLagDimCV = true;

switch opts_analysis.lag_cv.grid_type
    
    case 'uniform'
        opts_analysis.lag_cv.iv_lag = [1 2 3 4 5 8];
        opts_analysis.lag_cv.sub_dim = [2 3 4 6 8 10];
        
    case 'random'
        opts_analysis.lag_cv.numPars = 20;
        opts_analysis.lag_cv.max_lag = 10;
        opts_analysis.lag_cv.max_dim = 6;
             
end

opts_analysis.lag_cv.metric = 'mse';
opts_analysis.smooth_cv.alpha_all = [1e-1 1e0 1e1 5e1 1e2 5e2 1e3 5e3 1e4 1e5 1e6 1e8];
opts_analysis.smooth_cv.n_cv = 10;
opts_analysis.smooth_cv.n_folds = 5;
opts_analysis.smooth_cv.combine_align = true;

opts_analysis.final_fit.lag_cv_metric = 'mse';
opts_analysis.final_fit.use_common_lagdim_across_conditions = true;
if(~opts_analysis.doSmoothCV)
    opts_analysis.final_fit.smoothalpha = {1e4};
end
opts_analysis.final_fit.do_cross_validate_noise_model = true;
opts_analysis.final_fit.alpha_cv_metric = 'mse';

%% Create network parameters

self_wt_all = [0.20: 0.04: 0.40]; % weight values of self recurrence
cross_wt_all = [0:0.04:0.24]; % weight values of fbk/fwd connections

[self_wts, cross_wts] = meshgrid(self_wt_all' , cross_wt_all');
total_params = length(self_wt_all)*length(cross_wt_all);
act_fn.a = 270;
act_fn.b = 108;
act_fn.c = 0.154;
act_fn.gamma = 0.641;
net_pars = repmat(struct('tau_s',0.06, 'J_ext',5.2e-4,'tau_ampa', 0.002,...
    'J_t',[0.28387 0;0 0.28387],'act_fn',act_fn),size(self_wts));

%% Simulate data and fit the residual dynamics pipeline

data = repmat({struct()},size(self_wts));
areas = {'pfc','ppc','both'};
analysis = repmat({struct()},[size(self_wts) length(areas)]);
tmax = cell(size(self_wts));
do_shuffle = false(size(net_pars,1),size(net_pars,2));


fit_opts_.hankel_criterion = 'sem';
fit_opts_.hankel_type = 'combined';

fit_opts_.bootstrp.n_boot = 1000;
fit_opts_.bootstrp.use_parallel = true;
fit_opts_.bootstrp.type = 'direct'; % confidence intervals


fit_opts_.lag = 3;
for iw = 1:size(net_pars,1)
    for ic = 1:size(net_pars,2)
        
        
        
        net_pars(iw,ic).self_ppc = self_wts(iw,ic);
        net_pars(iw,ic).self_pfc = self_wts(iw,ic);
        net_pars(iw,ic).fwd_ppc_pfc = cross_wts(iw,ic);
        
        if(sim_pars.doFbk)
            net_pars(iw,ic).fbk_pfc_ppc = cross_wts(iw,ic);
        else
            net_pars(iw,ic).fbk_pfc_ppc = 0;
        end
        
       if(self_wts(iw,ic) == 0.36 & cross_wts(iw,ic) == 0.08 & ~sim_pars.doFbk)
           do_shuffle(iw,ic) = true;  
       end
        
        sim_pars.do_shuffle = do_shuffle(iw,ic);
        [data{iw,ic}] = generate_responses(input_pars, net_pars(iw,ic), sim_pars, obs_pars);
        [~,~,choice_labels] = unique(data{iw,ic}.task_variable.targ_dir);
        
        for ia = 1: length(areas)
            
            switch areas{ia}
                
                case 'ppc'
                    
                    mask_dim = [true(sim_pars.obs_dim,1); false(sim_pars.obs_dim,1)];
                    fit_opts_.max_dim = 2;
                case 'pfc'
                    
                    mask_dim = [false(sim_pars.obs_dim,1); true(sim_pars.obs_dim,1)];
                    fit_opts_.max_dim = 2;
                    
                case 'both'
                    
                    mask_dim = [true(sim_pars.obs_dim,1) ; true(sim_pars.obs_dim,1)];
                    fit_opts_.max_dim = 4;
            end
            
            [analysis{iw,ic,ia}] = run_pipeline(data{iw,ic}.residuals(mask_dim,:,:), ...
                data{iw,ic}.time, data{iw,ic}.time_rel, data{iw,ic}.time_iev, choice_labels, opts_analysis, rseed);
            
            [analysis{iw,ic,ia}.resdyn] = extract_residualdynamics_twoarearnn(data{iw,ic}.residuals(mask_dim,:,:),...
                data{iw,ic}.time, data{iw,ic}.time_rel, data{iw,ic}.time_iev, choice_labels, analysis{iw,ic,ia}, fit_opts_);
            
            
            if(strcmp(areas{ia},'ppc') | strcmp(areas{ia},'pfc'))
                
                choice_mode_obs = obs_pars.C_obs * data{iw,ic}.task_modes.(areas{ia}).choice.direction;
                [analysis{iw,ic,ia}.resdyn_along_choice] = extract_residualdynamics_twoarearnn(data{iw,ic}.residuals(mask_dim,:,:),...
                data{iw,ic}.time, data{iw,ic}.time_rel, data{iw,ic}.time_iev, choice_labels, analysis{iw,ic,ia}, fit_opts_,choice_mode_obs(mask_dim,:)); 
            end
            
            
            
            if(isfield(data{iw,ic},'shuffled'))
                
                [~,~,choice_labels_shuf] = unique(data{iw,ic}.shuffled.task_variable.targ_dir);
                
                [analysis{iw,ic,ia}.shuffled] = run_pipeline(data{iw,ic}.shuffled.residuals(mask_dim,:,:), ...
                    data{iw,ic}.shuffled.time, data{iw,ic}.shuffled.time_rel, data{iw,ic}.shuffled.time_iev,...
                    choice_labels_shuf, opts_analysis, rseed); 
                
                [analysis{iw,ic,ia}.shuffled.res_dyn] = extract_residualdynamics_twoarearnn(data{iw,ic}.shuffled.residuals(mask_dim,:,:),...
                    data{iw,ic}.shuffled.time, data{iw,ic}.shuffled.time_rel, data{iw,ic}.shuffled.time_iev,...
                    choice_labels_shuf, analysis{iw,ic,ia}.shuffled, fit_opts_);
            end
                    
        end
        fprintf('Done with self = %.2f, cross = %.2f\n', self_wts(iw,ic), cross_wts(iw,ic)); 
    end
    
end
sim_pars.do_shuffle = do_shuffle;

%% Saving cross-validated fit results to disk

result.sim_pars = sim_pars;
result.obs_pars = obs_pars;
result.input_pars = input_pars;
result.net_pars = net_pars;
result.data = data;
result.areas = areas;
result.analysis = analysis;

if(do_save_vars)
    
    save_file_path ='./'
    file_name = sprintf('%s%d%s%d%s%.2f%s%.4f%s',' TwoAreaModel_weightSweepsModelEstimates_mu0=',...
        input_pars.mu0,'Hz_zeroCohBinsize=',obs_pars.bin_size,'_sigmaLatent=',input_pars.sigma_noise,...
        '_sigmaNoise=',obs_pars.sigma_obs_noise,'.mat');
    
    save(fullfile(save_file_path, file_name),'analysis','-v7.3');

end
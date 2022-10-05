
%{
---------------------------------------------------------------------------
This script produces the simulations of time-varying dynamical systems 
with either a gaussian or poisson observation process and fits the residual
dynamics pipeline to the generated data. See Supplementary Figure 3.
---------------------------------------------------------------------------
Author: Aniruddh Galgali
%}

%% clear workspace
clearvars -except DIRS
clc
close all

%% Setting seed and model ID

rseed = 24321; %24321;
rng(rseed)
do_save_vars = false;
save_file_path = './data/simulations';

%% Define the basic parameters of LDS simulations

% Set random seed
options.rseed = rseed;

% Specify whether to simulate a continuous time or discrete time dynamical
% system and a poisson or gaussian observation process
options.dynamics_model = 'continuous';
options.observation_model = 'gaussian';
options.time_varying_latent_noise = false;
% Setting the latent/ob dims and number of trial/time steps and the
% time-step
options.T = 1500;

switch options.observation_model
    
    case 'poisson'
        options.obs_dim = 100;
        
    case 'gaussian'
        options.obs_dim = 20;
        
end

options.lat_dim = 3;
options.num_trials = 5000;
options.dt = 0.001; % simulation step size
options.bin_size = 45; % bin-size for analysis


%% Generate dynamics matrices

% There are 3 different models that one can simulate based on the value of
% the variable modelID, which is described below:
% {'1': change of ev, '2': change of basis, '3': change to rotational'}

modelID = 3;

% Define a time-varying dynamics for the selected model type

if(modelID == 1)
    
    T_switch_dyn = 800;
    options.modelType = 'changeLambda';
    options.dyn_type = 'normal';
    
    eig_vals1 = [-1;-3;-5];
    eig_vals2 = [-2;-4;-6];
    
    eig_vals = cat(2, eig_vals1.*ones(options.lat_dim,T_switch_dyn),...
        eig_vals2.*ones(options.lat_dim,options.T - T_switch_dyn));

    V = randn(options.lat_dim,options.lat_dim);
    
    switch options.dyn_type
        
        case 'normal'
    
            V = myqr(V);
    end
    
    A = NaN(options.lat_dim,options.lat_dim,options.T);
    
    for tt = 1:options.T
        
        A(:,:,tt) =   V * diag(eig_vals(:,tt)) * inv(V);
    
    end

    options.A = A;
  
elseif(modelID == 2)
    
    T_switch_dyn = 800;
    options.modelType = 'changeBasis';
    options.dyn_type = 'normal';
    eig_vals = [-1;-3;-5];
    eig_vals = eig_vals.*ones(options.lat_dim,options.T);
    
    V = randn(options.lat_dim,options.lat_dim);
    V1 = myqr(V);
    V2 = normc(V);
    
    A = NaN(options.lat_dim,options.lat_dim,options.T);
    
    for tt = 1:options.T
        
        if(tt <= options.T_switch_dyn)
            A(:,:,tt) =   V1 * diag(eig_vals(:,tt)) * inv(V1);
        else
            A(:,:,tt) =   V2 * diag(eig_vals(:,tt)) * inv(V2);
        end
    
    end
    
    options.A = A;

    
elseif(modelID == 3)
    
    T_switch_dyn = 800;
    options.modelType = 'rotation';
    rot_freq = 1;
    eig_vals1 = [-1;-3;-5];
    eig_vals2 = [-2 + 1i*(2*pi*rot_freq) ;...
        -2 - 1i*2*pi*rot_freq ; -5];
    eig_vals = cat(2,eig_vals1.*ones(options.lat_dim,T_switch_dyn),...
        eig_vals2.*ones(options.lat_dim,options.T - T_switch_dyn));
    
    V = randn(options.lat_dim,options.lat_dim);
    V1 = myqr(V);
    V2 = V1;
    V2(:,1) = (1/sqrt(2)).*(V1(:,1) + 1i*V1(:,2));
    V2(:,2) = (1/sqrt(2)).*(V1(:,1) - 1i*V1(:,2));
    
    A = NaN(options.lat_dim,options.lat_dim,options.T);
    
    for tt = 1:options.T
        
        if(tt <= T_switch_dyn)
            A(:,:,tt) =   V1 * diag(eig_vals(:,tt)) * inv(V1);
        else
            A(:,:,tt) =   V2 * diag(eig_vals(:,tt)) * inv(V2);
        end
        
        A(:,:,tt) = real(A(:,:,tt));
    
    end

    options.A = A;
    
end

% Obtained a binned-version of the ground-truth dynamics matrices for a
% fair comparison to the estimated dynamics matrices. This step is
% necessary since our simulation step size (dt) does not match the bin-size
% used to compute residuals.

[options.binned.dynamics.A] = compute_binned_groundtruthdynamics(options.A, ...
    options.dt, options.bin_size);

%% Generate initial conditions, latent noise, and mean driving input

% Specify the mean of the initial condition for each trial. (Typically set
% to zero in these simulations)

options.x0 = zeros(options.lat_dim, options.num_trials);


options.tau_input = [1;1;1].*options.dt;
% Specify the covariance of the latent noise
switch options.observation_model
    
    case 'poisson'
        
        sigma_noise_var = 12500; % This value is appropriate only for poisson observations 
        options.sigma_input_noise = sigma_noise_var.*diag(ones(options.lat_dim,1));

    case 'gaussian'
        
        % For the gaussian observation models, we either keep the latent
        % noise fixed over time OR vary it across time by switching it at
        % the time of the switch of the dynamics
        
        if(~options.time_varying_latent_noise)
            
            sigma_noise_var = 0.5;
            options.sigma_input_noise = sigma_noise_var.*diag(ones(options.lat_dim,1));

        else
            
            T_switch_noise = T_switch_dyn; % time at which you want to switch the noise statisticw
            sigma_noise_varpre = 0.5; %(pre-switch)  
            sigma_noise_varpost = 0.3; % post-switch
            options.sigma_input_noise = cat(3,repmat(sigma_noise_varpre.*diag(ones(options.lat_dim,1)),...
                [1 1 T_switch_noise]),repmat(sigma_noise_varpost.*diag(ones(options.lat_dim,1)),...
                [1 1 options.T - T_switch_noise]));
        end

end

% Specify the auto-correlation (integration) time constants of the generative 
% process that produces the latent noise. The value is set to 1 (scaled by dt)
% if you want white gaussian noise (i.e temporally uncorrelated noise),
% else > 1 (for temporally correlated latent noise)


% Generate the latent noise sequence for the simulation
[noise_seq] = generateCorrelatedInputNoise(options.sigma_input_noise, options.tau_input,...
    options.T,options.num_trials,options.dt, rseed);

% Generate the mean input drive
options.mean_input = [2;2;2];

% Noisy inputs (u_t) defined as the sum of input drive and noise (which may
% be either temporally correlated or uncorrelated depending on the setting
% of tau_input.

inputs = options.mean_input + noise_seq; 

% Use the augmented dynamics (if noise is temporally correlated) to obtain
% the steady-state value of the covariance of initial condition using the
% lyapunov equation.

[C_ss0, A_aug] = compute_augmentedLDS_correlatedInputNoise(options.A(:,:,1),...
    options.dt, options.tau_input./options.dt, options.sigma_input_noise);
C_ss0_xx = real(C_ss0(1:options.lat_dim,1:options.lat_dim));
options.sigma0 = 0.5*(C_ss0_xx + C_ss0_xx);

%% Setting observation model

% Generate observation matrix
C_pre = RandOrthMat(options.obs_dim); 
options.C = repmat(C_pre(:,1:options.lat_dim),[1 1 options.T]);

% Generate baseline mean and observation noise covariance 
% These parameters vay depending on whether you use gaussain observations 
% or poisson observations 
switch options.observation_model
    
    case 'poisson'
        
        options.d = 0.7*rand(options.obs_dim, 1) - 4.8; % baseline
        options.R = [];
        
        
    case 'gaussian'
        options.d =  8*rand(options.obs_dim, 1) - 2; % baseline
        options.R = diag(0.01.*rand(options.obs_dim,1));
        
end

%% Generate data for fitting

% Simulate Model
[states] = simulateLDS(options, inputs, rseed);

% Process simulated observations to obtain binned responses

data.response = states.y;
data.time = states.time;
use_sqrt_transform = false;
use_rates = false;
data_binned = binSpikeCounts(data, options.bin_size,(1./options.dt),use_sqrt_transform,use_rates);

% Compute residuals
switch options.observation_model
    
    % square-root transform spike counts

    case 'poisson'

        sqrt_transformed_data_binned = sqrt(data_binned.response);
        
        % Compute lowD residuals for computing the metrics that measure the
        % extent of variability. We reduce the dimensionality of the observations 
        % since these variability metrics scale with the number of dimensions.
        % Therefore, a fair comparison with the data, requires a match in
        % diemnsionality (i.e 20 aligned dimensions)
        
        n_dims_lowD = 20;
        [coeff,~,~,~,explained,mu] = pca(mean(sqrt_transformed_data_binned,3)');
        sqrt_transformed_data_binned_lowD = reshape(coeff(:,1:n_dims_lowD)'*(sqrt_transformed_data_binned(:,:) - mu'),[n_dims_lowD, ...
            size(sqrt_transformed_data_binned,2) size(sqrt_transformed_data_binned,3)]);
        residuals_lowD_binned = sqrt_transformed_data_binned_lowD - mean(sqrt_transformed_data_binned_lowD,3);
        [l0,l1,r,p_var] = compute_prop_ratio(residuals_lowD_binned);
        
        % Print the variability metrics
        fprintf(' l0 = %.2f\n', l0);
        fprintf(' l1 = %.2f\n', l1);
        fprintf(' pvar = %.2f\n', p_var);
        
        
        % Compute actual residuals used for fitting
        residuals_binned = sqrt_transformed_data_binned - mean(sqrt_transformed_data_binned,3);
        
    case 'gaussian'
        
        residuals_binned = data_binned.response - mean(data_binned.response,3);
       
        % Print the variability metrics
        [l0,l1,r,p_var] = compute_prop_ratio(residuals_binned);
        fprintf(' l0 = %.2f\n', l0);
        fprintf(' l1 = %.2f\n', l1);
        fprintf(' pvar = %.2f\n', p_var);
        
end

%% Running the residual dynamics pipeline

% Specify the various settings and hyper-parameters of the pipeline
opts_analysis.doFinalFit = true;
opts_analysis.doSmoothCV = true;

opts_analysis.hankel_order = 5;
opts_analysis.hankel_cv.hankel_ranks = [1 2 3 4 5];
opts_analysis.hankel_cv.n_cv = 20;
opts_analysis.hankel_cv.criterion = 'min';
opts_analysis.hankel_cv.type= 'combined';

opts_analysis.lag_cv.n_cv = 1;
opts_analysis.lag_cv.n_folds = 10;
opts_analysis.lag_cv.grid_type = 'uniform';
opts_analysis.lag_cv.do_display = true;
switch opts_analysis.lag_cv.grid_type
    
    case 'uniform'
        opts_analysis.lag_cv.iv_lag =  [1 2 3 4 5 8];
        opts_analysis.lag_cv.sub_dim = [1 2 3 4 5];
        
    case 'random'
        opts_analysis.lag_cv.numPars = 20;
        opts_analysis.lag_cv.max_lag = 10;
        opts_analysis.lag_cv.max_dim = 6;
             
end

opts_analysis.smooth_cv.alpha_all = [1e-3 1e-1 1e2 1e4 1e5 1e6 1e8];
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

% Fit the pipeline
time_labels = ones(length(data_binned.time),1);
trial_labels = ones(options.num_trials,1);
[analysis] = run_pipeline(residuals_binned, data_binned.time, data_binned.time,...
    time_labels, trial_labels, opts_analysis, rseed);

analysis.sim_opts = options;
%% Compute OLS fits (for comparison to 2SLS) for the gaussian observation simulations

if(strcmp(options.observation_model, 'gaussian'))

    unique_trial_labels = unique(trial_labels);
    num_conds = length(unique_trial_labels);
    
    % Fit the ordinary least squares dynamics matrix to residuals
    for icond = 1 : num_conds
        
        idxs_to_fit = analysis.result.fit_pars.lag{icond} + 1 : size(residuals_binned,2); 
        
        [data_fit] = applydynamicsSubspaceProjection(residuals_binned(:,idxs_to_fit,...
            trial_labels == unique_trial_labels(icond)),...
            analysis.result.U_dyn(:,1:analysis.result.fit_pars.dim{icond}));
        
        A_ols = cell(1, length(analysis.result.fit_pars.alpha));
        
        for ialign = 1: length(analysis.result.fit_pars.alpha)
            [A_ols{ialign}] = estimateTVRegularizedOLS(data_fit(:,2:end,:),...
                data_fit(:,1:end-1,:),analysis.result.fit_pars.alpha{icond}(ialign));
        end
        
        A_ols = cat(2,A_ols{:});
        
        [A_ols_stats] = computeDerivedDynamicsQuantities_v2(A_ols,'abs','right', ones(size(A_ols,3),1),...
            analysis.result.final_model{icond}.time(2) - analysis.result.final_model{icond}.time(1));

        stat_names = fieldnames(A_ols_stats);
        for ff = 1: length(stat_names)

            analysis.result.final_model_ols{icond}.(stat_names{ff}) = A_ols_stats.(stat_names{ff});

        end
        analysis.result.final_model_ols{icond}.A = A_ols;

    end 
    
end
%% Plot the results 

cond_idx_to_plot = 1;

% You only need this to reproduce plots if the analysis variable 
% is already stored to disk (see data subfolder?)

if(strcmp(options.observation_model,'poisson'))
    base_plot_pars.marker_style = '^';
elseif(strcmp(options.observation_model,'gaussian') & ~options.time_varying_latent_noise)
    base_plot_pars.marker_style = 'o';
elseif(strcmp(options.observation_model,'gaussian') & options.time_varying_latent_noise) 
    base_plot_pars.marker_style = 's';
end

base_plot_pars.col_map_type = 'seq';
if(cond_idx_to_plot == 1)
    base_plot_pars.col_map_cols = 'GnBu';
else
    base_plot_pars.col_map_cols = 'OrRd';
end
base_plot_pars.markersize = 4;
base_plot_pars.linewidth = 2;
base_plot_pars.do_plot_markers = true;
base_plot_pars.line_style = {'-'};
base_plot_pars.marker_face_color = [1 1 1];
base_plot_pars.do_plot_ols = false;
plot_evsvpanel(analysis, cond_idx_to_plot, base_plot_pars)

%% Saving the final results to disk

% Add generated data to variable that needs saving
analysis.data = data_binned;
analysis.residuals = residuals_binned;
analysis.data.trial_labels = trial_labels;

if(do_save_vars)
    file_name = sprintf('%s%s%s%d%s%d%s%d%s','SimulatedTVLDS',options.observation_model,'_binSize=',options.bin_size,'_modelID=',modelID,'_sigmaNoise=',sigma_noise_var,'.mat');
    save(fullfile(save_file_path, file_name),'analysis','-v7.3');
end
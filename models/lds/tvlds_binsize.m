%{
---------------------------------------------------------------------------
This script produces the simulations of a time-invariant dynamical system 
with gaussian observations and fits residual dynamics to observations binned
using different bin sizes. This is done in order to understand how the bin-size 
parameter contributes to bestiamtor bias. See Supplementary Figure 9.
---------------------------------------------------------------------------
Author: Aniruddh Galgali
%}

clearvars -except DIRS
clc

%% Define the basic parameters of LDS simulations

rseed = 2461;
options.rseed = rseed;
rng(rseed)
save_fig_path = './figures/simulations';

%%

options.dynamics_model = 'continuous';
options.observation_model = 'poisson';
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

%% Generate normal dynamics matrix

options.dyn_type = 'normal';
eig_vals = [-2;-4;-6];
eig_vals = eig_vals.*ones(options.lat_dim,options.T);
[options.A] =  generateDynamicsMatrices(options.lat_dim,options.dt,...
    options.T,eig_vals,options.dyn_type,rseed);

%% Specify the initial conditions, input drive and correlation time-scale of latent noise

% Specify the mean of the initial condition for each trial. (Typically set
% to zero in these simulations)
options.x0 = zeros(options.lat_dim, options.num_trials);

% Specify the auto-correlation (integration) time constants of the generative 
% process that produces the latent noise. The value is set to 1 (scaled by dt)
% if you want white gaussian noise (i.e temporally uncorrelated noise),
% else > 1 (for temporally correlated latent noise)
options.tau_input = [1;1;1].*options.dt;   % integration time constant of input noise

% Generate the mean input drive
options.mean_input = [2;2;2];
%%
% Generate observation matrix
C_pre = RandOrthMat(options.obs_dim); 
options.C = repmat(C_pre(:,1:options.lat_dim),[1 1 options.T]);

% Generate baseline mean and observation noise covariance 
% These parameters vay depending on whether you use gaussain observations 
% or poisson observations 
switch options.observation_model
    
    case 'poisson'
        
        options.d = 0.7*rand(options.obs_dim, 1) - 4.8; % baseline
        options.R = diag(0.5.*rand(options.obs_dim,1));
        
        
    case 'gaussian'
        options.d =  8*rand(options.obs_dim, 1) - 2; % baseline
        options.R = diag(0.05.*rand(options.obs_dim,1));

        
end


%% Generate input noise, initial noise and mean driving input

% Final fit hyperparameters of the residual dynamics. Here, we don't optimize
% the hyperparams using cross-validation, since the objective of the
% analyses are somewhat different. Nonetheless, we choose sensible values
% for these hyperparameters that match the ground truth of the simulation.

hankel_order = 5;
opts_final.do_cross_validate_noise_model = false;
opts_final.order_across_alignment = true;
opts_final.hankel_order = hankel_order;
opts_final.hankel_rank = {options.lat_dim}; % hardcoded this
opts_final.dim = {options.lat_dim}; % hardcoded this
opts_final.alpha = {1e6};

% bin-size parameters
reference_bin_size = 40; % in ms
bin_sizes = [2 3 5 10 15 30 40 60];
opt_lags = round(max(bin_sizes)./bin_sizes);

% Latent variances to test
latent_noise_variances = logspace(-0.5,0.7,6);

result_final = cell(length(latent_noise_variances), length(bin_sizes));
options = repmat(options,length(latent_noise_variances),1);

% additional params that specify how to preprocess simulate observations
use_sqrt_transform = false;
use_rates = false;

for ivar = 1: length(latent_noise_variances)
    
    % Setting the value of the latent noise covariance
    options(ivar).sigma_input_noise = latent_noise_variances(ivar).*diag(ones(options(ivar).lat_dim,1));

    % Generating the latent noise sequence
    [noise_seq] = generateCorrelatedInputNoise(options(ivar).sigma_input_noise, options(ivar).tau_input,...
        options(ivar).T,options(ivar).num_trials,options(ivar).dt, rseed);

    % Noisy inputs (u_t) defined as the sum of input drive and noise (which may
    % be either temporally correlated or uncorrelated depending on the setting
    % of tau_input.
    
    inputs = options(ivar).mean_input + noise_seq; 

    % Use the augmented dynamics (if noise is temporally correlated) to obtain
    % the steady-state value of the covariance of initial condition using the
    % lyapunov equation.
    
    [C_ss0, A_aug] = compute_augmentedLDS_correlatedInputNoise(options(ivar).A(:,:,1), options(ivar).dt, ...
        options(ivar).tau_input./options(ivar).dt, options(ivar).sigma_input_noise);
    C_ss0_xx = real(C_ss0(1:options(ivar).latent_dim,1:options(ivar).latent_dim));
    options(ivar).sigma0 = 0.5*(C_ss0_xx + C_ss0_xx);
    
    % Simulate Model & binning/computing residuals

    [states] = simulateLDS(options(ivar), inputs, rseed);
    data.response = states.y;
    data.time = states.time;
    data_binned = binSpikeCounts(data, options.bin_size,(1./options.dt),use_sqrt_transform,use_rates);
    residuals_binned = data_binned.response - mean(data_binned.response,3);
    
    % Fitting residual dynamics
    for ib = 1 : length(residuals_binned)

        opts_final.lag = {opt_lags(ib)};

        [result_final{ivar,ib}] = fit_residual_dynamics(residuals_binned{ib}, ...
            data_binned{ib}.time, data_binned{ib}.time,...
            ones(1,length(data_binned{ib}.time)), ones(size(residuals_binned{ib},3),1), opts_final);

        
        fprintf('Done with bin-size = %d ms and variance=%.5f\n',bin_sizes(ib),latent_noise_variances(ivar));
    end

end

% Extracting eigenvalues and averaging them over time
ev_rescaled = cell(size(result_final));

for ivar = 1:size(result_final,1)
    for ib = 1:size(result_final,2)
        for icond = 1:length(result_final{ivar,ib}.final_model)
        
            [ev_rescaled{ivar,ib, icond}] = ...
                rebin_eigenvalue_dynamics(result_final{ivar,ib}.final_model{icond}.A, ...
                bin_sizes(ib)*options(ivar).dt, reference_bin_size*options(ivar).dt,'standard');
        end
    end
end
%%
% This bins the ground truth matrices that have been defined at a finer
% temporal resolution
[binned_A] = compute_binned_groundtruthdynamics(options(1).A, options(1).dt, reference_bin_size);
[~,~,eigenvalue_groundtruth] = my_eigenshuffle(binned_A,'abs');

%% Plotting rescaled eigenvalues as function of bin-size

figure; cols = flipud(cbrewer('seq','Greys',10));set(gcf,'Position',[1416  233 1140 393])
do_save_fig = true;

for idim = 1:3
    
    ah = subplot(1,options(ivar).lat_dim,idim); hold(ah); set(ah,'plotbox',[1 1 1]);
    ev_rescaled_temporalmean = cellfun(@(x) mean(x(idim,:)), ev_rescaled,'uni', false);
    ev_rescaled_temporalmean = abs(cell2mat(ev_rescaled_temporalmean));
    hh = plot(ah, bin_sizes, ev_rescaled_temporalmean,'-');
    leg_labels = cell(length(hh),1);
    for ihh = 1:length(hh)
       set(hh(ihh),'Color',cols(ihh,:));
       set(hh(ihh),'Parent',ah);
       leg_labels{ihh} = sprintf('%s%.3f','latent-var=',latent_noise_variances(ihh));
    end
    
    plot([0 bin_sizes(end) + 1],[mean(abs(eigenvalue_groundtruth(idim,:)),2)...
        mean(abs(eigenvalue_groundtruth(idim,:)),2)], '--k');
    
    set(ah,'xlim',[0 bin_sizes(end) + 1],'xtick',[0: 10: bin_sizes(end)],...
        'xticklabel',[0: 10: bin_sizes(end)],'ylim',[0 1.1],'ytick',[0:0.2:1],...
        'yticklabel',[0:0.2:1]);
    
    for ib = 1: length(latent_noise_variances)
        
        plot(bin_sizes, ev_rescaled_temporalmean(ib,:),'o','markerfacecolor',[1 1 1],...
            'markeredgecolor',cols(ib,:));
    end
    
    if(idim == 1)
        legend(leg_labels,'Location','southeast');
    end
    
    xlabel(ah,'bin-size (ms)');
    ylabel(ah,sprintf('%s%d', 'rescaled eigenvalue',idim))

end

if(do_save_fig)
    
    figName = sprintf('%s%s%s%d%s%s',options.observation_model,'TVLDS_BinSizeExpts_latentNoiseSweeps','_dtref=',reference_bin_size);
    export_fig (fullfile(save_fig_path,figName),'-painters','-transparent','-pdf')
    
end
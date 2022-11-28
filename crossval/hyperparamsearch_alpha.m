function [cv_smoothness,T_alpha] = hyperparamsearch_alpha(X,U,lag,sub_dim,trial_labels,time_labels,pars,rseed)
%{ Function performs grid search for different settings (grid search) 
%  of the alpha hyperparameter (smoothness in second stage of 2SLS).
% Input
% - X (3d tensor)- of size n_dim x n_times x n_trials - Typically the
% residual data matrix.
% - U (matrix) - of size n_obs_dim x n_dyn_dim, indicating the dynamics
%   subspace
% - lag (scalar) - the number of lags (instruments) for the first stage
%   of 2SLS (hyperparameter denoted by l)
% - sub_dim (scalar) - dimensionality of dyn. sub. for the first stage
%   of 2SLS (hyperparameter denoted by l)
% - time_labels (array) - of size n_times x 1, each element is a
%   numerical assignment of the time bin to a specific task epoch.
% - trial_labels (array) - of size n_trials x 1, each element is a
%   numerical assignment of the trial to a given task condition (for e.g
%   choice)
% - rand_state - matlab random seed
% - pars - hyperparameters and associated settings
% Output
% - cv_smoothness(struct) - containing the cross-validated fit errors and the
%   associated statistics
% - T_alpha (table) - contains the optimal alpha for the
%   specified fit error criteria
% Author - Aniruddh Galgali (Jul 2018, Heavily modified May 2020)
%}
unique_trial_labels = unique(trial_labels);
unique_time_labels = unique(time_labels);
nconds = length(unique_trial_labels);
nalign = length(unique_time_labels);
alpha_all = pars.alpha_all;

assert(length(lag) == nconds,'length of lag array should match number of residual condtions');
assert(length(sub_dim) == nconds,'length of lag array should match number of residual condtions');


cv_smoothness = cell(nconds,1);
T_alpha = cell(nalign,nconds);

for icond = 1:nconds

    cv_smoothness{icond} = cell(length(alpha_all),1);

    for ii = 1:length(alpha_all)

        [cv_smoothness{icond}{ii}] = cross_validate_smoothness(X(:,:,trial_labels == unique_trial_labels(icond)),...
            U(:,1:sub_dim(icond)), lag(icond), alpha_all(ii),time_labels,pars.n_cv,pars.n_folds,rseed);

    end
    if(pars.combine_align)
       err_test_mse = cellfun(@(x) x.mse_test_full, cv_smoothness{icond},'uni',false);
       err_train_mse = cellfun(@(x) x.mse_train_full, cv_smoothness{icond},'uni',false);
       err_test_pvar = cellfun(@(x) x.pvar_test_full, cv_smoothness{icond},'uni',false);
       err_train_pvar = cellfun(@(x) x.pvar_train_full, cv_smoothness{icond},'uni',false);
    else
       err_test_mse = cellfun(@(x) x.mse_test_align, cv_smoothness{icond},'uni',false);
       err_train_mse = cellfun(@(x) x.mse_train_align, cv_smoothness{icond},'uni',false);
       err_test_pvar = cellfun(@(x) x.pvar_test_align, cv_smoothness{icond},'uni',false);
       err_train_pvar = cellfun(@(x) x.pvar_train_align, cv_smoothness{icond},'uni',false);
    end

    

    for ialign = 1:nalign

        err_test_mse_foldmean = cell2mat(cellfun(@(x) mean(x(:,:,ialign),2)',err_test_mse,'uni',false));
        err_train_mse_foldmean = cell2mat(cellfun(@(x) mean(x(:,:,ialign),2)',err_train_mse,'uni',false));
        
        err_test_pvar_foldmean = cell2mat(cellfun(@(x) mean(x(:,:,ialign),2)',err_test_pvar,'uni',false));
        err_train_pvar_foldmean = cell2mat(cellfun(@(x) mean(x(:,:,ialign),2)',err_train_pvar,'uni',false));
        
        err_test_mse_runmean = mean(err_test_mse_foldmean,2);
        err_train_mse_runmean = mean(err_train_mse_foldmean,2);

        
        err_test_mse_runstd = std(err_test_mse_foldmean,0,2);
        err_train_mse_runstd = std(err_train_mse_foldmean,0,2);

        err_test_pvar_runmean = mean(err_test_pvar_foldmean,2);
        err_train_pvar_runmean = mean(err_train_pvar_foldmean,2);
        
        err_test_pvar_runstd = std(err_test_pvar_foldmean,0,2);
        err_train_pvar_runstd = std(err_train_pvar_foldmean,0,2);
        
        [~,min_idx_mse] = min(err_test_mse_runmean);
        uB_mse = err_test_mse_runmean(min_idx_mse) + err_test_mse_runstd(min_idx_mse);  
        [valid_idx_mse] = err_test_mse_runmean <= uB_mse;
        
        [~,max_idx_pvar] = max(err_test_pvar_runmean);
        uB_pvar = err_test_pvar_runmean(max_idx_pvar) - err_test_pvar_runstd(max_idx_pvar);  
        [valid_idx_pvar] = err_test_pvar_runmean <= uB_pvar;

        mse_test_valid_mean = err_test_mse_runmean(valid_idx_mse);
        mse_test_valid_std = err_test_mse_runstd(valid_idx_mse);
        mse_train_valid_mean = err_train_mse_runmean(valid_idx_mse);
        mse_train_valid_std = err_train_mse_runstd(valid_idx_mse); 
        alpha_mse = alpha_all(valid_idx_mse)';
                
        pvar_test_valid_mean = err_test_pvar_runmean(valid_idx_pvar);
        pvar_test_valid_std = err_test_pvar_runstd(valid_idx_pvar);
        pvar_train_valid_mean = err_train_pvar_runmean(valid_idx_pvar);
        pvar_train_valid_std = err_train_pvar_runstd(valid_idx_pvar); 
        alpha_pvar = alpha_all(valid_idx_pvar)';

        T_mse = table(mse_test_valid_mean,mse_test_valid_std,mse_train_valid_mean,mse_train_valid_std,alpha_mse);
        T_mse.Properties.VariableNames = {'mean_cv_test','std_cv_test','mean_cv_train','std_cv_train','alpha'};
        
        T_pvar = table(pvar_test_valid_mean,pvar_test_valid_std,pvar_train_valid_mean,pvar_train_valid_std,alpha_pvar);
        T_pvar.Properties.VariableNames = {'mean_cv_test','std_cv_test','mean_cv_train','std_cv_train','alpha'};
       
        
        T_alpha{ialign,icond}.T_mse = sortrows(T_mse,{'alpha','mean_cv_test'},{'descend','ascend'});
        T_alpha{ialign,icond}.T_pvar = sortrows(T_pvar,{'alpha','mean_cv_test'},{'descend','descend'});

    end

%     if(n_align == 1)
% 
%         T_alpha = T_alpha{1};
% 
%     end
    
end

end

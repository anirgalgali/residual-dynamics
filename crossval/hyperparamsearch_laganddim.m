function [cv_results,cv_opt,hyper_params] = hyperparamsearch_laganddim(data,hankel_order, hankel_rank, time_abs, trial_labels, time_labels, rand_state, pars)
%{ Function performs cross-validation for different settings (grid search) 
%  of the hyperparams (dimensionality and lag).
% Input
% - data (3d tensor)- of size n_dim x n_times x n_trials - Typically the
% residual data matrix.
% - hankel_order (scalar) - the hankel order (q)
% - hankel_rank (cell array) - where, weach element specifies the optimal hankel
%   rank for a specific task condition within the data specified by X
% - time_abs (array) - of size n_times x 1, idnicates the absolute times
%   for all time-bins in the data 
% - time_labels (array) - of size n_times x 1, each element is a
%   numerical assignment of the time bin to a specific task epoch.
% - trial_labels (array) - of size n_trials x 1, each element is a
%   numerical assignment of the trial to a given task condition (for e.g
%   choice)
% - rand_state - matlab random state
% - pars - hyperparameters and associated settings
% Output
% - cv_results(struct) - containing the cross-validated fit errors and the
%   associated statistics
% - cv_opt(struct) - contains the optimal lag and dim params for the
%   specified fit error criteria
% - hyperparams (table) - table of tested hyperparameter settings
% Author - Aniruddh Galgali (Jul 2018, Heavily modified May 2020)
%}
[hyper_params] = sample_hyperparams(pars);
cv_results = cell(size(hyper_params,1),1);

for ii = 1: size(hyper_params,1)
    
    n_discard = max(hyper_params.lag) - hyper_params.lag(ii);
    
    [cv_results{ii}] = cross_validate_laganddim(data, hankel_order,hankel_rank, ...
        hyper_params.lag(ii),hyper_params.dim(ii),time_abs,time_labels,trial_labels,...
        n_discard, pars.n_cv, pars.n_folds,rand_state);

   if(pars.do_display)
        
        fprintf('Done with set %d\n',ii);
   end
    
end

[cv_opt] = extract_optimal_lagdim(cv_results,hyper_params);

end

function [hyper_params] = sample_hyperparams(pars)


switch pars.grid_type
    
    case 'uniform'
        
        iv_lag = pars.iv_lag; 
        sub_dim = pars.sub_dim;

        hyp_grid = allcomb(iv_lag,sub_dim);
        hyper_params = array2table(hyp_grid,'VariableNames',{'lag';'dim'});
       
  
    case 'random'
        
        iv_lag = randi([pars.min_lag pars.max_lag],1,pars.numPars)';
        sub_dim = randi([pars.min_dim pars.max_dim],1,pars.numPars)';
        hyper_params = table(iv_lag,sub_dim,'VariableNames',{'lag';'dim'});



end

end

function [cv_opt] = extract_optimal_lagdim(cv_joint,hyper_pars_all)


err_test_mse = cellfun(@(x) x.mse_test,cv_joint,'uni',false);
err_test_pvar = cellfun(@(x) x.pvar_test,cv_joint,'uni',false);

[n_cv,n_folds,n_conds] = size(err_test_mse{1});

cv_opt = cell(n_conds,1);

for i_cond = 1:n_conds
    
    
    mse_test_cond_mean = cell2mat(cellfun(@(x) mean(x(:,:,i_cond),2)', err_test_mse,'uni',false));
    mse_test_cond_std = cell2mat(cellfun(@(x) std(x(:,:,i_cond),0,2)', err_test_mse,'uni',false));
    mse_test_cond_sem = (1./sqrt(n_folds)).* mse_test_cond_std;
    
    pvar_test_cond_mean = cell2mat(cellfun(@(x) mean(x(:,:,i_cond),2)', err_test_pvar,'uni',false));
    pvar_test_cond_std = cell2mat(cellfun(@(x) std(x(:,:,i_cond),0,2)', err_test_pvar,'uni',false));
    pvar_test_cond_sem = (1./sqrt(n_folds)).* pvar_test_cond_std;
    T_mse = cell(n_cv,1);
    T_pvar = cell(n_cv,1);
    
    for i_cv = 1: n_cv
       
        [~,min_idx_mse] = min(mse_test_cond_mean(:,i_cv));
        uB_mse = mse_test_cond_mean(min_idx_mse,i_cv) + mse_test_cond_sem(min_idx_mse,i_cv);  
        [valid_idx_mse] = mse_test_cond_mean(:,i_cv) <= uB_mse;
        
        mse_test_valid_mean = mse_test_cond_mean(valid_idx_mse);
        mse_test_valid_std = mse_test_cond_std(valid_idx_mse);
        lags_mse = hyper_pars_all.lag(valid_idx_mse);
        dims_mse = hyper_pars_all.dim(valid_idx_mse);
        num_pars_mse = lags_mse.*(dims_mse.^2);
        
        T_mse{i_cv} = table(mse_test_valid_mean,mse_test_valid_std,lags_mse,dims_mse,num_pars_mse);
        T_mse{i_cv}.Properties.VariableNames = {'mean_cv','std_cv','lags','dims','num_pars'};
        T_mse{i_cv} = sortrows(T_mse{i_cv},{'num_pars'}, {'ascend'});

        [~,max_idx_pvar] = max(pvar_test_cond_mean(:,i_cv));  
        lB_pvar = pvar_test_cond_mean(max_idx_pvar,i_cv) - pvar_test_cond_sem(max_idx_pvar,i_cv);   
        [valid_idx_pvar] = pvar_test_cond_mean(:,i_cv) >= lB_pvar;
        
        pvar_test_valid_mean = pvar_test_cond_mean(valid_idx_pvar);
        pvar_test_valid_std = pvar_test_cond_std(valid_idx_pvar);
        lags_pvar = hyper_pars_all.lag(valid_idx_pvar);
        dims_pvar = hyper_pars_all.dim(valid_idx_pvar);
        num_pars_pvar = lags_pvar.*(dims_pvar.^2);
        T_pvar{i_cv} = table(pvar_test_valid_mean,pvar_test_valid_std,lags_pvar,dims_pvar,num_pars_pvar);
        T_pvar{i_cv}.Properties.VariableNames = {'mean_cv','std_cv','lags','dims','num_pars'};
        T_pvar{i_cv} = sortrows(T_pvar{i_cv},{'num_pars'}, {'ascend'});
    end
    
    if(n_cv == 1)
        
        T_mse = T_mse{1};
        T_pvar = T_pvar{1};
       
    end
    
    cv_opt{i_cond}.mse_opt = T_mse;
    cv_opt{i_cond}.pvar_opt = T_pvar;

end

end
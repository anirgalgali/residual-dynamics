function [cv_lagdim] = cross_validate_laganddim(X,h_order,h_rank,iv_lag,sub_dim,time_abs,time_labels,trial_labels,n_discard,n_cv,n_folds,rand_state)
%{ Function that runs the cross-validation of the lag and dimensionality of 
%  the first stage of the 2SLS procedure.
% Input
% - X (3d tensor)- of size n_dim x n_times x n_trials - Typically the
% residual data matrix.
% - h_order (scalar) - the hankel order (q)
% - h_rank (cell array) - where, weach element specifies the optimal hankel
%   rank for a specific task condition within the data specified by X
% - iv_lag (scalar) - the number of lags (instruments) for the first stage
%   of 2SLS (hyperparameter denoted by l)
% - sub_dim (scalar) - dimensionality of the dynamics (hyperparameter denoted by d) 
% - time_abs (array) - of size n_times x 1, idnicates the absolute times
%   for all time-bins in the data 
% - time_labels (array) - of size n_times x 1, each element is a
%   numerical assignment of the time bin to a specific task epoch.
% - trial_labels (array) - of size n_trials x 1, each element is a
%   numerical assignment of the trial to a given task condition (for e.g
%   choice)
% - n_discard (scalar) - number of time-points to discard while coputing fit error
% - n_cv (scalar) - number of repeats of n-fold CV
% - n_folds (scalar) - number of folds of CV (tpyically = 5)
% - rand_state - matlab random seed
% Output
% - cv_lagdim (struct) - containing the cross-validated fit errors and the
%   associated statistics
%
% Author - Aniruddh Galgali (Jul 2018, Heavily Modified May 2020)
%}
rng(rand_state);
assert(iscell(h_rank),'hankel_rank should be a cell type, with number of elements equal to number of task conditions');
[n_dim,n_time,n_trials] = size(X);
assert(n_time == length(time_labels),'number of times do not match');

unique_trial_labels = unique(trial_labels);
unique_time_labels = unique(time_labels);

num_conds = length(unique_trial_labels);
num_align = length(unique_time_labels);

mse_test = NaN(n_cv,n_folds,num_conds);
mse_train = NaN(n_cv,n_folds,num_conds);

pvar_test = NaN(n_cv,n_folds,num_conds);
pvar_train= NaN(n_cv,n_folds,num_conds);


for i_cv = 1:n_cv
    
    cv_inds = crossvalind('Kfold',n_trials,n_folds);
    
    for ifold = 1:n_folds
        
        data_train = X(:,:,cv_inds~= ifold);
        data_test = X(:,:,cv_inds == ifold);
        
        trial_labels_train = trial_labels(cv_inds ~= ifold);
        trial_labels_test = trial_labels(cv_inds == ifold);
        
        [U_dyn] = compute_dynamics_subspace(data_train,trial_labels_train,...
            time_abs,time_labels,h_order,h_rank);
        
        unique_cond_labels_train = unique(trial_labels_train);
        unique_cond_labels_test = unique(trial_labels_test);
        
        for i_cond = 1:num_conds
            
            Xtrain_raw = data_train(:,:,trial_labels_train == unique_cond_labels_train(i_cond));
            Xtest_raw = data_test(:,:,trial_labels_test == unique_cond_labels_test(i_cond));
            
            time_idxs_predicted = [];
            X_pred_test = cell(num_align,1);
            X_pred_train = cell(num_align,1);
            
            for ialign = 1:num_align
                
                time_labels_align = find(time_labels == unique_time_labels(ialign));
                X_pred_test{ialign} = NaN(n_dim,length(time_labels_align),size(Xtest_raw,3));
                X_pred_train{ialign} = NaN(n_dim,length(time_labels_align),size(Xtrain_raw,3));
                
            end
            
            for ialign = 1:num_align
                
                time_labels_align = find(time_labels == unique_time_labels(ialign));
                
                Xtrain = Xtrain_raw (:,time_labels == unique_time_labels(ialign),:);
                Xtest = Xtest_raw(:,time_labels == unique_time_labels(ialign),:);
                
                % subtracting the mean (since data has been split)
                Xtrain_mean =  mean(Xtrain,3);
                
%                 Xtest_mean =  mean(Xtrain,3);
                
                Xtrain = bsxfun(@minus, Xtrain,Xtrain_mean);
                Xtest = bsxfun(@minus, Xtest, Xtrain_mean);
                
                [Xtrain_fit] = applydynamicsSubspaceProjection(Xtrain, U_dyn(:,1:sub_dim));
                [Xtest_dyn] = applydynamicsSubspaceProjection(Xtest, U_dyn(:,1:sub_dim));
                
                [betas_train] = getFirstStageCoefficients(Xtrain_fit,iv_lag);
                
                [X_pred_dyn_test,idxs_predicted] = predictFromPastLags(Xtest_dyn,betas_train);
                [X_pred_dyn_test] = applydynamicsSubspaceProjection(X_pred_dyn_test, U_dyn(:,1:sub_dim)');
                X_pred_dyn_test = bsxfun(@plus,X_pred_dyn_test,Xtrain_mean(:,idxs_predicted));
                
                
                [X_pred_dyn_train,idxs_predicted] = predictFromPastLags(Xtrain_fit,betas_train);
                [X_pred_dyn_train] = applydynamicsSubspaceProjection(X_pred_dyn_train, U_dyn(:,1:sub_dim)');
                X_pred_dyn_train = bsxfun(@plus,X_pred_dyn_train,Xtrain_mean(:,idxs_predicted));
                
                X_pred_test{ialign}(:,idxs_predicted,:) = X_pred_dyn_test ;
                X_pred_train{ialign}(:,idxs_predicted,:) = X_pred_dyn_train ;
                
                
                idxs_align_predicted = time_labels_align(idxs_predicted);
                idxs_align_predicted(1:n_discard) = [];
                time_idxs_predicted = cat(2,time_idxs_predicted,idxs_align_predicted);
                
            end
            
            X_pred_test_all = cat(2,X_pred_test{:});
            X_pred_train_all = cat(2,X_pred_train{:});
            
            model_res_train = Xtrain_raw (:,time_idxs_predicted,:) - X_pred_train_all(:,time_idxs_predicted,:) ;
            model_res_test = Xtest_raw(:,time_idxs_predicted,:) - X_pred_test_all(:,time_idxs_predicted,:) ;

            
            mse_train(i_cv,ifold,i_cond) = mean(sum(abs(model_res_train(:,:)).^2,1));
            mse_test(i_cv,ifold,i_cond) = mean(sum(abs(model_res_test(:,:)).^2,1));
            
            pvar_train(i_cv,ifold,i_cond) = percvar(Xtrain_raw(:,time_idxs_predicted,:),X_pred_train_all(:,time_idxs_predicted,:));
            pvar_test(i_cv,ifold,i_cond) = percvar(Xtest_raw(:,time_idxs_predicted,:),X_pred_test_all(:,time_idxs_predicted,:));
            
            
        end
        
    end
    
end
cv_lagdim.mse_train = mse_train;
cv_lagdim.mse_test = mse_test;
cv_lagdim.pvar_train = pvar_train;
cv_lagdim.pvar_test = pvar_test;

end
            

function [cv_smooth] = cross_validate_smoothness(X,U,iv_lag,alpha,time_labels,n_cv,n_folds,rand_state)
%{ Function that runs the cross-validation of smoothnes shyperparams (alpha) 
%  of the second stage of the 2SLS procedure. 
% Input
% - X (3d tensor)- of size n_dim x n_times x n_trials - Typically the
% residual data matrix.
% - U (matrix) - of size n_obs_dim x n_dyn_dim, indicating the dynamics
%   subspace
% - iv_lag (scalar) - the number of lags (instruments) for the first stage
%   of 2SLS (hyperparameter denoted by l)
% - time_labels (array) - of size n_times x 1, each element is a
%   numerical assignment of the time bin to a specific task epoch.
% - n_cv (scalar) - number of repeats of n-fold CV
% - n_folds (scalar) - number of folds of CV (tpyically = 5)
% - rand_state - matlab random seed
% Output
% - cv_smooth (struct) - containing the cross-validated fit errors and the
%   associated statistics
%
% Author - Aniruddh Galgali (First written Jul 2018, Heavily Modified May 2020)
%}


% Build first stage predictions
X = applydynamicsSubspaceProjection(X,U);
unique_time_labels = unique(time_labels);
numAlignments = length(unique_time_labels);
X2 = cell(1,numAlignments);
Y2 = cell(1,numAlignments);
time_idxs_predicted = [];
for ialign = 1: numAlignments
    
    time_labels_align = time_labels == unique_time_labels(ialign);
    Xres = X(:,time_labels_align,:);
    [betas_res] = getFirstStageCoefficients(Xres,iv_lag);
    [X_fs,idxs_predicted] = predictFromPastLags(Xres,betas_res);
        
    X2{ialign} = X_fs(:,1:end-1,:);
    Y2{ialign} = Xres(:,idxs_predicted(2:end),:);
    time_labels_align = find(time_labels_align);
    idxs_align_predicted = time_labels_align(idxs_predicted(1:end-1));
    time_idxs_predicted = cat(2,time_idxs_predicted,idxs_align_predicted);
    
end

X2 = cat(2,X2{:});
Y2 = cat(2,Y2{:});
time_labels_predictors = time_labels(time_idxs_predicted);

rng(rand_state);
[n_dim,~,n_trials] = size(X2);

mse_train_full = NaN(n_cv,n_folds);
mse_test_full = NaN(n_cv,n_folds);
pvar_train_full = NaN(n_cv,n_folds);
pvar_test_full = NaN(n_cv,n_folds);
mse_train_align = NaN(n_cv,n_folds,numAlignments);
mse_test_align = NaN(n_cv,n_folds,numAlignments);
pvar_train_align = NaN(n_cv,n_folds,numAlignments);
pvar_test_align = NaN(n_cv,n_folds,numAlignments);
mse_ratio_improv_test = NaN(n_cv,n_folds);
mse_ratio_improv_train = NaN(n_cv,n_folds);
% test_loss = NaN(n_cv,n_folds,numAlignments);
% train_loss = NaN(n_cv,n_folds,numAlignments);

for i_cv = 1:n_cv
    

    cv_inds = crossvalind('Kfold',n_trials,n_folds);
    
    for i_fold = 1:n_folds
        
        data_Ytest = Y2(:,:,cv_inds == i_fold);
        data_Ytrain = Y2(:,:,cv_inds ~= i_fold);
        
        data_Xtest = X2(:,:,cv_inds == i_fold);
        data_Xtrain = X2(:,:,cv_inds ~= i_fold);
        
        X_pred_test_filtered = cell(numAlignments,1);
        X_pred_train_filtered = cell(numAlignments,1);
        X_pred_test_fs = cell(numAlignments,1);
        X_pred_train_fs = cell(numAlignments,1);
        
        for ialign = 1:numAlignments
            
            time_labels_align = find(time_labels_predictors == unique_time_labels(ialign));
            X_pred_test_filtered{ialign} = NaN(n_dim,length(time_labels_align),sum(cv_inds == i_fold));
            X_pred_train_filtered{ialign} = NaN(n_dim,length(time_labels_align),sum(cv_inds ~= i_fold));
            X_pred_test_fs{ialign} = NaN(n_dim,length(time_labels_align),sum(cv_inds == i_fold));
            X_pred_train_fs{ialign} = NaN(n_dim,length(time_labels_align),sum(cv_inds ~= i_fold));
        end
        
    
        for ialign = 1:numAlignments
        
            time_labels_align = find(time_labels_predictors == unique_time_labels(ialign));
            
            Xtrain = data_Xtrain(:,time_labels_align,:);
            Ytrain = data_Ytrain(:,time_labels_align,:);
            
            [A] = estimateTVRegularizedOLS(Ytrain,Xtrain,alpha);
            
            Xtest = data_Xtest(:,time_labels_predictors == unique_time_labels(ialign),:);
            Xtest_pred = propagateTestData(Xtest,A);
            Xtrain_pred = propagateTestData(Xtrain,A);
                       
            X_pred_test_filtered{ialign} = Xtest_pred;
            X_pred_train_filtered{ialign} = Xtrain_pred;            
            
            X_pred_test_fs{ialign} = Xtest;
            X_pred_train_fs{ialign} = Xtrain; 
            
        end
        
        
        X_test_pred_filtered_all = cat(2,X_pred_test_filtered{:});
        X_train_pred_filtered_all = cat(2,X_pred_train_filtered{:});
        
        X_test_pred_fs_all = cat(2,X_pred_test_fs{:});
        X_train_pred_fs_all = cat(2,X_pred_train_fs{:});

        model_res_test_filtered = data_Ytest - X_test_pred_filtered_all ;
        model_res_train_filtered = data_Ytrain - X_train_pred_filtered_all ;
        
        model_res_test_fs = data_Ytest - X_test_pred_fs_all ;
        model_res_train_fs = data_Ytrain - X_train_pred_fs_all ;
  
        mse_train_full(i_cv,i_fold) = (1./(numel(model_res_train_filtered))).*(sum(model_res_train_filtered(:).^2));
        mse_test_full(i_cv,i_fold) = (1./(numel(model_res_test_filtered))).*(sum(model_res_test_filtered(:).^2));
        pvar_train_full(i_cv,i_fold) = percvar(data_Ytrain, X_train_pred_filtered_all);
        pvar_test_full(i_cv,i_fold) = percvar(data_Ytest, X_test_pred_filtered_all);
        
        mse_ratio_improv_test(i_cv,i_fold)  = mean(sum(abs(model_res_test_filtered(:,:)).^2,1))./ mean(sum(abs(model_res_test_fs(:,:)).^2,1));
        mse_ratio_improv_train(i_cv,i_fold)  = mean(sum(abs(model_res_train_filtered(:,:)).^2,1))./ mean(sum(abs(model_res_train_fs(:,:)).^2,1));
        
        for ialign = 1:numAlignments
            
            time_labels_align = find(time_labels_predictors == unique_time_labels(ialign));
            
            model_res_train_filtered_align = model_res_train_filtered(:,time_labels_align,:);
            mse_train_align(i_cv,i_fold,ialign) = (1./(numel(model_res_train_filtered_align))).*(sum(model_res_train_filtered_align(:).^2));
            
            pvar_train_align(i_cv,i_fold,ialign) = percvar(data_Ytrain(:,time_labels_align,:), X_train_pred_filtered_all(:,time_labels_align,:));
             
            model_res_test_filtered_align = model_res_test_filtered(:,time_labels_align,:);
            mse_test_align(i_cv,i_fold,ialign) = (1./(numel(model_res_test_filtered_align))).*(sum(model_res_test_filtered_align(:).^2));
            
            pvar_test_align(i_cv,i_fold,ialign) = percvar(data_Ytest(:,time_labels_align,:), X_test_pred_filtered_all(:,time_labels_align,:));

            
        end
       
  
    end
  
    
end

cv_smooth.mse_train_full = mse_train_full;
cv_smooth.mse_test_full = mse_test_full;
cv_smooth.pvar_train_full = pvar_train_full;
cv_smooth.pvar_test_full = pvar_test_full;
cv_smooth.mse_train_align = mse_train_align;
cv_smooth.mse_test_align = mse_test_align;
cv_smooth.pvar_train_align = pvar_train_align;
cv_smooth.pvar_test_align = pvar_test_align;
% cv_results.test_loss = test_loss;
% cv_results.train_loss = train_loss;
cv_smooth.mse_ratio_test = mse_ratio_improv_test;
cv_smooth.mse_ratio_train = mse_ratio_improv_train;
end

function [data_pred] = propagateTestData(data,A)
data_pred = NaN(size(data));
for tt = 1:size(A,3)
    
    data_pred(:,tt,:) =  A(:,:,tt)*squeeze(data(:,tt,:));
    
end

end

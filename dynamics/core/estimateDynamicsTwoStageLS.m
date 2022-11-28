function [A,varargout] = estimateDynamicsTwoStageLS(X,iv_lag,lambda,time_labels)
%{ 

This function performs a 2SLS (two-stage least squares) regression to 
obtain an "unbiased" parameter estimate of the dynamics matrix that maps 
X_{t} to X_{t+1}

Inputs : X - 3D data tensor of size n_dim x n_time x n_trials (time-series
             for which you want to estimate dynamics)
         
         iv_lag - instrumental variable lag (specifies the number of
         "past-lags" that need to be used as instruments for the first
         stage of the regression
         
         lambda - total variation regularization parameter for second stage
         regression.
        
         cond_labels - vector of size n_trials x 1 that specifies the condition 
         id corresponding to each trial. 
  
         time_labels = vector of n_time x 1 that specifies which alignment
         the time-series belongs to

Outputs: A - Time-varying dynamics matrix. 3D tensor of size
         n_dim x n_dim x n_time - 1
        
        varargout{1} - time indices corresponding to the estimated dynamics
        matrices

The fit of the dynamics matrix can be made separately for time-series
belonging to different alignments and codnitions. If you want to ignore this 
just pass a "vector of ones" as argument for cond_labels and time_labels.  

        Author : Aniruddh Galgali, Dec 2018
        Edited - Feb 2019 - the fit no longer includes the condition flags.
        Just do it separately for the two conditions

%} 

unique_time_labels = unique(time_labels);
numAlignments = length(unique_time_labels);
time_idxs_predicted = [];
time_idxs_predicted_stage1 = [];
A_align = cell(numAlignments,1);
betaFS = cell(numAlignments,1);
X_res_fs_all = cell(numAlignments,1);
if(numel(lambda) == 1)
    lambda_align = lambda.*ones(numAlignments,1);
else
    lambda_align = lambda;
end
for ialign = 1: numAlignments
    
    time_labels_align = find(time_labels == unique_time_labels(ialign));

    Xres = X(:,time_labels == unique_time_labels(ialign),:);
    [betas_res] = getFirstStageCoefficients(Xres,iv_lag);
    [X_res_fs,idxs_predicted] = predictFromPastLags(Xres,betas_res);
    [A_align{ialign},loss_align{ialign}] = estimateTVRegularizedOLS(Xres(:,idxs_predicted(2:end),:),X_res_fs(:,1:end-1,:),lambda_align(ialign));
        
     betaFS{ialign} = betas_res; 
     X_res_fs_all{ialign} = X_res_fs;

    idxs_align_predicted = time_labels_align(idxs_predicted(1:end-1));
    time_idxs_predicted = cat(2,time_idxs_predicted,idxs_align_predicted);
    time_idxs_predicted_stage1 = cat(2,time_idxs_predicted_stage1,time_labels_align(idxs_predicted(1:end)));
end

A = cat(3,A_align{:});
loss = cell2mat(loss_align);
betas = cat(3,betaFS{:}); 
varargout{1} = time_idxs_predicted;
varargout{2} = X_res_fs_all;
varargout{3} = time_idxs_predicted_stage1;
varargout{4} = betas;
varargout{5} = loss;

end
function [X_pred,varargout] = predictFromPastLags(X,betas)

%%
%{ 

% This function uses the beta coefficients estimated in the first
% stage regression to predict the next time point from the past

Inputs : X - 3D data tensor of size n_dim x n_time x n_trials
         betas - Time-varying beta coefficients. 3D tensor of size
         n_dim*iv_lag x n_dim x n_time_betas

Outputs: X_pred :  3D data tensor of size n_dim x n_time x n_time_betas (We
         can start predicting only from the (l + 1)th time point onwards since we
         have to use the past 'l' lags for prediction.
         
        varargout : arg1 - time_indices of predicted data points
    
     Author : Aniruddh Galgali, Dec 2018
%}


X = permute(X,[3 1 2]);

[n_trials,n_dim,n_time] = size(X);

[n1,n2,n_time_betas] = size(betas);

assert(n_dim == n2,'betas are of incorrect size')

lag = n1/n2;

idxs_predicted = lag + 1 :  n_time;

assert(length(idxs_predicted) == n_time_betas,'bleh')

X_pred = NaN(n_trials,n_dim,n_time_betas);

for tt = 1 : n_time_betas 
   
     X_past = X(:,:,idxs_predicted(tt)-1:-1:idxs_predicted(tt)-lag);
     X_pred(:,:,tt) = (X_past(:,:) * betas(:,:,tt));
    
end

X_pred = permute(X_pred,[2 3 1]);
varargout{1} = idxs_predicted;

end
function [betas, varargout] = getFirstStageCoefficients(X,iv_lag)
%%
%{ 

% This function computes the beta coefficients corresponding to the first
% stage of the regression.

Inputs : X - 3D data tensor of size n_dim x n_time x n_trials
         iv_lag - number of past lags to use as predictors

Outputs: betas - Time-varying beta coefficients. 3D tensor of size
         n_dim*iv_lag x n_dim x n_time
         

        Author : Aniruddh Galgali, Dec 2018
            
Edited - Feb 2019 - returns degrees of freedom and residuals sume of
squares as optional output arguments
%}
           
X = permute(X,[3 1 2]);
[~, n_dim, n_time] = size(X);
idxs_predicted = iv_lag + 1 : n_time;
betas = NaN(n_dim*iv_lag,n_dim,length(idxs_predicted));
RSS = NaN(length(idxs_predicted),1);
doF = NaN(length(idxs_predicted),1);
t_count = 0;    
for tt = iv_lag + 1 : n_time
    
    t_count = t_count + 1;
    
    yP = X(:,:,tt-1:-1:tt-iv_lag);
    yF = X(:,:,tt);
    betas(:,:,t_count) = yP(:,:)\yF(:,:) ;
    
    RSS(t_count) = sum((yF(:) - vec(yP(:,:)*betas(:,:,t_count))).^2);
    doF(t_count) = numel(yF(:)) - numel(betas(:,:,t_count));
end


varargout{1} = idxs_predicted;
varargout{2} = RSS;
varargout{3} = doF;
end
function [pvar,varargout] = compute_variance_explained_by_projection(X, U, time_mask, trial_mask)
% Function to compute varaince explained by a projection matrix
% Input 
% - X (3d tensor) - of size n_dim x n_times x n_trials - typically a data
%  matrix containing either single trials or residuals.
% - U(matrix) - of size n_dim x p, where p is the dimenisonality of the
%   subspace defiend by the columns if U.
% - time_mask (boolean array) - indicating the times in X to consider for
%   the projection
% - trial_mask (boolean array) - indicating the trials in X to consider for
%   the projection
% 
% Outputs
% - pvar (scalar) - percent variance explained
% - varargout
%  - X_masked_msub - projection of X into subspace defined by U
%  - origVar - original variance explained in the space of X.
% Author: Aniruddh Galgali
if(isempty(trial_mask))
   trial_mask = true(size(X,3),1); 
end

if(isempty(time_mask))
   time_mask = true(size(X,2),1); 
end

X_masked = X(:,time_mask,trial_mask);
X_masked_msub = bsxfun(@minus, X_masked(:,:), mean(X_masked(:,:),2));
origVar = trace(cov(X_masked_msub'));
pvar =  sum(diag(U'*cov(X_masked_msub')* U)./ origVar);
varargout{1} = reshape(X_masked_msub, [size(X_masked)]);
varargout{2} = origVar;
end
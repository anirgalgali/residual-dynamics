function [X_out] = applydynamicsSubspaceProjection(X, U)
%{ A helper function that projects data in X into the column space defined
%  by the columns of U. Typically U is the dynamics subspace and X is the
%  data matrix
%  Input:
%  X - data matrix of size n_dim x n_times x n_trials (or n_conds)
%  depending on twhether data corresponds to single trials or condition
%  averages
%  U - projection matrix of size n_dim x n_dyn_dim, where n_dyn_dim is the 
%  dimensionality of the subspace
%  
%  Output
%  Xout - projected data (structured similary to X).
%}
% warning('The projection matrix U should be orthonormal')
assert(size(U,1) == size(X,1), 'projection matrix of incorrect dimension. The number of rows of U should match the first dimension of X')
X_out =  reshape(U'*X(:,:),[size(U,2) size(X,2) size(X,3)]);

end
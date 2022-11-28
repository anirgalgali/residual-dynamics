function [O_t,varargout] = computeObservabilityMatrixfromHankel(H_t,r)
%{ 

Computes a sequence of time-varying observability matrices O_t using the provided 
sequence of hankel matrices in H_t. The constructed O_t are of rank r(t), where
r(t) specifies the rank of O(t). The observability matrices are constructed
using what is known as a "balanced realization" in system identification
theory.

Inputs : H_t - 3D data tensor of size (n_dim*q) x (n_dim*q) x
                n_times_hankel
         r - 1D array of length 1 x n_times_hankel, which specifies the rank
             to be used to construct I_t from its corresponding H_t

Outputs: O_t - A cell array of length 1 x n_times_hankel. Each element of the 
               cell array specifies the obsevability matrix at time 't', which
               is a matrix of size (n_dim*q) x r(t) 
         OPTIONAL (in-order):
            1) "S" - A matrix of size n_times_hankel x n_dim_hankel, where
                     n_dim_hankel is the dimensionality of the hankel matrix
                     (= number of dimensions of data matrix x hankel_order
                     (q)). Rach column of S_t provides the singular values
                     of H_t at that time instant

        Author : Aniruddh Galgali, Dec 2017

%}

n_dim_hankel = size(H_t,1);
n_times_hankel = size(H_t,3);
assert(length(r) == n_times_hankel, 'need a rank defined for each hankel matrix')
O_t = cell(1,n_times_hankel);
S = NaN(n_times_hankel, n_dim_hankel);
for tt = 1:n_times_hankel
    
      [U,ss,~] = svd(H_t(:,:,tt),'econ');
      O_t{tt} = U(:,1:r(tt))*sqrt(ss(1:r(tt),1:r(tt))); %balanced realization
      S(tt,:) = diag(ss);
end

varargout{1} = S;
end
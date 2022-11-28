function [A,varargout] = estimateTVRegularizedOLS(Y,X,lambda)
    
%{ 

This function does a least squares regression between two adjacent time 
points of a time-series resulting in a time-varying regression parameter.
The regression naturally employs a regularization involving total-variation
type norm to ensure a "smoothness" between regression parameters corresponding
to adjacent time-points. The loss function that is minimized is as follows:  
\sum_{t} || Y_{t} - A_{t}*X_{t} ||_{F}^{2} + lambda*|| A_{t+1 - A_{t} ||_{F}^{2} 
The function vectorizes the above loss function and solves for the entire
series of A matrices in one go.

Inputs : Y - 3D data tensor of size n_dim x n_time x n_trials
         X - 3D data tensor of size n_dim x n_time x n_trials
         lambda - total variation regularization parameter

Outputs: A - Time-varying regulaized OLS coefficients. 3D tensor of size
         n_dim x n_dim x n_time
         

        Author : Aniruddh Galgali, Dec 2018
                 Modified : Apr 2019 (bug fixes)    


%} 
       

[p,T,K] = size(Y);


D = spdiags(ones(p.^2*(T-1),1)*[-1 1],[0 p.^2],p.^2*(T-1),p.^2*(T));
Y_v = NaN(p*K*T,1);
for tt = 1:size(X,2)
    
    B{tt} = sparse(kron(squeeze(X(:,tt,:))', eye(p)));
    Y_v((tt-1)*p*K + 1 : tt*p*K) = vec(Y(:,tt,:));
    
end

Ztilde = blkdiag(B{:});
if(p~=1)
    A = ((Ztilde'*Ztilde) + lambda*(D'*D)) \ (Ztilde'*Y_v(:)) ;
   loss = norm(Y_v(:) - Ztilde*A).^2 ;
else
    A = ((Ztilde*Ztilde') + lambda*(D'*D)) \ (Ztilde*Y_v(:)) ;
    loss = norm(Y_v(:) - Ztilde'*A).^2 ;
end

% loss = norm(Y_v(:) - Ztilde*A).^2 ;

A = reshape(A,[p p T]);

varargout{1} = loss;

end

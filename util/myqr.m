function [xx,rr] = myqr(kk)
% myqr orthogonalizes sets of vectors with qr decomposition
% 
% Inputs:
%  kk: raw vectors [ndim nvec]
%
% Output:
%  xx: orthogonalized vectors [ndim nvec]
%
% xx = myqr(kk)


% Number of regressors
[ndim,nvec] = size(kk);

% Orthogonalize
[qq,rr] = qr(kk);

% Keep signs consistent
sgn = sign(dot(qq(:,1:nvec),kk));
xx = qq(:,1:nvec).*repmat(sgn,[ndim 1]);

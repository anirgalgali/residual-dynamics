function [R, Q0, Q_iso_fixed,Q_aniso_fixed,Q_tv, Q_tv_iso] = computeLatentCovarianceParametersMLE(A_t, C, Y, X)

%{ This function computes the maximum-likelihood estimates of the latent and 
%  observation noise covariance given the estimates of the 
%   i) time-varying residual dynamics (A_t).
%   ii) the dynamics subspace (C)
%   iii) The observed residuals (Y)
%   iv) The output of the first stage of 2SLS i.e the latent residual
%   states (X)
%
%  The fucntion outputs several different estimates of the
%  latent/observation noise covariances conditioned on the type of
%  covariance estimates (i.e isotropic/anisotropic or
%  time-invariant/time-varying).
%
% Outputs
% R - observation noise covariance estimate
% Q0 - initial state noise covariance estimate
% Q_iso_fixed - isotropic, time-invaraint latent noise covariance estimate
% Q_aniso_fixed - anisotropic, time-invaraint latent noise covariance estimate
% Q_tv - anisotropic, time-varying latent noise covariance estimate
% Q_tv_iso - isotropic, time-varying latent noise covariance estimate
%
% Author : Aniruddh Galgali (Written : Apr 2019)
%}


[lat_dim, num_time_pts, num_trials] = size(X);

C_x =C*X(:,:);
R = diag(diag((1/(num_time_pts * num_trials)).*((Y(:,:) - C_x)*(Y(:,:) - C_x)')));
Q0 = (1/num_trials).*(squeeze(X(:,1,:))*squeeze(X(:,1,:))');

res_t = NaN(lat_dim,num_time_pts-1,num_trials);
Q_tv = NaN(lat_dim,lat_dim,num_time_pts-1);
Q_tv_iso = NaN(lat_dim,lat_dim,num_time_pts-1);
for tt = 1:size(A_t,3)
   
    res_t(:,tt,:) = squeeze(X(:,tt+1,:)) - A_t(:,:,tt)*squeeze(X(:,tt,:));
    Q_tv(:,:,tt) = (1/num_trials).*(squeeze(res_t(:,tt,:)) * squeeze(res_t(:,tt,:))');
    Q_tv_iso_t = (1/(num_trials*lat_dim))*(trace(squeeze(res_t(:,tt,:))*squeeze(res_t(:,tt,:))'));
    Q_tv_iso(:,:,tt) = Q_tv_iso_t.*eye(lat_dim);
end

Q_iso_fixed = (1/(num_trials*(num_time_pts - 1)*lat_dim))*(trace(res_t(:,:)*res_t(:,:)'));
Q_iso_fixed = repmat(Q_iso_fixed.*eye(lat_dim),[1 1 size(Q_tv,3)]);


Q_aniso_fixed = (1/(num_trials*(num_time_pts - 1)))*(res_t(:,:)*res_t(:,:)');
Q_aniso_fixed = repmat(Q_aniso_fixed,[1 1 size(Q_tv,3)]);
end
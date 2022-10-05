function [out_labels] = compute_analyticalsolutions_correlatedinputs(x, latent_dim, dt, V_dyn, V_inp, varargin)

% Computes the analytical solution to the 2SLS estimator for a simple, linear 
% dynamical system driven by correlated "input" noise.

% INPUTS

%   x - parameter vector (size : (3*latent_dim) x 1),

%   default value of latent_dim is always equal to 2 i.e the dynamical system
%   in each "area" is 2-dimensional. The individual elements of the vector x 
%   is as follows: x = [x1, x2, x3, x4, x5, x6], where,
%   x1 - eigenvalue of dimension-1 of the recorded area (used to define A)
%   x2 - eigenvalue of dimension-2 of the recorded area (used to define A)
%   x3 - eigenvalue of dimension-1 of the input area (used to define Phi)
%   x4 - eigenvalue of dimension-2 of the input area (used to define Phi)
%   x5 - variance of latent noise that drives the input area, along 
%        dimension 1 (used to define Q)
%   x6 - variance of latent noise that drives the input area, along 
%        dimension 2 (used to define Q)
%
%   latent_dim - latent dimensionality (always set to 2 anyway)
%   
%   V_eig - matrix (size : latent_dim x latent_dim), where each column 
%           corresponds to the eigenvector of the dynamics of the recorded
%           area. Typically set to eye(2)
%   V_inp - matrix (size : latent_dim x latent_dim), where each column 
%           corresponds to the eigenvector of the dynamics of the input
%           area. Typically set to eye(2)
%
%   varargin - (i) lag (l) used to compute 2sls solution

% OUTPUTS

%   out_labels- is a structure containing the
%   eigenvalues/time-constants/rotation-freqs that constitute the
%   input arguments to this function, as well as the eigenvalues/
%   time-constants/rotation-freqs of the analytical 2SLS estiamte.


if(nargin == 5)    
    lag = 1;
elseif(nargin == 6)
    lag = varargin{1};
end

diagA = x(1:latent_dim); % dynamics of observed area
A = V_dyn*diag(diagA)*inv(V_dyn);
A = real(A);

diagPhi = x(latent_dim+1: 2*latent_dim); % dynamics of inputs
Phi = V_inp*diag(diagPhi)*inv(V_inp);
Phi = real(Phi);

noiseCovs = x(2*latent_dim + 1 : 3*latent_dim); % latent noise covariance matrix
Q = diag(noiseCovs);

tau_dyn = convert_originaleigenvalues_to_timeconstants(abs(diagA), dt);
tau_noise = convert_originaleigenvalues_to_timeconstants(abs(diagPhi), dt);
rot_ang_dyn = uniquetol(abs(angle(diagA)),1e-6);
rot_freq_dyn =  rot_ang_dyn./(2*pi*dt);
rot_ang_noise = uniquetol(abs(angle(diagPhi )),1e-6);
rot_freq_noise =  rot_ang_noise./(2*pi*dt);


out_labels.eigval_dyn_dim1 = abs(diagA(1));
out_labels.eigval_dyn_dim2 = abs(diagA(2));
out_labels.eigval_noise_dim1 = abs(diagPhi(1));
out_labels.eigval_noise_dim2 = abs(diagPhi(2));


out_labels.rot_angle_dyn = rot_ang_dyn;
out_labels.rot_freq_dyn = rot_freq_dyn;
out_labels.rot_angle_noise = rot_ang_noise;
out_labels.rot_freq_noise = rot_freq_noise;

out_labels.tau_dyn_dim1 = tau_dyn(1);
out_labels.tau_dyn_dim2 = tau_dyn(2);
out_labels.tau_noise_dim1 = tau_noise(1);
out_labels.tau_noise_dim2 = tau_noise(2);
out_labels.log_latentnoise_dim1 = log10(noiseCovs(1));
out_labels.log_latentnoise_dim2 = log10(noiseCovs(2));
out_labels.lin_latentnoise_dim1 = noiseCovs(1);
out_labels.lin_latentnoise_dim2 = noiseCovs(2);


% Computing the analytical solution
[A_2sls] = compute_2sls(A, Phi, Q, lag);


% effective eigenvalues
[vv, dd] = eig(A_2sls);
[~,ss,~] = svd(A_2sls);

sval = diag(ss); % singular values
dd = diag(dd); % eigen values


% having a consistent alignment for the eigenvectors of the analytical
% soln.

fl_dotP = abs(real((vv' * V_dyn(:,1))));
if(fl_dotP(1) >= fl_dotP(2))
    idx_fl = [1 2];
else
    idx_fl = [2 1];
end

dd = dd(idx_fl);
vv = vv(:,idx_fl);

if(isreal(vv))
    S = squeeze(real(sum(vv.*V_dyn,1))) < 0; % making signs consistent
    vv(:,S) =  - vv(:,S);   
end

if(any(imag(dd) >= 1e-8))
    
    rot_angle = uniquetol(abs(angle(dd)),1e-6);
    rot_freq = rot_angle./(2*pi*dt);
else
   
    rot_angle = 0;
    rot_freq = 0;
end


eigEff = abs(dd);
tauEff = convert_originaleigenvalues_to_timeconstants(eigEff, dt);


out_labels.eigval_eff_dim1 = eigEff(1);
out_labels.eigval_eff_dim2 = eigEff(2);
out_labels.tau_eff_dim1 = tauEff(1);
out_labels.tau_eff_dim2 = tauEff(2);
out_labels.max_tau_eff = max(tauEff(1),tauEff(2));
out_labels.min_tau_eff = min(tauEff(1),tauEff(2));
out_labels.max_tau_dyn = max(tau_dyn(1),tau_dyn(2));
out_labels.min_tau_dyn = min(tau_dyn(1),tau_dyn(2));
out_labels.max_tau_noise = max(tau_noise(1),tau_noise(2));
out_labels.min_tau_noise = min(tau_noise(1),tau_noise(2));
out_labels.max_sval_eff = sval(1);
out_labels.min_sval_eff = sval(2);
out_labels.rot_angle_eff = rot_angle;
out_labels.rot_freq_eff = rot_freq;

end

function [A_2sls] = compute_2sls(A, Phi, Q, m)

%{ This function computes the analytical solution to the two stage least squares (2SLS). 
%  estimate of a linear dynamical system driven by correlated input noise, for a specific. 
%  lag 'l'.

%  For a discrete-time linear system defined as 
%  (i) x_{t+1} = A*x_{t} + \epsilon_{t}
%  (ii) epsilon_{t} = \Phi*epsilon_{t-1} + xeta_{t} 
%  where, epsilon_{t} describes a 'correlated; latent noise process. 

%  This function uses the augmented dynamical system defined by:
%  p(t+1) = A_aug*p_{t} + /nu_{t}, where p(t) = [x(t)' epsilon_{t}']'
%  to obtain the 2SLS estimate for a lag = l. See Supplementary Math Note B 

%  Inputs - A - matrix 
%           (latent_dim x latent_dim):
%           matrix specifying the dynamics of the "recorded" area. 
%           
%           Phi (matrix)
%           (latent_dim x latent_dim):
%           matrix specifying the dynamics of the "input" area.
%
%           Q (matrix)
%           matrix specifying the covariance of the noise that drives the
%           "input" area
%       
%           m - scalar
%           The lag parameter used for computation of A_2sls (default = 3)

%
%  Outputs - A_2sls - matrix
%            (latent_dim x latent_dim)
%            The analytical 2SLS estimate
%
% Author : Aniruddh Galgali
%}

% Define the augmented higher-dim dynamical system
latent_dim = size(A,1);
Q_aug = [zeros(latent_dim ,latent_dim) zeros(latent_dim ,latent_dim);...
    zeros(latent_dim , latent_dim) Q];

A_aug = [A eye(latent_dim); zeros(latent_dim,latent_dim) Phi];

% steady-state covariance of the augmented system at lag= 0
Css_p_0 = dlyap(A_aug,Q_aug);

% steady-state covariance of the augmented system at each lag l between 1
% and m
Css_p_l = cell(m+1,1);
for ii = 1:(m+1)
    Css_p_l{ii} = (A_aug^(ii)) * Css_p_0;
end

% Compute matrices P0, P1 and P2 for computation of analytical estimate
P0 = zeros(latent_dim*m,latent_dim*m); 
for ii = 1:m
    for jj = 1:m
    
        idx_row = (ii-1)*latent_dim + 1 : ii*latent_dim;
        idx_col = (jj-1)*latent_dim  + 1: jj*latent_dim;
        if(ii == jj)
            
            P0(idx_row,idx_col) = Css_p_0(1:latent_dim,1:latent_dim);
            
        elseif(jj > ii)
            
            P0(idx_row,idx_col) = Css_p_l{abs(ii-jj)}(1:latent_dim,1:latent_dim);
            
        else
            
            P0(idx_row,idx_col) = Css_p_l{abs(ii-jj)}(1:latent_dim,1:latent_dim)';
            
        end
    end
end

P1 = [];
P2 = [];

for ii = 1:m
   
    P1 = cat(2,P1,Css_p_l{ii}(1:latent_dim,1:latent_dim));
    P2 = cat(2,P2,Css_p_l{ii + 1}(1:latent_dim,1:latent_dim));
end

% Compute the analytical estimate
N = P2*inv(P0)'*P1';
D = P1*inv(P0)'*P1';

A_2sls = N/D;

end
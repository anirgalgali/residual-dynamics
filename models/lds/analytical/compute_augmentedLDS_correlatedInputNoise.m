function [C_ss0, A_aug, varargout] = compute_augmentedLDS_correlatedInputNoise...
                                     (A, dt, tau_noise_input, sigma_noise_input)

%{ This function casts an autonomous, continuous-time linear dynamical process
%  driven by correlated latent/input noise as a "higher-dimensional" augmented
% linear system driven by uncorrelated input noise. 

%  For a continuous-time linear system defined as the following two set of equations:
%  (i) x_dot = A*x_{t} + \epsilon_{t}, where \epsilon_{t} is a correlated 
% latent noise process determined by the equation below:
%  (ii) epsilon_{t} = \Phi*epsilon_{t-1} + xeta_{t} 

%  This function casts the above equation in terms of the discretized dynamics 
%  of a higher-dimensional state variable p(t) = [x(t)' epsilon_{t}']', where
%  p(t+1) = A_aug*p_{t} + /nu_{t}
%
%  Inputs - A - matrix (continuous-time dynamics!)
%           (latent_dim x latent_dim):
%           matrix specifying the dynamics matrix. 
%           
%           dt (scalar)
%           step-size of the dynamics
%
%           tau_noise_input - array
%           (latent_dim x 1)
%           time-constants that determine the correlation time scale of the
%           latent noise epsilon_{t}. These parameters are directly related 
%           to the dynamics matrix (\Phi, see above). Specified in units of
%           milliseconds!
%       
%           sigma_noise_input - matrix (continuous time representation!)
%           (latent_dim x latent_dim)
%           The covariance matrix of the latent noise process xeta_{t} (see
%           above)
%
%  Outputs - C_ss0 - matrix
%            (latent_dim x latent_dim)
%            The steady-state covariance of the initial conditions
%            
%            A_aug - matrix
%            (2*latent_dim x 2*latent_dim)
%            The dynamics matrix of the augmented higher-D linear system
%
%            varargout
%}

latent_dim = size(A,1);

% Discretizing the continuous-time dynamics matrix
Ad = discretize_dynamics_matrix(A,dt);

% Obtaining the discrete-time equivalent of the covariance matrix that
% spcifies the latent noise of the continuous time system
sigma_input_discrete = sigma_noise_input * (dt^2);

% Adding a correlation time-scale to the latent noise. The time constants
% of the auto-correlation (tau) determine the dynamics of the noise
% specified by matrix \Phi

[v_noise,d_noise] = eig(sigma_input_discrete);
phi_noise = eye(latent_dim) + v_noise * diag(-dt./(tau_noise_input.*dt))*v_noise';
Q_sub = v_noise * (2* d_noise * diag(dt./(tau_noise_input.*dt))) * v_noise';

for i_dim = 1:latent_dim
    if((tau_noise_input(i_dim)* dt) == dt)
        Q_sub(:,i_dim) = 0.5*Q_sub(:,i_dim);
    end
end

% Constructing the augmented, higher-dimensional system
Q_aug = [zeros(latent_dim ,latent_dim ) zeros(latent_dim ,latent_dim);...
    zeros(latent_dim , latent_dim) Q_sub];

A_aug = [Ad eye(latent_dim);zeros(latent_dim,latent_dim) phi_noise];
C_ss0 = dlyap(A_aug,Q_aug);

varargout{1} = Ad;
varargout{2} = sigma_input_discrete;
varargout{3} = Q_sub;

end
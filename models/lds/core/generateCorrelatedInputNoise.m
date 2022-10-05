function [noise_seq] = generateCorrelatedInputNoise(noise_cov, taus, num_time_steps, num_trials, dt, varargin)
%{
This function generates trajectories from an Ornstein-Uhlenbeck process, in
order to produce time-series that have a defined temporal autocorrrelation.  

Inputs : noise_cov - matrix/tensor: 
            (latent_dim x latent_dim) or (latent_dim x latent_dim x num_times)
            A matrix or tensor specifying the time-invariant or time-varying
            covariance matrix. This determines the covariance of the time
            series at lag = 0

         taus - array:
            (latent_dim x 1)
            Determines the time-scale of the auto-correlation for each
            dimension of the time series.
         
         num_time_steps - scalar:
            Number of time-points in the trajectory

         num_trials - scalar:
            Number of independent trials to simulate

         dt - scalar(double)
            The time-step of the simulation
        
        varargin
            (i) random_seed - seed for reproducibility
  
Outputs : latent_states - 3D tensor:
          (latent_dim x num_times x num_trials)
          The simulated latent states

Author : Aniruddh Galgali (Dec 2018)

%}

if(nargin == 5)
    rng('default')
elseif(nargin == 6)
    rand_state = varargin{1};
    rng(rand_state);
else
    error('invalid number of input arguments')
end

lat_dim = size(noise_cov,1);

noise_seq = zeros(lat_dim,num_time_steps,num_trials);
assert(~any(taus==0),'time-constants cannot be zero');

thetas = zeros(lat_dim,1);
sigmas = zeros(lat_dim,1);

if(ndims(noise_cov) == 2)
   noise_cov = repmat(noise_cov,[1 1 num_time_steps]); 
end

% The parameters that define the OU process
for i_dim = 1:lat_dim
    
    thetas(i_dim) = (1/taus(i_dim));
    sigmas(i_dim) = (1/taus(i_dim));
    
end

if(lat_dim > 1)
    
    thetas_t = cell(num_time_steps,1);
    sigmas_t = cell(num_time_steps,1);
    
    for tt = 1:num_time_steps
        [vv,dd] = eig(noise_cov(:,:,tt));
        thetas_t{tt} = vv*diag(thetas)*vv';
        sigmas_t{tt} = vv*sqrt(2*dd.*diag(sigmas))*vv';
    end
end

% Generate trajectories (using a multi-variate Ornstein-Uhlenbeck process)
for kk = 1: num_trials
    
    noise_seq(:,1,kk) = chol(noise_cov(:,:,1))'*randn(lat_dim,1);
    
    for tt = 1:num_time_steps - 1
        % The way the parameters are defined, one requires a hack to scale
        % the variance of the OU process when the dimensionality = 1. This
        % is a hack that works, but will need to be made cleaner.
        if(lat_dim == 1)
        
            d_eta = thetas.*(-1 * noise_seq(:,tt,kk) * dt) + ...
                sqrt(2*dt* noise_cov(:,:,tt) .* sigmas)*randn(lat_dim, 1);
        
        else
            
           
            d_eta = -1*thetas_t{tt}*(noise_seq(:,tt,kk).* dt) + sigmas_t{tt}.*sqrt(dt)*randn(lat_dim,1);

        end
        
        noise_seq(:,tt+1,kk) = noise_seq(:,tt,kk) + d_eta;
  
    end
            
end

for i_dim = 1:lat_dim
   
    if(taus(i_dim) == dt)
        noise_seq(i_dim,:,:) = (1/sqrt(2)).*noise_seq(i_dim,:,:);
    end
    
end
noise_seq = noise_seq(:,2:end,:); 
end
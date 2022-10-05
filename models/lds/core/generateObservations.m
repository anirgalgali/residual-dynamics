function [observed_states,observed_states_noiseless] = generateObservations(C, d, latent_states, R, mode, varargin)
%{
This function generates observations using trajectories obtained from a
time invariant/varying latent dynamical process. The observations are either 
generated using a gaussian noise model or a poisson noise model 

Inputs : C - matrix/tensor: 
            (obs_dim x latent_dim) or (obs_dim x latent_dim x num_times)
            A matrix or tensor specifying the time-invariant or time-varying
            observation matrix. 

         d - vector:
            (obs_dim x 1)
            A time/trial independent additive baseline
         
         latent_states - tensor:
            (latent_dim x num_times x num_trials)
            Latent states generated from a time-varying/time-invariant
            latent dynamical process

         R - matrix
            (obs_dim x obs_dim)
            Covariance matrix of observation noise. We do not specify this
            as a time-varying parameter.

         mode - string
             A string that specifies whether to simulate samples froma
             gaussian or poisson observation noise model
        
        varargin
            (i) random_seed - seed for reproducibility
  
Outputs : observed_states - 3D tensor:
          (obs_dim x num_times x num_trials)
          The simulated observed states (contaminated with obs. noise)

          observed_states_noiseless - 3D tensor:
          (obs_dim x num_times x num_trials)
          The simulated observed states (without obs. noise)
            

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

obs_dim = size(C,1);
num_time_steps = size(latent_states,2);
num_trials = size(latent_states,3);

observed_states = NaN(obs_dim, num_time_steps, num_trials);
observed_states_noiseless = observed_states;

if(isempty(d))
    d_off = 0;
else
    d_off = d;
end

for tt = 1: num_time_steps
    
    if(size(C,3) == num_time_steps) % if C has a time dimension
        observed_states_noiseless(:, tt, :) = C(:,:,tt) * squeeze(latent_states(:,tt,:)) + d_off; 
    else
        observed_states_noiseless(:, tt, :) = C * squeeze(latent_states(:,tt,:)) + d_off;
    end
    
end

switch mode
    
    case 'gaussian'
    
        assert(size(R,1) == obs_dim,'');
        for tt = 1: num_time_steps
            observed_states(:,tt,:) = squeeze(observed_states_noiseless(:,tt,:)) + chol(R)' * randn(obs_dim, num_trials);
        end
        
        
    case 'poisson'
        
        observed_states_noiseless = exp(observed_states_noiseless);
        for tt = 1: num_time_steps
            observed_states(:,tt,:) = poissrnd(squeeze(observed_states_noiseless(:,tt,:)));
        end
 
end

end
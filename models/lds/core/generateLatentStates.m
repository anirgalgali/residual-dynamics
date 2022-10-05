function [latent_states] = generateLatentStates(A, x_0, Q_0, u_t, mode, dt, varargin)

%{
This function generates latent states of a time invariant/varying latent
dynamical process determined by the dynamics matrix A. The function
generates trajectories either from a continuous-time or discrete-time
linear dynamical process

Inputs : A - matrix/tensor: 
            (latent_dim x latent_dim) or (latent_dim x latent_dim x num_times)
            A matrix or tensor specifying the time-invariant or time-varying
            dynamics matrix. 

         x_0 - matrix:
            (latent_dim x num_trials)
            The mean of the initial condition for each trial
         
         u_t - matrix/tensor:
            (latent_dim x num_times x num_trials) or (latent_dim x num_times)
            A matrix or tensor specifying the time-varying inputs for
            each trial

         mode - string
            A string that specifies whether to simulate sample using a
            'continuous-time' or 'discrete-time' dynamical process

         dt - scalar(double)
            The time-step of the simulation (relevant only for continuous
            time processes)
        
        varargin
            (i) random_seed - seed for reproducibility
  
Outputs : latent_states - 3D tensor:
          (latent_dim x num_times x num_trials)
          The simulated latent states

Author : Aniruddh Galgali (Dec 2018)

%}

if(nargin == 6)
    rng('default')
elseif(nargin == 7)
    rand_state = varargin{1};
    rng(rand_state);
else
    error('invalid number of input arguments')
end

assert(ndims(u_t) >= 2, 'need to provide a time-varyiing input for each trial')
assert(ndims(x_0) == 2, 'need to provide a mean initial condition for each trial')


latent_dim = size(A,1);
num_trials = size(x_0,2);
num_time_steps = size(A,3);

latent_states = zeros(latent_dim, num_trials, num_time_steps);

if(ismatrix(A))
    A = repmat(A,[1 1 num_time_steps]);
end


for kk = 1:num_trials

    s0 = chol(Q_0)'*randn(latent_dim,1);
    latent_states(:,kk,1) = x_0(:,kk) + s0;
end

for tt = 2:num_time_steps
    
    if(latent_dim>1)
        
        switch mode
            
            case 'continuous'
            
                s_dot = A(:,:,tt-1)*latent_states(:,:,tt-1) +  squeeze(u_t(:,tt-1,:));
                latent_states(:,:,tt) = latent_states(:,:,tt-1) + s_dot.*dt;
            
            case 'discrete'
                latent_states(:,:,tt) = A(:,:,tt-1)*latent_states(:,:,tt-1) +  squeeze(u_t(:,tt-1,:));
                
                
        end
                
    else
        
        switch mode

            case 'continuous'
                s_dot = (squeeze(A(:,:,tt-1))*latent_states(:,:,tt-1))' +  squeeze(u_t(:,tt-1,:));
                latent_states(:,:,tt) = latent_states(:,:,tt-1) + s_dot'.*dt;
                
            case 'discrete'
                
                latent_states(:,:,tt) = (squeeze(A(:,:,tt-1))*latent_states(:,:,tt-1))' +  squeeze(u_t(:,tt-1,:));
  
        end
                
    end
    
end

latent_states = permute(latent_states,[1 3 2]);

end
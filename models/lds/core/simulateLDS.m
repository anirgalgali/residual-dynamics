function [states] = simulateLDS(options, inputs, rseed)
%{
This function generates the data from latent-variable linear dynamical
system with either a time-invariant/time-varying dynamical process and a
gaussina/poisson observation process

Inputs : options - struct: 
            Contains all the relevant parameters of the simulation

         inputs - matrix/tensor:
            (latent_dim x num_times x num_trials) or (latent_dim x num_times)
            A matrix or tensor specifying the time-varying inputs for
            each trial
         
         rseed - rng:
            random seed for reproducibility
  
Outputs : states - struct:
          contains the simulated latent and observed states along with
          timing information

Author : Aniruddh Galgali (Dec 2018)

%}

[latent_states] = generateLatentStates(options.A, options.x0, options.sigma0, inputs, ...
    options.dynamics_model, options.dt, rseed) ;
[observed_states, rates] = generateObservations(options.C, options.d, latent_states,...
    options.R, options.observation_model, rseed);

    
states.x = latent_states;
states.y = observed_states;
states.yr = rates;
states.time = 0:options.dt:(options.T-1)*options.dt;

end

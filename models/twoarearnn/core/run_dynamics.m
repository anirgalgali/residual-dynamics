function [latent,rates,varargout] = run_dynamics(J, I_noise, I_stim, I0, x0, tau_nmda, dt, pars_i_o, shuffle_ppc_to_pfc, varargin)

%{ Runs the dynamics of the two-area RNN  
%  Input
%  - J - (4 x 4 matrix) recurrent connectivity matrix  2 dims per  area)
%  - I_noise (4  x num_trials x num_times) - OU noise as inputs to each population
%  - I_stim (4 x num_trials x num_times) - External input to each population
%  - I0 (scalar) - baseline input
%  - x0 (4 x num_trials) - Initial conditions of the pop. activation
%  - tau_nmda (scalar) - time constant of the nmda synapses
%  - dt (scalar) - time-step of simulation
%  - pars_i_o - struct containing parameters defining the act. function
%  - shuffle_ppc_to_pfc (logical) - True = shuffle ppc to pfc connection
%  - varargin:
%     1) trial_labels (n_trials x 1) - vector contianing the trial labels
%     determined using the unshuffled simulation. This argument is only
%     necessary when shuffle_ppc_to_pfc is True
%
%  Output
%  - latent - (4 x  num_trials x num_times) - activation of the 4
%  populations across both area. First two dimensions correspond to ppc.
%  - rates - (4 x  num_trials x num_times) - rates of the 4
%  populations across both area. Rates correspond to the latents passed
%  through the activation function
%  - input_to_pfc (2 x  num_trials x num_times)- input to pfc from ppc
%
%  Author: Aniruddh Galgali (Oct 2020)
%}
minArgs = 9;
maxArgs = 10;
narginchk(minArgs,maxArgs)

if(nargin == minArgs)
    assert(~shuffle_ppc_to_pfc, 'missing an additional input argument for shuffle analysis')
elseif(nargin == maxArgs)
    assert(shuffle_ppc_to_pfc, 'shuffle mode cannot be off')
    targ_dir = varargin{1};
    [~,~,targ_dir_labels] =  unique(targ_dir);
end

K = size(I_stim,2);
T = size(I_noise,3);
latent = NaN(size(J,1), K, T);
latent(:,:,1) = x0;
H = @(x) (pars_i_o.a.*x - pars_i_o.b)./(1 - exp(-pars_i_o.c.*(pars_i_o.a.*x - pars_i_o.b))); % input-output response

rates = NaN(size(J,1), K, T-1);
input_to_pfc = NaN(size(J,1)/2, K, T-1);
J_zero_fwd = J;
J_zero_fwd(3:4,1:2) = zeros(2,2);

for tt = 1 : T-1
    
    if(shuffle_ppc_to_pfc)
        
        
        weighted_input = J_zero_fwd*latent(:,:,tt); 
        un_targ_labels = unique(targ_dir_labels);
        
        for ii = 1:length(un_targ_labels)
            input_to_pfc_cond = J(3:4,1:2) * latent(1:2,targ_dir_labels == un_targ_labels(ii),tt);
            idx_perm = randperm(size(input_to_pfc_cond,2));
            input_to_pfc_cond_shuf = input_to_pfc_cond(:,idx_perm);
            input_to_pfc(:,targ_dir_labels == un_targ_labels(ii),tt) = input_to_pfc_cond_shuf;
        end
       
        weighted_input(3:4,:) = weighted_input(3:4,:) + input_to_pfc(:,:,tt);
        
        
    else
        
        weighted_input  = J*latent(:,:,tt);
        input_to_pfc(:,:,tt) = J(3:4,1:2)*latent(1:2,:,tt);
        
    end
    
    inp_curr = weighted_input +  I0 +  I_stim(:,:,tt) + I_noise(:,:,tt);
    
    s_dot = -(1/tau_nmda)*latent(:,:,tt) + (1 - latent(:,:,tt)).*(pars_i_o.gamma*H(inp_curr));
    latent(:,:,tt+1) = squeeze(latent(:,:,tt)) + s_dot.*dt;
    
    rates(:,:,tt) = H(inp_curr);
end
varargout{1} = input_to_pfc;
end
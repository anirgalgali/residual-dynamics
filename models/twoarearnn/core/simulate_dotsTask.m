function [data] = simulate_dotsTask(input_pars, net_pars, sim_pars)
%{ Simulates the dots task for a single two-area network confiuration with
%  the potential of also running shuffled simulations   
%  
%  Input
%  - input_pars - struct containing the parameters of the inputs
%  - net_pars - struct containing the parameters defining the network conn.
%  - sim_pars - struct containing other relevant simulation parameters
%
%  Output
%  - data - struct containing the activations of the network (both latents
%           and rates) and other relevant task variables. Additonally, data 
%           can contain a 'shuffled' field containing the responses of the
%           model when the input to pfc from pfc is shuffled across trials.
%
%
%  Author: Aniruddh Galgali (Oct 2020)
%}

rng(sim_pars.rseed);

time = 0:sim_pars.dt:(sim_pars.T-1)*sim_pars.dt;

i_noise = NaN(2*sim_pars.lat_dim , sim_pars.K, sim_pars.T+1);
i_noise(:,:,1) = zeros(2*sim_pars.lat_dim ,sim_pars.K);

if(length(input_pars.coh_set) > 1)
    dots_coh = randsample(input_pars.coh_set,sim_pars.K,true)';
else
    dots_coh = repmat(input_pars.coh_set,[sim_pars.K 1]);
end

dots_dir = randsample(input_pars.dir_stim_set,sim_pars.K,true)';
sign_coh = dots_dir.*dots_coh;

for tt = 2:sim_pars.T+1
    
    
    i_dot = (1/net_pars.tau_ampa)* (-squeeze(i_noise(:,:, tt-1)) + ...
        sqrt(net_pars.tau_ampa*(input_pars.sigma_noise^2))*randn(2*sim_pars.lat_dim,sim_pars.K));
    
    i_noise(:,:,tt) = i_noise(:,:,tt-1) + i_dot.*sim_pars.dt;
    
end
i_noise = i_noise(:,:,2:end);

% stimulus input

I_stim_curr = [(net_pars.J_ext* input_pars.mu0.* (1 + sign_coh./100))'; ...
    (net_pars.J_ext* input_pars.mu0.* (1 - sign_coh./100))' ; zeros(1,sim_pars.K); zeros(1,sim_pars.K)];


I_stim = cat(3,zeros(2*sim_pars.lat_dim,sim_pars.K,sim_pars.T_on),repmat(I_stim_curr,[1 1 sim_pars.T - sim_pars.T_on]));

if(isfield(input_pars,'perturb'))
    
    I_perturb = zeros(size(I_stim));
    I_perturb(input_pars.perturb.dim_mask,:,input_pars.perturb.times) = ...
        repmat(input_pars.perturb.magnitude*input_pars.perturb.direction,[1 sim_pars.K length(input_pars.perturb.times)]);
    I_stim = I_stim + I_perturb;

end

J_s = [net_pars.self_ppc net_pars.fbk_pfc_ppc; net_pars.fwd_ppc_pfc net_pars.self_pfc]; % [local recurrence in ppc feedback from pfc to ppc; feed-fwd from ppc to pfc local recurrence in pfc]
[J] = create_multiarea_weightMatrix(J_s, net_pars.J_t);

%IC
ics = input_pars.x0 + sqrt(net_pars.tau_ampa*(input_pars.sigma_noise^2))*randn(size(J,1),sim_pars.K);

[latent,rates,input_to_pfc] = run_dynamics(J, i_noise, I_stim, input_pars.I0,...
    ics, net_pars.tau_s, sim_pars.dt, net_pars.act_fn, false);


data.latents = latent;
data.rates = rates;
data.input_to_pfc = input_to_pfc;
data.I_noise = i_noise;
data.I_stim = I_stim;
data.x0 = ics;
data.task_variable.dots_coh = dots_coh;
data.task_variable.sign_coh = sign_coh;
data.task_variable.dots_dir = dots_dir;
data.time = time;
data.time_rel = data.time;
data.time_iev = 1.*ones(1,length(data.time_rel));
data.readout_pfc = squeeze(latent(3,:,end) - latent(4,:,end));
data.readout_ppc = squeeze(latent(1,:,end) - latent(2,:,end));

[data] = assign_condition_to_trials(data, sim_pars.area_to_decode);

if(sim_pars.do_shuffle)
    
    [latent_shuffled,rates_shuffled,input_to_pfc_shuffled] = run_dynamics(J, i_noise, I_stim, input_pars.I0,...
        ics, net_pars.tau_s, sim_pars.dt, net_pars.act_fn, sim_pars.do_shuffle, data.task_variable.targ_dir);
    
    shuffled.latents = latent_shuffled;
    shuffled.rates = rates_shuffled;
    shuffled.input_to_pfc = input_to_pfc_shuffled;
    shuffled.task_variable.dots_coh = dots_coh;
    shuffled.task_variable.sign_coh = sign_coh;
    shuffled.task_variable.dots_dir = dots_dir;
    shuffled.time = time;
    shuffled.time_rel = time;
    shuffled.time_iev = 1.*ones(1,length(data.time_rel));
    shuffled.readout_pfc = squeeze(latent_shuffled(3,:,end) - latent_shuffled(4,:,end));
    shuffled.readout_ppc = squeeze(latent_shuffled(1,:,end) - latent_shuffled(2,:,end));
    [data.shuffled] = assign_condition_to_trials(shuffled, sim_pars.area_to_decode);

    
else
    fprintf('Skipping shuffle analysis for this network\n')
end

end
function [data] = assign_condition_to_trials(data, area_to_decode)

switch area_to_decode
    
    case 'ppc'
        data.task_variable.targ_dir = sign(data.readout_ppc)';
    
    case 'pfc'
        
        data.task_variable.targ_dir = sign(data.readout_pfc)';
end

data.task_variable.correct = false(length(data.task_variable.targ_dir),1);
unique_dots_coh = unique(data.task_variable.dots_coh);

for ii = 1: length(unique_dots_coh)
    
    trial_idx = find(data.task_variable.dots_coh == unique_dots_coh(ii));
    
    if(unique_dots_coh(ii) == 0)
        for kk = 1: length(trial_idx)
            if(rand(1) >= 0.5)
                data.task_variable.correct(trial_idx(kk)) = true;
            else
                data.task_variable.correct(trial_idx(kk)) = false;
            end
            
        end
        
    else
        
        idx_correct = data.task_variable.targ_dir(trial_idx) ==  data.task_variable.dots_dir(trial_idx);
        trial_idx_correct = trial_idx(idx_correct);
        trial_idx_incorrect = trial_idx(~idx_correct);
        data.task_variable.correct(trial_idx_correct) = true;
        data.task_variable.correct(trial_idx_incorrect) = false;
        
    end
    
end

end
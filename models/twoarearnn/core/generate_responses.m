function [data] = generate_responses(input_pars, net_pars, sim_pars, obs_pars)
%{ Generates high-dimensional observations and computes residuals using the 
%  two area network simulations.
%  
%  Input
%  - input_pars - struct containing the parameters of the inputs
%  - net_pars - struct containing the parameters defining the network conn.
%  - sim_pars - struct containing other relevant simulation parameters
%  - obs_pars - struct containing parameters of the linear gaussian
%               osbervation model
%
%  Output
%  - data - struct containing the 'binned' activations of the network (both
%           latents and rates), the binned high-dimensional observations,
%           residuals and other relevant task variables. Additonally, data
%           can contain a 'shuffled' field containing the corresponding
%           quantities for the shuffled condition
%
%
%  Author: Aniruddh Galgali (Oct 2020)
%}    
    [data] = simulate_dotsTask(input_pars, net_pars, sim_pars);
    [data] = generate_binned_observations(data,obs_pars);
    [data] = extract_choice_and_time_modes(data);
    if(isfield(data, 'shuffled'))
       [data.shuffled] = generate_binned_observations(data.shuffled,obs_pars);
       [data.shuffled] = extract_choice_and_time_modes(data.shuffled);
    end
 
    %% Compute condition dependent residuals

    if(iscell(data))
        [data_res] = compute_residual_trajectories(data, obs_pars.residual_condition_vars);
    else
        [data_res] = compute_residual_trajectories(num2cell(data), obs_pars.residual_condition_vars);
    end
    
    if(isfield(data, 'shuffled'))
        if(iscell(data.shuffled))
            [data_shuffled_res] = compute_residual_trajectories(data.shuffled, obs_pars.residual_condition_vars);
        else
            [data_shuffled_res] = compute_residual_trajectories(num2cell(data.shuffled), obs_pars.residual_condition_vars);
        end
        data_shuffled_res = data_shuffled_res{1};
        data.shuffled.residuals = data_shuffled_res.response;
    end
    
    data_res = data_res{1};    
    data.residuals = data_res.response;
    

end

function [data_out] = generate_binned_observations(data,obs_pars)
    


    dt = data.time(2) - data.time(1);

    yDim = size(obs_pars.C_obs,1);
    xDim = size(data.latents,1);
    nBins = round(obs_pars.bin_size/(dt*1000));

    T    = floor(size(data.latents,3) / nBins);
    K = size(data.latents,2);

    fnames = fieldnames(data);

    c = cell(length(fnames),1);
    data_out = cell2struct(c,fnames);


    data_out.response = nan(yDim,T,K);
    data_out.latents = nan(xDim,T,K);
    data_out.rates = nan(xDim,T,K);
    data_out.time = nan(1,T);
    if(isfield(data_out,'time_rel'))
        data_out.time_rel = nan(1,T);
    end

    if(isfield(data_out,'time_rel'))
        data_out.time_iev = nan(1,T);
    end

    for tt = 1:T

        iStart = nBins * (tt-1) + 1;
        iEnd   = nBins * tt;
        
        data_out.latents(:,tt,:)= sum(data.latents(:,:, iStart:iEnd), 3);
        
        data_out.rates(:,tt,:)= mean(data.rates(:,:, iStart:iEnd), 3);
        data_out.time(tt) = median(data.time(iStart:iEnd));
        
        data_out.response(:,tt,:) = obs_pars.C_obs*squeeze(data_out.latents(:,tt,:)) + ...
            chol(nBins.*obs_pars.R)'*randn(yDim,K);
        
        if(isfield(data_out,'time_rel'))
            data_out.time_rel(tt) = median(data.time_rel(iStart:iEnd));
        end

        if(isfield(data_out,'time_iev'))
            data_out.time_iev(tt) = median(data.time_iev(iStart:iEnd));
        end
    end
    data_out.task_variable = data.task_variable;
    data_out.readout_ppc = data.readout_ppc;
    data_out.readout_pfc = data.readout_pfc;
    
    if(isfield(data_out,'shuffled'))
       data_out.shuffled = data.shuffled; 
    end
    if(isfield(data_out,'x0'))
        data_out = rmfield(data_out,'x0');
    end
    if(isfield(data_out,'I_stim'))
        data_out = rmfield(data_out,'I_stim');
    end
    if(isfield(data_out,'I_noise'))
        data_out = rmfield(data_out,'I_noise');
    end
end


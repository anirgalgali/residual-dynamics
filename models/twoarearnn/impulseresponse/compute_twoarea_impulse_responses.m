function [impulse_response_groundTruth,impulse_response_predicted, fit_result_all,perturb_dirs_latent,varargout] = ...
    compute_twoarea_impulse_responses_v2(input_pars,sim_pars, net_pars, obs_pars, areas, result_cv, perturb_pars, fit_pars)

%{ This function simulates the ground-truth impulse response and evaluates the predicted impulse respone for a given two-area network
%
%  INPUT
%  - input_pars (struct) - containing the parameters of the inputs
%  - net_pars   (struct) -  containing the parameters defining the network conn.
%  - sim_pars   (struct) - containing other relevant simulation parameters
%  - obs_pars   (struct) -  containing parameters of the linear gaussian
%               osbervation model
%  - areas (cell array) - size of (n_areas + 1) x 1, containing the string
%    labels that define each fit condition. The first n_areas correspond to 
%    each of the individual areas of the model, and the last element
%    corresponds to the string label associated with global dynamics.
%  - result_cv (struct) - structure containing the cross-validated model
%    information (such as the hyperparameters of the pipeline)
%  - perturb_pars (struct) - containing the perturbation params - such as
%    the perturbation times, perturbation amgnitudes, and perturbation areas.
%  - fit_pars (struct) - containing additional parameters used tofit
%    residual dynmics
%
%  OUTPUT
%  - impulse_response_groundTruth (cell array) - contianing the ground
%    truth impulse response for different perturbation conditions
%  - impulse_response_predicted (cell array) - contianing the ground
%    truth impulse response for different perturbation conditions
%  - fit_result_all (cell array) - of size (n_areas + 1) x 1 containing the
%    local (first n_areas elements) and the global (last element) residual
%    dynamics used for impulse response prediction.
%  - perturb_dirs_latent (cell array) - information about the applied
%    direction and magnitue of perturbation.
%
%  see sizes of cell arrays in comments below
%  Author : Aniruddh Galgali
%}
% need to set the shuffle option to false as we don't want shuffled
% networks
sim_pars.do_shuffle = false;

% Simulating the "unperturbed" response of the network
[data_unpert] = generate_responses(input_pars, net_pars, sim_pars, obs_pars);


% Simulating the "perturbed" response of the network and subtracting it
% from the "unperturbed" response to compute the ground truth impulse
% response.
n_perturbs = length(perturb_pars);
impulse_response_groundTruth = cell(length(areas),n_perturbs,2); % size :  n_areas x n_perturbs x n_perturb_dirs
perturb_dirs_latent = cell(length(areas),n_perturbs,2,2); % size :  n_areas x n_perturbs x n_perturb_dirs x n_choices
t0_idx = NaN(length(areas),n_perturbs,2,2); % size :  n_areas x n_perturbs x n_perturb_dirs x n_choices
impulse_response_timeaxes = data_unpert.time_rel;

for ia = 1: length(areas)
    
    input_pars_perturb = input_pars;
    
    switch areas{ia}
        case 'pfc'
            
            input_pars_perturb.perturb.dim_mask = [false;false;true;true];
            perturb_dirs = {[1;-1];[1;1]}; %{local_pfc_diff; locl_pfc_sum};
            
        case 'ppc'
            
            input_pars_perturb.perturb.dim_mask = [true;true;false;false];
            perturb_dirs = {[1;-1];[1;1]}; %{local_ppc_diff; locl_ppc_sum};
            
        case 'both'
            
            input_pars_perturb.perturb.dim_mask = [true;true;true;true];
            perturb_dirs = {[1;-1;1;-1];[1;1;1;1]}; %{global_diff; global_sum};
            
    end
    
    for idir = 1: length(perturb_dirs)
        
        input_pars_perturb.perturb.direction = perturb_dirs{idir};
        
        for ipert = 1:n_perturbs
            
            input_pars_perturb.perturb.magnitude = perturb_pars{ipert}.perturb_mag;
            input_pars_perturb.perturb.times = perturb_pars{ipert}.perturb_times;
            
            [data_pert] = generate_responses(input_pars_perturb, net_pars, sim_pars, obs_pars);
            
            [impulse_response_groundTruth{ia,ipert,idir}, perturb_dirs_all, start_idxs_all] = ...
                compute_true_impulse(data_unpert, data_pert, input_pars_perturb.perturb.times, obs_pars, areas,sim_pars.dt);
            
            for icond = 1: size(perturb_dirs_all,2)
                
                perturb_dirs_latent{ia,ipert,idir,icond} = perturb_dirs_all{1,icond};
                t0_idx(ia,ipert,idir,icond) = start_idxs_all;
            end
        end
    end
end



% Predicting the impulse response by first fitting the residual dynamcis and
% then using it to obtain the prediction.
impulse_response_predicted = cell(length(areas),n_perturbs,2); % size :  n_areas x n_perturbs x n_perturb_dirs
[~,~,trial_labels] = unique(data_unpert.task_variable.targ_dir);
fit_result_all = cell(length(areas),1);

for ia = 1: length(areas)
    fit_pars_ = fit_pars;
    switch areas{ia}
        
        case 'ppc'
            
            mask_dim = [true(sim_pars.obs_dim,1); false(sim_pars.obs_dim,1)];
            fit_pars_.max_dim = 2;
            
        case 'pfc'
            
            mask_dim = [false(sim_pars.obs_dim,1); true(sim_pars.obs_dim,1)];
            fit_pars_.max_dim = 2;
            
        case 'both'
            
            mask_dim = [true(sim_pars.obs_dim,1) ; true(sim_pars.obs_dim,1)];
            fit_pars_.dim = 4; % We always fit a full 4 dimensional model to the global dyn
    end
    
    [fit_result] = extract_residualdynamics_twoarearnn(data_unpert.residuals(mask_dim,:,:),...
        data_unpert.time, data_unpert.time_rel, data_unpert.time_iev, trial_labels,...
        result_cv{ia}, fit_pars_);
    
    for iap = 1:size(perturb_dirs_latent,1)
        
        for ipert = 1:size(perturb_dirs_latent,2)
            
            for idir = 1:size(perturb_dirs_latent,3)
                
                for icond = 1:size(perturb_dirs_latent,4)
                    
                    n_dim_fit = size(fit_result.final_model{icond}.A,1);
                    
                    true_perturb_dir = perturb_dirs_latent{iap,ipert,idir,icond};
                    
                    dirs_to_prop = fit_result.U_dyn(:,1:n_dim_fit)'* true_perturb_dir(mask_dim);
                    
                    tRefs.time = data_unpert.time_rel(t0_idx(iap,ipert,idir,icond));
                    tRefs.time_labels = data_unpert.time_iev(t0_idx(iap,ipert,idir,icond));
                    
                    [predicted_response] = computeTimeDependentImpulseResponse(fit_result.final_model{icond}.A, ...
                        fit_result.final_model{icond}.time_rel, dirs_to_prop ,tRefs,'time_labels',fit_result.final_model{icond}.time_iev);
                    
                    propagated_predicted_impulse = [true_perturb_dir(mask_dim) ...
                        fit_result.U_dyn(:,1:n_dim_fit)*squeeze(predicted_response.prop_v)];
                    impulse_response_predicted{iap,ipert,idir}{ia,icond}.full = propagated_predicted_impulse;
                    
                    if(strcmp(areas{ia},'ppc'))
                        
                        impulse_response_predicted{iap,ipert,idir}{ia,icond}.ppc_norm = sqrt(sum(abs(propagated_predicted_impulse).^2,1));
                        impulse_response_predicted{iap,ipert,idir}{ia,icond}.pfc_norm = NaN;
                        
                    elseif(strcmp(areas{ia},'pfc'))
                        
                        impulse_response_predicted{iap,ipert,idir}{ia,icond}.ppc_norm = NaN;
                        impulse_response_predicted{iap,ipert,idir}{ia,icond}.pfc_norm = sqrt(sum(abs(propagated_predicted_impulse).^2,1));
                        
                    elseif(strcmp(areas{ia},'both'))
                        
                        impulse_response_predicted{iap,ipert,idir}{ia,icond}.ppc_norm = sqrt(sum(abs(propagated_predicted_impulse...
                            (1:ceil(size(obs_pars.C_obs,1)/2),:)).^2,1));
                        impulse_response_predicted{iap,ipert,idir}{ia,icond}.pfc_norm = sqrt(sum(abs(propagated_predicted_impulse...
                            (ceil(size(obs_pars.C_obs,1)/2) + 1 :end,:)).^2,1));
                        
                    end
                    
                    impulse_response_predicted{iap,ipert,idir}{ia,icond}.time = [
                        predicted_response.time predicted_response.time(end) + ...
                        (predicted_response.time(end) - predicted_response.time(end-1))];
                    
                end
                
            end
            
        end
        
    end
    
    fit_result_all{ia} = fit_result;
    
end

varargout{1} = impulse_response_timeaxes;
end


function [impulse_response_groundTruth, perturb_dirs_all, t0_idxs] = ...
    compute_true_impulse(data_unpert, data_pert, times_perturb, obs_pars, areas, dt)

%{ This function computes the ground truth impulse response.
%
%  INPUT
%  - data_unpert (struct) - "unperturbed" simulated responses
%  - data_pert (struct) - "perturbed" simulated responses
%  - times_perturb (2 x 1 array- [t_start t_end]) - indicating the stard
%    and end times of the applied perturbation
%  - obs_pars   (struct) -  containing parameters of the linear gaussian
%               osbervation model
%  - areas (cell array) - size of (n_areas + 1) x 1, containing the string
%    labels that define each fit condition. The first n_areas correspond to 
%    each of the individual areas of the model, and the last element
%    corresponds to the string label associated with global dynamics.
%  - dt (scalar) - the time-step of the simulation

%
%  OUTPUT
%  - impulse_response_groundTruth (cell array) - contianing the ground
%    truth impulse response for a specific perturbation condition, measured
%    in each area and for all task conditions.
%  - perturb_dirs_all(cell array) - information about the applied
%    direction and magnitude of perturbation for each task condition.   
%  - t0_idxs - time indices within the trial corresponding to the perturbation.    
%   
%}
pert_times = times_perturb * dt; % This is the step size of the network

[~,t0_idxs] = min(abs(data_unpert.time_rel - median(pert_times)));

[~,~,trial_labels_unpert] = unique(data_unpert.task_variable.targ_dir);
unique_trial_labels = unique(trial_labels_unpert);
nconds = length(unique_trial_labels);
impulse_response_groundTruth = cell(length(areas),nconds);
perturb_dirs_all = cell(1,nconds);

for icond = 1:nconds
        

        cond_trial_idxs = trial_labels_unpert == unique_trial_labels(icond);
        unperturb_mean = mean(data_unpert.response(:,:,cond_trial_idxs),3);
        
        residuals_perturb = data_pert.response(:, :, trial_labels_unpert == unique_trial_labels(icond)) - ...
            unperturb_mean;
        
        impulse_response = squeeze(mean(residuals_perturb,3));
        
        for ia = 1: length(areas)
            
            switch areas{ia}
                
                case 'ppc'
                    
                    mask_dim = [true(size(obs_pars.C_obs,1)/2,1); ...
                        false(size(obs_pars.C_obs,1)/2,1) ];
                    
                case 'pfc'
                    
                    mask_dim = [false(size(obs_pars.C_obs,1)/2,1); ...
                        true(size(obs_pars.C_obs,1)/2,1) ];
                    
                case 'both'
                    
                    mask_dim = [true(size(obs_pars.C_obs,1),1)];
                    
            end
            
            impulse_response_groundTruth{ia,icond}.full = impulse_response(mask_dim,:);
            impulse_response_groundTruth{ia,icond}.norm = sqrt(sum(abs(impulse_response(mask_dim,:)).^2,1));
            impulse_response_groundTruth{ia,icond}.time = data_unpert.time;
            
        
        end
        
        perturb_dirs_all{icond} =  impulse_response(:,t0_idxs);
    

end

end


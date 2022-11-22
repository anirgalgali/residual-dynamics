%{ This script computes the ground truth and predicted impulse response for a specified
%  two-area model. Typically, you want o save all the analyses of the
%  2-area model (data/fits) before running this script.
%}
%% clear workspace

clearvars -except DIRS
clc
close all

%% Loading stored results

file_path = './data/twoarearnn/';
file_name = 'nofeedbackmodels.mat'; % choose either 'withfeedback' or 'nofeedback' here
do_save_vars = false;

% Specifying which models and choice goes into the final plot
% This specifies all the relevant model configurations for both types of
% networks (i.e with or without feedback)
model_pars = table([0.3600;0.2000],...
                   [0.0800;0.2400],...
                   [0;0.2400],...
             'VariableNames',{'self' 'fwd' 'fbk'});

[~,~,~,~,model_result] = extract_model_analyses(file_path,file_name,model_pars);

%% Computing the impulse response

% Declaring the fit parameres for the residual dynamics
fit_pars.hankel_criterion = 'sem';
fit_pars.hankel_type = 'combined';
fit_pars.lag = 3;
 
% Defining perturbation paraemters

% times in trial at which to aply perturbation
perturb_time_all = {[2701:2702];[4051:4052];[5401:5402];[6751:6752];[8101:8102];[9451:9452]};

% magnitude of perturbation
perturbation_magnitude = {0.001};
all_perturb_combs = allcomb(perturbation_magnitude,perturb_time_all);
perturb_pars = cell(length(perturb_time_all),1);
for kk = 1: size(all_perturb_combs,1)
    
    perturb_pars{kk}.perturb_mag = all_perturb_combs{kk,1};
    perturb_pars{kk}.perturb_times = all_perturb_combs{kk,2} ;% should be a perfect multiple of 45

end

% Computing impulse responses
if(length(model_result) == 1)
   
    fprintf('Only one model to be analyzed\n')
    
    [analysis.impulse_response_groundTruth,analysis.impulse_response_predicted,...
        analysis.fit_result, analysis.measured_perturbs] = ...
        compute_twoarea_impulse_responses(model_result{1}.input_pars,...
    model_result{1}.sim_pars, model_result{1}.net_pars, model_result{1}.obs_pars,...
    model_result{1}.areas, model_result{1}.result_cv, perturb_pars, fit_pars);
 
    analysis.perturb_pars = perturb_pars;
    analysis.net_pars = model_result{1}.net_pars;
    analysis.input_pars = model_result{1}.input_pars;
    analysis.sim_pars = model_result{1}.sim_pars;
    analysis.obs_pars = model_result{1}.obs_pars;
    analysis.fit_pars = fit_pars;
    analysis.areas = model_result{1}.areas;
    analysis.root_file_name = file_name;
    
    % Saving results to disk 
    if(do_save_vars)
        
        save_file_name = sprintf('%s%.2f%s%.2f%s%s%.3f%s',....
            'ImpulseResponse_self=',model_result{1}.net_pars.self_pfc,...
            '_cross=',model_result{1}.net_pars.fwd_ppc_pfc,strrep(file_name,...
            'models.mat',[]),'_multiplePerturbsmag=',perturbation_magnitude{1},'.mat');
        
        save(fullfile(file_path,save_file_name),'analysis','-v7.3')
    end
    
else
    
    fprintf('Multiple models to be analyzed\n')
    analysis = cell(length(model_result),1);
    for imdl = 1: length(model_result)
        
        [analysis{imdl}.impulse_response_groundTruth,analysis{imdl}.impulse_response_predicted,...
            analysis{imdl}.fit_result, analysis{imdl}.measured_perturbs] = ...
            compute_twoarea_impulse_responses(model_result{imdl}.input_pars,...
            model_result{imdl}.sim_pars, model_result{imdl}.net_pars, model_result{imdl}.obs_pars,...
            model_result{imdl}.areas, model_result{imdl}.result_cv, perturb_pars, fit_pars);
        
        analysis{imdl}.perturb_pars = perturb_pars;
        analysis{imdl}.model_result = model_result{imdl};
        analysis{imdl}.fit_pars = fit_pars;
        analysis.root_file_name = file_name;
  
    end
 
end



%% Plotting the results

close all

n_look_ahead = 3;%  Number of future time-points (relative to perturbation) to plot
all_fit_conditions = model_result{1}.areas;
local_areas = all_fit_conditions(~ismember(all_fit_conditions, 'both'));

%dimensions along which we apply perturbation
perturb_dim_labels = {'choice','time'};
num_perturb_dirs = length(perturb_dim_labels);

% colors for ground truth
cols_mat_gt = [0.0430    0.3555    0.1016;  % choice...
    0.4609    0.1797    0.6836]; % time
assert(size(cols_mat_gt,1) == num_perturb_dirs,'number of colors should be equal to number of perturbation directions');

% colors for prediction
cols_mat_pred = [0.5 0.5 0.5;
    0.5 0.5 0.5;
    0  0  0];
assert(size(cols_mat_pred,1) == length(all_fit_conditions),'number of colors should be equal to number of fit conditions');

cidx = 1; % task condition for which we want to plot results (cidx = 1 implies choice-1 condition
figure;  

for ia1 = 1: length(local_areas) % indexes the perturbed area
    for ia3 = 1: length(local_areas) % indexes the measured area
        
        for idir= 1:num_perturb_dirs % indexes the perturbation direction
            
            ah = subplot(length(local_areas),length(local_areas) * num_perturb_dirs,...
                (ia1 - 1)*(length(local_areas) * num_perturb_dirs) + ((ia3 - 1)*num_perturb_dirs) + idir);
            hold(ah);set(ah,'plotbox',[1 1 1]);
            
            if(ia3 == 1)
                set(ah,'ylim',[0 0.025]);
            else
                set(ah,'ylim',[0 0.00725]);
            end
            
            imp_gts = NaN(length(analysis.perturb_pars),n_look_ahead);
            imp_preds = NaN(length(analysis.perturb_pars),n_look_ahead,length(all_fit_conditions));
            
            for kk = 1: length(analysis.perturb_pars)
                
                t0 = analysis.impulse_response_predicted{ia1,kk,idir}{1,cidx}.time(1);
                
                if(kk < length(analysis.perturb_pars))
                    
                    t0_next = analysis.impulse_response_predicted{ia1,kk+1,idir}{1,cidx}.time(1);
                    time_idxs_gt = find(analysis.impulse_response_groundTruth{ia1,kk,idir}{ia3,cidx}.time >= t0 & ...
                        analysis.impulse_response_groundTruth{ia1,kk,idir}{ia3,cidx}.time < t0_next);
                    
                    
                else
                    time_idxs_gt = find(analysis.impulse_response_groundTruth{ia1,kk,idir}{ia3,cidx}.time >= t0);
                end
                
                time_axes_gt = analysis.impulse_response_groundTruth{ia1,kk,idir}{ia3,cidx}.time(time_idxs_gt);
                imp_gt = analysis.impulse_response_groundTruth{ia1,kk,idir}{ia3,cidx}.norm(time_idxs_gt);
                
                plot(ah,time_axes_gt(1:n_look_ahead),imp_gt(1:n_look_ahead),...
                    '-','Color',cols_mat_gt(idir,:));
                plot(ah,time_axes_gt(1:n_look_ahead),imp_gt(1:n_look_ahead),...
                    'o','markeredgecolor',cols_mat_gt(idir,:),'markerfacecolor',[1 1 1],'markersize',4);

                
                imp_gts(kk,:) = imp_gt(1:n_look_ahead);
                
                
                for ia2 = 1: length(all_fit_conditions) % indexes the fit type
                    
                    if(kk < length(analysis.perturb_pars))
                        
                        time_idxs_pred = find(analysis.impulse_response_predicted{ia1,kk,idir}{ia2,cidx}.time >= t0 & ...
                            analysis.impulse_response_predicted{ia1,kk,idir}{ia2,cidx}.time < t0_next);
                    else
                        
                        time_idxs_pred = find(analysis.impulse_response_predicted{ia1,kk,idir}{ia2,cidx}.time >= t0);
                        
                    end
                    time_axes_pred = analysis.impulse_response_predicted{ia1,kk,idir}{ia2,cidx}.time(time_idxs_pred);
                    
                    if(~isnan(analysis.impulse_response_predicted{ia1,kk,idir}{ia2,cidx}.([all_fit_conditions{ia3} '_norm'])))
                        imp_p = analysis.impulse_response_predicted{ia1,kk,idir}{ia2,cidx}.([all_fit_conditions{ia3} '_norm'])(time_idxs_pred);
                        plot(ah,time_axes_pred(1:n_look_ahead),imp_p(1:n_look_ahead),'-','color',cols_mat_pred(ia2,:));
                    end

                    imp_preds(kk,:,ia2) = imp_p(1:n_look_ahead);
                        
                    if(ia1 == 2)
                       xlabel('time (s)') 
                    end
                    
                    if(ia3 == 1)
                        ylabel('norm (a.u)')
                    end
                    
                    title(ah,sprintf('%s%s%s%s%s%s','perturb-',local_areas{ia1},', measure-',local_areas{ia3}, ', ', perturb_dim_labels{idir}),...
                        'FontSize', 8)
                end
            end
         
            
        end

    end
    
end
pos = get(gcf,'Position');
scale_fac = 2;
set(gcf,'Position',[pos(1) pos(2) scale_fac*pos(3) scale_fac*pos(4)])
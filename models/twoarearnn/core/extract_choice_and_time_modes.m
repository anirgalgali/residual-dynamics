function [data] = extract_choice_and_time_modes(data)
%{ This function extracts activity time-courses along the choice and time 
%  modes of the network.
%  INPUT
%  - data : struct containing the relevant responses of the model
%  OUTPUT
%  - data : same as above but wiith additional 'task_modes' field.
%  Author: Aniruddh Galgali (Nov 2020)
%}
areas = {'ppc';'pfc'};

u_choice_ppc = normc([1;-1;0;0]);
u_choice_pfc  = normc([0;0;1;-1]);

u_time_ppc = normc([1;1;0;0]);
u_time_pfc = normc([0;0;1;1]);


if(strcmp(areas{1},'ppc'))
    u_choice = {u_choice_ppc, u_choice_pfc};
    u_time = {u_time_ppc, u_time_pfc};
elseif(strcmp(areas{1},'pfc'))
    u_choice = {u_choice_pfc, u_choice_ppc};
    u_time = {u_time_pfc, u_time_ppc};
else
    error('invalid area names') 
end
    

trial_labels = data.task_variable.targ_dir;
unique_trial_labels = unique(trial_labels);


for iarea = 1: length(areas)
     for icond = 1:length(unique_trial_labels)

        data.task_modes.(areas{iarea}).choice.proj(icond,:) = u_choice{iarea}'*squeeze(mean(data.rates(:,:,trial_labels == unique_trial_labels(icond)),3));
        data.task_modes.(areas{iarea}).time.proj(icond,:) = u_time{iarea}'*squeeze(mean(data.rates(:,:,trial_labels == unique_trial_labels(icond)),3));
        data.task_modes.(areas{iarea}).choice.direction = u_choice{iarea};
        data.task_modes.(areas{iarea}).time.direction = u_time{iarea};
     
     end

end
       
end
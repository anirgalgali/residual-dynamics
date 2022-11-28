function [target_params] = get_saccadetarget_configurations(expt, target_discrepancy_threshold)
%{ Extracts the locations of the saccade targets and the assignment of the 
%  task configuration for a single experiment. 
% Input
% - expt(struct) - contians the data and associated variables for a single
%   recording session 
% - target_discrepancy_threshold (scalar) -a threshold used to determine
%   how to assign target locations to configurations, when tagret locations
%   vary randomly from session to session.
% Output
% - target_params (struct) - contains information about the target
%   locations and configuration
% Author : Aniruddh Galgali (Jan 2022)
%}
evnames = fieldnames(expt.align);
target_params = [];

if(isfield(expt.targInfo,'screenOffsetX') || isfield(expt.targInfo,'screenOffsetY'))    
    truePosX_FP = expt.targInfo.screenOffsetX;
    truePosY_FP = expt.targInfo.screenOffsetY;
else
    truePosX_FP = 0;
    truePosY_FP = 0;
end

truePosX_T1 = expt.targInfo.T1x;
truePosX_T2 = expt.targInfo.T2x;
truePosY_T1 = expt.targInfo.T1y;
truePosY_T2 = expt.targInfo.T2y;

relPosX_T1 = truePosX_T1 - truePosX_FP;
relPosX_T2 = truePosX_T2 - truePosX_FP;
relPosY_T1 = truePosY_T1 - truePosY_FP;
relPosY_T2 = truePosY_T2 - truePosY_FP;

T1_angle = cart2pol(relPosX_T1,relPosY_T1);
T2_angle = cart2pol(relPosX_T2,relPosY_T2);


if(T1_angle < 0 & T1_angle > -pi)
    
    T1_angle = wrapTo2Pi(T1_angle);
    
end

if(T2_angle < 0 & T2_angle > -pi)
    
    T2_angle = wrapTo2Pi(T2_angle);
    
end

[allTargetAngles_rad,ind_rad]  = sort([T1_angle;T2_angle]);
[allTargetAngles_deg,ind_deg] = sort([round(rad2deg(T1_angle));round(rad2deg(T2_angle))]);


target_params.truePos_FP = [truePosX_FP truePosY_FP];
target_params.truePos_T1 = [truePosX_T1 truePosY_T1];
target_params.truePos_T2 = [truePosX_T2 truePosY_T2];

target_params.relPos_T1 = [relPosX_T1 relPosY_T1];
target_params.relPos_T2 = [relPosX_T2 relPosY_T2];
target_params.targetAngles = allTargetAngles_deg;
labels = {'T1';'T2'};

target_params.labels = labels;

% Purging NaNs
if(isfield(expt.align.(evnames{1}).data.task_variable,'choicesactheta'))
    
    idx_valid_theta = ~isnan(expt.align.(evnames{1}).data.task_variable.choicesactheta);
    sac_angles_rad =  mod(deg2rad(expt.align.(evnames{1}).data.task_variable.choicesactheta(idx_valid_theta)),2*pi);
    ichoice_valid = expt.align.(evnames{1}).data.task_variable.i_choicecd(idx_valid_theta);
    targdir_valid = expt.align.(evnames{1}).data.task_variable.targ_dir(idx_valid_theta);
    correct_valid = expt.align.(evnames{1}).data.task_variable.correct(idx_valid_theta);
    dotscoh_valid = expt.align.(evnames{1}).data.task_variable.dots_coh(idx_valid_theta) ~= 0;
    
    targ_idx_theta = NaN(length(allTargetAngles_rad),1);
    targ_var_theta = NaN(length(allTargetAngles_rad),1);
    
    for ii = 1:length(allTargetAngles_rad)
        
        if(allTargetAngles_rad(ii) == 0)
            
            idx_theta = sac_angles_rad <= allTargetAngles_rad(ii) + target_discrepancy_threshold | sac_angles_rad >= 2*pi - target_discrepancy_threshold;
            
        else
            
            idx_theta = sac_angles_rad <= allTargetAngles_rad(ii) + target_discrepancy_threshold & sac_angles_rad >= allTargetAngles_rad(ii) - target_discrepancy_threshold;
            
        end
        
        idx_theta_correct = idx_theta & correct_valid & dotscoh_valid ;
        targ_idx_theta(ii) = unique(ichoice_valid(idx_theta_correct));
        targ_var_theta(ii) = unique(targdir_valid(idx_theta_correct));
        
%         fprintf('Angle = %d, targ_idx = %d \n',allTargetAngles_deg(ii), targ_idx_theta(ii));
   
    end
    
    if(allTargetAngles_deg(1) == 0 & allTargetAngles_deg(2) == 180)
        
        targetConfig = 2;
        pref_idx = targ_idx_theta(allTargetAngles_deg == 0);
        anti_idx = targ_idx_theta(allTargetAngles_deg == 180);
        
    elseif(allTargetAngles_deg(1) >= 0 & allTargetAngles_deg(1) < 90 & allTargetAngles_deg(2) <= 270 & allTargetAngles_deg(2) > 180)
        
        targetConfig = 3;
        pref_idx = targ_idx_theta(allTargetAngles_deg >= 0 & allTargetAngles_deg < 90);
        anti_idx = targ_idx_theta(allTargetAngles_deg <= 270 & allTargetAngles_deg > 180);
    
    
    elseif(allTargetAngles_deg(1) > 0 & allTargetAngles_deg(1) <= 90 & allTargetAngles_deg(2) <= 360 & allTargetAngles_deg(2) > 270)
    
        targetConfig = 4;
        pref_idx = targ_idx_theta(allTargetAngles_deg > 0 & allTargetAngles_deg <= 90);
        anti_idx = targ_idx_theta(allTargetAngles_deg <= 360 & allTargetAngles_deg > 270);
        
        
    elseif(allTargetAngles_deg(1) > 90 & allTargetAngles_deg(1) <= 180 & allTargetAngles_deg(2) <= 360 & allTargetAngles_deg(2) > 270)
        
        targetConfig = 1;
        pref_idx = targ_idx_theta(allTargetAngles_deg <= 360 & allTargetAngles_deg > 270);
        anti_idx = targ_idx_theta(allTargetAngles_deg > 90 & allTargetAngles_deg <= 180);
        
        
    elseif(allTargetAngles_deg(1) > 0 & allTargetAngles_deg(1) <= 90 & allTargetAngles_deg(2) <= 180 & allTargetAngles_deg(2) > 90)
        
        targetConfig = 2;
        pref_idx = targ_idx_theta(allTargetAngles_deg > 0 & allTargetAngles_deg <= 90);
        anti_idx = targ_idx_theta(allTargetAngles_deg > 90 & allTargetAngles_deg <= 180);
        
    elseif(allTargetAngles_deg(1) >= 180 & allTargetAngles_deg(1) < 270 & allTargetAngles_deg(2) <= 360 & allTargetAngles_deg(2) > 270)
        
        targetConfig = 2;
        pref_idx = targ_idx_theta(allTargetAngles_deg > 270 & allTargetAngles_deg <= 360);
        anti_idx = targ_idx_theta(allTargetAngles_deg >= 180 & allTargetAngles_deg < 270);

    end
    
else
    
    warning('choice sac theta info not available')
    target_params.targdir_vals = NaN;
    target_params.targdir_idx  = NaN;
    target_params.pref_idx = NaN;
    target_params.anti_idx = NaN;
    target_params.config = NaN;
    return;
    
end

target_params.targdir_vals = targ_var_theta;
target_params.targdir_idx = targ_idx_theta;

target_params.config = targetConfig;
target_params.pref_idx = pref_idx;
target_params.anti_idx = anti_idx;
end

function extractBehaviour(animal,doSegment)

%% Loading datasets

pexp_preP = ProjectLoad_v2('array_reppasdotsTask_datasegmented_binsize=45ms');
pexp_data = ProjectLoad('array_TDR_dotsTask_data');
pexp_beh = ProjectCreate(strcat('array_dotsTask_behv2','_',animal));
%% Extracting data indices relevant to a particular animal

idx_animal = false(1,length(pexp_data.idlist));
for iexp = 1:length(pexp_data.idlist)
    
    idx_animal(iexp) = any(strfind(pexp_data.idlist{iexp},[animal]));
    
end
idlist_animal = pexp_data.idlist(idx_animal);

%% Run Behaviour Analysis

nExps = length(idlist_animal);
analysis_out = [];

for iExp = 1:nExps
    
    analysis = AnalysisLoad(pexp_data.name,idlist_animal{iExp},'fast');
    [analysis_preprocess] = AnalysisLoad_v2(pexp_preP.name,idlist_animal{iExp});
    [targetParams] = extractTargetandConfigVariables(analysis,'stimoncd');

    analysis_out.experiment = analysis.experiment;
    analysis_out.thisFile = analysis.thisFile;
    analysis_out.targInfo = targetParams;
    analysis_out.segments.trialNum = analysis_preprocess.session_segmentation_info.trialNumSegment;
    analysis_out.segments.idxs = analysis_preprocess.session_segmentation_info.idxsSegment;
    analysis_out.segments.valid = analysis_preprocess.session_segmentation_info.segment_mask;
    
    if(~isnan(targetParams.pref_idx))
        [perChoiceToPref,dotCoh] = computePsychometricData(analysis,analysis_preprocess,targetParams,'stimoncd',doSegment);
%         segmentFlag_all(s_count : s_count + length(perChoiceToPref) - 1) = analysis_preprocess.segmentsToInclude;
%         exp_cor_all(s_count : s_count + length(perChoiceToPref) - 1) = perChoiceToPref ;
%         dot_coh_all(s_count : s_count + length(perChoiceToPref) - 1) = dotCoh ;
%         s_count = s_count + length(perChoiceToPref);
        analysis_out.behaviour.percChoiceToPref = perChoiceToPref; 
        analysis_out.behaviour.signed_coh = dotCoh;
        
    else
        analysis_out.behaviour.percChoiceToPref = {NaN};
        analysis_out.behaviour.signed_coh = {NaN};
        
    end
    
    AnalysisSave(pexp_beh.name,idlist_animal{iExp},analysis_out)
    fprintf('Done with Exp - %d\n',iExp);
end


end


function [targetParams] = extractTargetandConfigVariables(analysis,evname)

targetParams = [];
discrepancyThreshold = 0.4;
if(isfield(analysis.targInfo,'screenOffsetX') || isfield(analysis.targInfo,'screenOffsetY'))
    
    truePosX_FP = analysis.targInfo.screenOffsetX;
    truePosY_FP = analysis.targInfo.screenOffsetY;
    
else
    
    truePosX_FP = 0;
    truePosY_FP = 0;
    
end

truePosX_T1 = analysis.targInfo.T1x;
truePosX_T2 = analysis.targInfo.T2x;
truePosY_T1 = analysis.targInfo.T1y;
truePosY_T2 = analysis.targInfo.T2y;

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


targetParams.truePos_FP = [truePosX_FP truePosY_FP];
targetParams.truePos_T1 = [truePosX_T1 truePosY_T1];
targetParams.truePos_T2 = [truePosX_T2 truePosY_T2];

targetParams.relPos_T1 = [relPosX_T1 relPosY_T1];
targetParams.relPos_T2 = [relPosX_T2 relPosY_T2];
targetParams.targetAngles = allTargetAngles_deg;
labels = {'T1';'T2'};

targetParams.labels = labels;

% Purging NaNs
if(isfield(analysis.align.(evname).data.task_variable,'choicesactheta'))
    
    idx_valid_theta = ~isnan(analysis.align.(evname).data.task_variable.choicesactheta);
    sac_angles_rad =  mod(deg2rad(analysis.align.(evname).data.task_variable.choicesactheta(idx_valid_theta)),2*pi);
    ichoice_valid = analysis.align.(evname).data.task_variable.i_choicecd(idx_valid_theta);
    targdir_valid = analysis.align.(evname).data.task_variable.targ_dir(idx_valid_theta);
    correct_valid = analysis.align.(evname).data.task_variable.correct(idx_valid_theta);
    dotscoh_valid = analysis.align.(evname).data.task_variable.dots_coh(idx_valid_theta) ~= 0;
    for ii = 1:length(allTargetAngles_rad)
        
        if(allTargetAngles_rad(ii) == 0)
            
            idx_theta = sac_angles_rad <= allTargetAngles_rad(ii) + discrepancyThreshold | sac_angles_rad >= 2*pi - discrepancyThreshold;
            
        else
            
            idx_theta = sac_angles_rad <= allTargetAngles_rad(ii) + discrepancyThreshold & sac_angles_rad >= allTargetAngles_rad(ii) - discrepancyThreshold;
            
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
    targetParams.targdir_vals = NaN;
    targetParams.targdir_idx  = NaN;
    targetParams.pref_idx = NaN;
    targetParams.anti_idx = NaN;
    targetParams.config = NaN;
    
    
    
    
    
    return;
    
end

targetParams.targdir_vals = targ_var_theta;
targetParams.targdir_idx = targ_idx_theta;

targetParams.config = targetConfig;
targetParams.pref_idx = pref_idx;
targetParams.anti_idx = anti_idx;

end


function [expCor,dotCoh] = computePsychometricData(analysis,analysis_preprocess,targetPars,evname,doSegment)

pref_targ_idx = targetPars.targdir_idx == targetPars.pref_idx;

if(doSegment)
    
    idxsAllSegments = analysis_preprocess.session_segmentation_info.idxsSegment;
    expCor = cell(length(idxsAllSegments),1);
    
    for ii = 1:length(idxsAllSegments)
        
        getpars.trial_get = idxsAllSegments{ii};
        data_segment = tdrGetTrials(analysis.align.(evname).data,getpars);
        
        [dotCoh{ii},expCor{ii}] = psychfunction(data_segment.task_variable.sign_coh,data_segment.task_index.targ_dir == targetPars.targdir_idx(pref_targ_idx));
        
        if(targetPars.targdir_vals(pref_targ_idx) == -1)
           
            dotCoh{ii} = -1.*dotCoh{ii};
            
        end
    
    end
    
else
    
    
    [dotCoh,expCor] = psychfunction(analysis.align.(evname).data.task_variable.sign_coh,...
        analysis.align.(evname).data.task_index.targ_dir == targetPars.targdir_idx(pref_targ_idx));

end

  
end

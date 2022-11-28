function [analysis] = extract_alignedtaskplanes(data, trial_labels, jPC_pars, choice_time_pars,alignment_labels)
%{ This function extracts all the task-planes using the aligned neural responses
%  Input
%  - data (struct) - contains the single trial data and associated inforrmation
%  - trial_labels (array) - of size 1 x n_trials,  indicating the task condition
%    associated with each trial. n_trials -> total number of trials in data.
%  - jPC_pars (struct) - contains parameters used for computing jPCA
%  - choice_time_pars (struct) - contains parameters used for computing
%    choice/time planes
%  - alignment_labels(cell array) - of size n_Align x 1. Each element of
%    this cell array contains a string label indicating the different 
%    task "epochs" in the data.
%  Output
%  - analysis(struct) - contains all the information about the JPCs,
%    choice/time planes and the respective projection of the
%    condition-averaged trajectories along these planes
%
%  Author: Aniruddh Galgali (Jan 2021)
%}
time_rel = data.time_rel;
time_labels = data.time_iev;
unique_time_labels = unique(time_labels);
unique_trial_labels = unique(trial_labels);
nalign = length(unique_time_labels);
nconds = length(unique_trial_labels); 
assert(length(alignment_labels) == nalign,'invalid number of alignments');    
% Extract jPCA subspaces and projections

analyze_time_mask = cell(nalign,1);    
for ialign = 1:nalign
    time_rel_ev = time_rel(time_labels == unique_time_labels(ialign));
    analyze_time_mask{ialign} = time_rel_ev >= jPC_pars.time_win_analyze{ialign}(1) & ...
        time_rel_ev <= jPC_pars.time_win_analyze{ialign}(end);
end
    
if(jPC_pars.do_pre_smooth)
    filt_pars = [];
    filt_pars.filt_width = 4*(time_rel(2) - time_rel(1));
    
    [jPCSummary] = computejPCA(data.response, time_rel, trial_labels, time_labels,...
        analyze_time_mask, jPC_pars.n_jpc_planes, filt_pars);
    
else  
    [jPCSummary] = computejPCA(data.response, time_rel, trial_labels, time_labels,...
        analyze_time_mask, jPC_pars.n_jpc_planes);
end

for ialign = 1:nalign
    for ipair = 1: size(jPCSummary(ialign).jPCs_fullD, 2) / 2
        
        vi1 = 1+2*(ipair-1);
        vi2 = 2*ipair;
        analysis.(alignment_labels{ialign}).directions{ipair} = jPCSummary(ialign).jPCs_fullD(:,[vi1 vi2]);
        
        for icond = 1:nconds
            analysis.(alignment_labels{ialign}).projections{ipair}(:,:,icond) = jPCSummary(ialign).Projection(icond).projAllTimes(:,[vi1 vi2])';
        end

        analysis.(alignment_labels{ialign}).direction_labels{ipair} = sprintf('%s%d%d','jPC',vi1,vi2);
        analysis.(alignment_labels{ialign}).variance_explained(ipair) = jPCSummary(ialign).varCaptEachPlane(ipair);
        analysis.(alignment_labels{ialign}).rot_evals{ipair} = jPCSummary(ialign).evals([vi1 vi2]); 
    end
    analysis.(alignment_labels{ialign}).time = time_rel(time_labels  == ialign);
end

% Extract choice-time subspace    
[D_cond] = sort_trials_by_condition(data, choice_time_pars.condition_vars);
for ialign = 1:nalign
    
    Xavg = [];    
    for icond = 1: nconds
        
        Xavg(:,:,unique_trial_labels(icond)) = ...
            mean(D_cond{unique_trial_labels(icond)}.response(:, ...
            time_labels == unique_time_labels(ialign),:),3);
    end
    
    CD_t = squeeze(Xavg(:,:, choice_time_pars.pref_idx) - Xavg(:,:,choice_time_pars.anti_idx));
    CD_i = 0.5*(squeeze(Xavg(:,:, choice_time_pars.pref_idx) + Xavg(:, :, choice_time_pars.anti_idx)));
    
    [coeff_t,~,~,~,~,~] = pca(CD_t(:,:)');
    [coeff_i,~,~,~,~,~] = pca(CD_i(:,:)');
    
    [expvar_t, Xavg_msub] = compute_variance_explained_by_projection(Xavg, coeff_t(:,1:2), [], []);
    expvar_i = compute_variance_explained_by_projection(Xavg, coeff_i(:,1:2), [], []);
    
    Xavg_proj_t = reshape((Xavg_msub(:,:)' * coeff_t(:,1:2))',[2 size(Xavg,2) size(Xavg,3)]);
    Xavg_proj_i = reshape((Xavg_msub(:,:)' * coeff_i(:,1:2))',[2 size(Xavg,2) size(Xavg,3)]);
    Xavg_proj_all = {Xavg_proj_t, Xavg_proj_i};
    
    dirs_all = {coeff_t(:,1:2),coeff_i(:,1:2)};
    dirs_label_all = {'choice';'time'};
    explained_var_all = {expvar_t, expvar_i};
    
    for ipair = 1: length(Xavg_proj_all)
    
        analysis.(alignment_labels{ialign}).projections = cat(2,analysis.(alignment_labels{ialign}).projections,Xavg_proj_all{ipair});
        analysis.(alignment_labels{ialign}).directions = cat(2,analysis.(alignment_labels{ialign}).directions, dirs_all{ipair});
        analysis.(alignment_labels{ialign}).direction_labels = cat(2,analysis.(alignment_labels{ialign}).direction_labels, dirs_label_all{ipair});
        analysis.(alignment_labels{ialign}).variance_explained = cat(2,analysis.(alignment_labels{ialign}).variance_explained,explained_var_all{ipair}(end));
        analysis.(alignment_labels{ialign}).rot_evals= cat(2,analysis.(alignment_labels{ialign}).rot_evals, {NaN});
    end

    analysis.(alignment_labels{ialign}).Xavg = Xavg_msub;
end

end
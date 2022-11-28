function [analysis] = compute_alignment_and_residuals(data_dir, animal, pars)

[valid_session_list] = get_session_list(data_dir, animal);
data = {};
dsp_info = {};
expt_name = {};
for isession = 1: length(valid_session_list)
    
    expt = AnalysisLoad_v2(data_dir.name, valid_session_list{isession},'fast');
    data = cat(2,data,expt.data);
    dsp_info = cat(2,dsp_info, expt.dspInfo);
    expt_name = cat(2,expt_name, repmat({expt.thisFile},[1 length(expt.data)]));
end

num_trials_per_session = cell2mat(cellfun(@(x) size(x.response,3),data,'uni',false));
zero_trial_exps = num_trials_per_session == 0;
data = data(~zero_trial_exps);
dsp_info = dsp_info(~zero_trial_exps);
fprintf('Analyzing a total of %d experiments in %s\n',length(data),animal);


if(~isempty(pars.config_vars))
    config_variable = cell2mat(cellfun(@(x) unique(x.task_variable.(pars.config_vars)),data,'uni',false));
else
    config_variable = ones(length(data),1); % trivial, all sessions are lumped into a single "configuration"
end

[unique_configs] = unique(config_variable);
num_unique_configs = length(unique_configs);
fprintf('Number of task configurations to consider in %s = %d\n',animal, num_unique_configs);

aligned_single_trials = cell(num_unique_configs,1);
aligned_residuals = cell(num_unique_configs,1);
expt_names_in_config = cell(num_unique_configs,1);
align_projs = cell(num_unique_configs,1);
task_conds_proj = cell(num_unique_configs,1);
align_stats = cell(num_unique_configs,1);
excluded_trialidxs_from_residuals = cell(num_unique_configs,1);
for icfg = 1: num_unique_configs
    
    idx_cfg = config_variable == unique_configs(icfg);
    expt_names_in_config{unique_configs(icfg)} = expt_name(idx_cfg);
    switch pars.alignment.type 
        
        case 'full'
            
            if(icfg == 1)
                fprintf('Computing the alignment\n')
            end
            
            [aligned_single_trials{unique_configs(icfg)}, align_projs{unique_configs(icfg)},...
            task_conds_proj{unique_configs(icfg)}, align_stats{unique_configs(icfg)}] = ...
            compute_session_alignment(data(idx_cfg),pars.align.condition_vars,pars);
        
        case 'precomputed'
            
            if(icfg == 1)
                fprintf('Using a precomputed alignment\n')
            end
            assert(isfield(pars,'precomputed_alignment'), 'need to provide the alignment explicitly')
            [aligned_single_trials{unique_configs(icfg)}] = ...
                compute_prealigned(data(idx_cfg),pars.precomputed_alignment.U_orth{unique_configs(icfg)},...
                pars.precomputed_alignment.X_mu{unique_configs(icfg)});

        case 'none'
            
            if(icfg == 1)
                fprintf('Foregoing session alignment\n')
            end
            aligned_single_trials{unique_configs(icfg)} = data(idx_cfg);
    
    end
    fprintf('Number of experiments for task configuration %d = %d\n', unique_configs(icfg), sum(idx_cfg));
        
    [aligned_residuals{unique_configs(icfg)}, excluded_trialidxs_from_residuals{unique_configs(icfg)}] = ...
        compute_residual_trajectories(aligned_single_trials{unique_configs(icfg)},...
        pars.residual.condition_vars,pars.residual.exclude_trials);
    
    [aligned_single_trials{unique_configs(icfg)}] = filter_excluded_trials...
        (aligned_single_trials{unique_configs(icfg)}, excluded_trialidxs_from_residuals{unique_configs(icfg)});
    
    if(~strcmp(pars.alignment.type, 'none'))
        % Only merging data across sessions if alignment was done
        [aligned_single_trials{unique_configs(icfg)}] = mergeDataAcrossSessions(aligned_single_trials{unique_configs(icfg)});
        [aligned_residuals{unique_configs(icfg)}] = mergeDataAcrossSessions(aligned_residuals{unique_configs(icfg)});
        
    end
end
    
analysis.aligned_single_trials = aligned_single_trials;
analysis.residuals = aligned_residuals;
analysis.expt_names = expt_names_in_config;
analysis.alignment_info.stats = align_stats;
analysis.alignment_info.projs = align_projs;
analysis.alignment_info.task_conds_proj = task_conds_proj;
analysis.dspInfo = dsp_info;
analysis.params = pars;
analysis.metadata.project = data_dir;
analysis.metadata.animal = animal;

end

function [D_align] = compute_prealigned(data, U_orth, mu)

num_sessions = length(data);
D_align = cell(1, num_sessions);
            
for i_session = 1:num_sessions
    
    D_align{i_session} = data{i_session};
    [~,num_times,num_trials] = size(data{i_session}.response);
    meanSub_response_full = (data{i_session}.response(:,:) - mu{i_session});
    aligned_response = U_orth{i_session}' * meanSub_response_full;
    D_align{i_session}.response = reshape(aligned_response,[size(U_orth{i_session},2) num_times num_trials]);
    D_align{i_session}.U = U_orth{i_session};
    D_align{i_session}.mu_subtracted = mu{i_session}';

    
end

end

function [single_trials] = filter_excluded_trials(single_trials, excluded_idxs)

if(~isempty(excluded_idxs))
    
    num_sessions = length(excluded_idxs);
    for i_session = 1:num_sessions
       getpars.trial_get = ~excluded_idxs{i_session};
       data_ = tdrGetTrials(single_trials{i_session},getpars);
       excluded_fields = setdiff(fieldnames(single_trials{i_session}), fieldnames(data_));

       if(~isempty(excluded_fields))

           for iexf = 1:length(excluded_fields)
               data_.(excluded_fields{iexf}) = single_trials{i_session}.(excluded_fields{iexf});
           end

           single_trials{i_session} = data_;
       end
    end
else
    return
end
end


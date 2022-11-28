%{These set of functions extract the raw spike counts for each session and 
% performs a few important preprocessing steps on the raw spike counts to 
% to make it amenable for the analyses in the residual dynamics pipeline.
% In particular, these set of function perform the following operations.
% 1) Removes any specified sessions that we do not want to analyze
% 2) Assign "task configurations" to each recorded session based on the
%    approx. locations of the saccade targets
% 3) Divides a session into "smaller experiments" to ensure that units are
%    "stable" within an experiment. A change-point analysis is done to 
%    obtain the segmentation fo individual sessions into experiments
% 4) Removes "bad units" - bad units are either those that do not fire at
%    all (i.e potentially indicating faulty spike-sorting) or neurons that
%    show large non-stationarities in the temporally-averaged firing rates 
%    within an experiment. 
% 5) Finally, computes the single-trial trajectories binned at a specific
%    resolution for each experiment
% Author : Aniruddh Galgali (Jan 2018)
%}
function extract_and_preprocess_data(data_dir, animal, save_project, pars)
%{ Top-level function to preprocess data
%  INPUT
%  - data_dir (struct) - contains info about path to the directory 
%    containing the raw neural data
%  - animal (string) - animal name to analyze
%  - save_project (struct) - info about where to store
%    preprocessed data
%  - pars(struct) - preprocessing parameters which determine which units 
%    get romoved, how to determine the segmentation and what bin-size to use 
%}
[valid_session_list, excluded_session_list] = get_session_list(data_dir, animal, pars.session.exclude_list);
pars.session.exclude_list = excluded_session_list;
for isession = 1: length(valid_session_list)
    
    expt = AnalysisLoad(data_dir.name, valid_session_list{isession},'fast');
    [expt.targInfo] = get_saccadetarget_configurations(expt, pars.session.target_discrepancy_threshold);   
    
    if(isnan(expt.targInfo.config))
        pars.session.sessions_to_exclude = cat(2,pars.session.exclude_list,valid_session_list{isession}); 
    	continue
    else
        expt = get_session_data(expt, pars);
        AnalysisSave(save_project.loader.name, valid_session_list{isession}, expt);
        fprintf('Done with %s\n',valid_session_list{isession})
    end
end
save(fullfile(save_project.root_dir, save_project.name,'metadata',strcat(animal, '_parameters.mat')),'pars');

end

function [expt] = get_session_data(expt, pars)
%{ Helper function that calls all the subroutines
%  INPUT
%  - expt (struct) - contains the data and associated info.
%  - pars(struct) - preprocessing parameters
%  OUTPUT
%  - expt (struct) - modified data structure
%}
    [expt] = extract_alignments(expt, pars.align);
    
    [expt.unit_exclusion_info, expt.session_segmentation_info] = ...
        extract_segments_and_bad_units(expt, pars.session.segmentation, pars.unit_exclusion);

    [expt] = extract_trajectories(expt, pars.single_trial);
end


function [expt] = extract_alignments(expt, align_pars)
%{ Retains only specified epoch alignments
%  INPUT
%  - expt (struct) - contains the data and associated info
%  - align_pars(struct) - specifies which epoch alignments to retain
%  OUTPUT
%  - expt (struct) - modified data structure
%}
    all_alignments = fieldnames(expt.align);
    aligns_to_keep = align_pars.aligns_to_keep;
    aligns_to_remove = all_alignments(~ismember(all_alignments, aligns_to_keep));
    for ialign = 1: length(aligns_to_remove)
        expt.align = rmfield(expt.align,aligns_to_remove{ialign});
    end

    for ialign = 1:length(aligns_to_keep)

        analysisTimeWindow.tStart = align_pars.pre_analysis_window.tLims{ialign}(1);
        analysisTimeWindow.tEnd = align_pars.pre_analysis_window.tLims{ialign}(end); 
        [D] = extractAnalysisTimeWindow(expt.align.(aligns_to_keep{ialign}).data, analysisTimeWindow);
        [D] = removeFsMultiplier(D);      
        expt.align.(aligns_to_keep{ialign}).data = D;
   
    end
            
end

function [unit_exclusion_info, session_segmentation_info] = extract_segments_and_bad_units(expt, session_segmentation_pars, unit_exclusion_pars)
%{ Identifies the bad units and performs change-point analyses to determine 
%  how to split a session into "stable" segments.
%  INPUT
%  - expt (struct) - contains the data and associated info.
%  - session_segmentation_pars(struct) - specifies information
%    relevant to the change-point algorithm used for identifying stable
%    segments
%  - unit_exclusion_pars(struct) - specifies criteria for unit
%    removal
%  OUTPUT
%  - unit_exclusion_info (struct) - information about which units get
%    removed (mostly indices)
%  - session_segmentation_info (struct) - information about session
%    segments (indices indicating where to split the session)
%}
    all_alignments = fieldnames(expt.align);

    for ialign = 1:length(all_alignments)
        D(ialign) = expt.align.(all_alignments{ialign}).data;
    end

    [Dmerged] = mergeDataAlignments(D);

    [idxs_silent_units, tempaveraged_FR_per_trial] = excludeSilentUnits(Dmerged, unit_exclusion_pars.quiescentUnits.threshold);

    [idxs_lowchoiceselective] = computeChoiceSelectivity(Dmerged,unit_exclusion_pars.choiceSel.threshold,...
    unit_exclusion_pars.choiceSel.numTimeSegments);

    if(unit_exclusion_pars.do_remove_quiescent && unit_exclusion_pars.do_remove_lowchoiceselective)
        idxs_remove = idxs_silent_units & idxs_lowchoiceselective;
    elseif(unit_exclusion_pars.do_remove_quiescent && ~unit_exclusion_pars.do_remove_lowchoiceselective)
        idxs_remove = idxs_silent_units;
    elseif(~unit_exclusion_pars.do_remove_quiescent && unit_exclusion_pars.do_remove_lowchoiceselective)
        idxs_remove = idxs_lowchoiceselective;
    else
        idxs_remove = false(length(idxs_silent_units),1);
    end

    unit_exclusion_info.silent_units = idxs_silent_units;
    unit_exclusion_info.low_choiceselectivity_units = idxs_lowchoiceselective;
    unit_exclusion_info.remove_pre_final = idxs_remove;
    tempaveraged_FR_per_trial(idxs_remove,:) = [];

    %% Finding change points for each session

    session_segmentation_info = extract_session_segments(tempaveraged_FR_per_trial, session_segmentation_pars);


    %% Finding the non-stationary units for each segment

    num_segments = length(session_segmentation_info.idxsSegment);
    segments_to_include = true(length(session_segmentation_info.idxsSegment),1);
    nonStat_out = cell(num_segments,1);
    for ii = 1:num_segments

        idxsSegment = session_segmentation_info.idxsSegment{ii}; 
        segment_response = tempaveraged_FR_per_trial(:,idxsSegment);  
        [nonStat_out{ii},idxsNonStat] = removeNonStationaryUnits(segment_response,unit_exclusion_pars.nonStat);

        if(isempty(nonStat_out{ii}))
          segments_to_include(ii) = false;
          unit_exclusion_info.non_stationary_units{ii} = [];
        else
          unit_exclusion_info.non_stationary_units{ii} = idxsNonStat;
        end
      
    end
    session_segmentation_info.segment_mask = segments_to_include;
    session_segmentation_info.tempaveraged_FR_per_trial = tempaveraged_FR_per_trial;
end


function [expt] = extract_trajectories(expt, trajectory_pars)
%{ Extracts the single-trial trajectories for each experiments after applying 
%  all the pre-processing criteria
%  INPUT
%  - expt (struct) - contains the data and associated info.
%  - trajectory_pars(struct) - specifies information about timw-window,
%    bin-size and other properites determining the single trials
%  OUTPUT
%  - expt (struct) - modified data structure
%}
    all_alignments = fieldnames(expt.align);
    
    if any(ismember(all_alignments,{'outfixcd'}))
        [binned_delay_lengths,idx_delay_lengths] = binDelayLengths(expt.align.outfixcd.data, trajectory_pars.delay_length_bin_boundaries);
        [reaction_times] = abs(expt.align.outfixcd.data.task_event.fpoffcd);
        [~,~,idx_RT_bins] = histcounts(reaction_times,[0 prctile(reaction_times,[25 50 75]) max(reaction_times)]);
    end
    unique_binned_delay_lengths = unique(binned_delay_lengths);
    
    for ialign = 1:length(all_alignments)
        
        D = expt.align.(all_alignments{ialign}).data;
        analysisTimeWindow.tStart = trajectory_pars.analysis_window.tLims{ialign}(1);
        analysisTimeWindow.tEnd = trajectory_pars.analysis_window.tLims{ialign}(end); 
        [D] = extractAnalysisTimeWindow(D, analysisTimeWindow);
        current_sampling_freq = (1./(D.time(2) - D.time(1)));
        
        [D] = binSpikeCounts(D, trajectory_pars.bin_size, current_sampling_freq, ...
            trajectory_pars.use_sqrt_transform, trajectory_pars.use_rates);
          
        D.task_variable.delay_length = binned_delay_lengths;
        D.task_index.delay_length = idx_delay_lengths;
        nvl = length(unique_binned_delay_lengths);
        D.task_variable_index.delay_length = [unique_binned_delay_lengths (1:nvl)'];
    
        D.task_variable.reaction_time = reaction_times;
        D.task_index.reaction_time = idx_RT_bins;
        
        D.task_variable.target_config = expt.targInfo.config.*ones(size(D.response,3),1);
        D.task_index.target_config = expt.targInfo.config.*ones(size(D.response,3),1);
        
        D_out(ialign) = D;
   
    end
    
    expt = rmfield(expt,'align');
    
    [Dmerged] = mergeDataAlignments(D_out);
    [Dmerged] = applyUnitRemoval(Dmerged, expt.unit_exclusion_info.remove_pre_final);
    valid_array_locs = expt.dspInfo(~expt.unit_exclusion_info.remove_pre_final); 

    segment_indices = expt.session_segmentation_info.idxsSegment;
    segment_masks = expt.session_segmentation_info.segment_mask;
    idxsNonStatUnits = expt.unit_exclusion_info.non_stationary_units;
    
    expt.data = {};
    expt.dspInfo = {};
    c_count = 1;
    for ii = 1:length(segment_indices)

        if(segment_masks(ii))
            
            getpars.trial_get = segment_indices{ii};
            expt.data{c_count} = tdrGetTrials(Dmerged, getpars);
            [expt.data{c_count}] = applyUnitRemoval(expt.data{c_count},idxsNonStatUnits{ii});
            expt.dspInfo{c_count} = valid_array_locs(~idxsNonStatUnits{ii});
            c_count = c_count + 1;
        else
            continue;
        end
    end
 
end
        

